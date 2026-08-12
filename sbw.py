"""Shared SBW pipeline: instrument constants, sample generation, ConvNext model.

Data-generation recipe distilled from SBW.ipynb: build a scene image
(exozodi ring + planet pixel; star is all zeros), then project it through the
fixed single-Bracewell transmission-map scenery to get the (wavelength, time)
signal `signal_both_sbw`. The scenery depends only on instrument constants, so
we precompute it once and every sample is a single matmul.
"""
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from transmission_maps import point_signal_SBW

# --- Instrument constants (from SBW.ipynb) ---
baseline = 10.0          # [m]
q = 6                    # imaging/nulling baseline ratio (unused for SBW)
lambda_ref = 10e-6       # [m]
fov_as = 0.4             # [arcsec]
Npix = 80
alpha = np.linspace(-fov_as, fov_as, Npix)
beta = np.linspace(-fov_as, fov_as, Npix)
ALPHA, BETA = np.meshgrid(alpha, beta)

time_total = 30 * 3600   # [s]
read_time = 1000         # [s]
N = int(time_total / read_time)                     # 108 timesteps
t = np.linspace(0, time_total, N)
Omega = (2 / time_total) * 2 * np.pi                # [rad/s]

R = 20
lambda_min, lambda_max = 4e-6, 18.5e-6
lambdas = np.logspace(np.log10(lambda_min), np.log10(lambda_max),
                      int(R * np.log(lambda_max / lambda_min) + 0.5), base=10)
NL = len(lambdas)                                   # 31 wavelengths

sig2 = 0.1 ** 2          # exozodi ring width^2

# --- Parameters: name, low, high, is_int ---
PARAMS = [
    ("idx1", 0, Npix, True),
    # idx2 restricted to one half-plane: SBW can't distinguish a planet from
    # its 180-deg reflection through the star, so we keep one representative.
    ("idx2", 0, Npix // 2, True),
    ("Amplitude_planet", 100.0, 5000.0, False),
    ("Amplitude_exozodi", 10.0, 100.0, False),
    ("inc", 0.0, np.pi / 2, False),
    ("phi", 0.0, np.pi, False),
]
NAMES = [p[0] for p in PARAMS]
LOWS = np.array([p[1] for p in PARAMS], dtype=np.float64)
HIGHS = np.array([p[2] for p in PARAMS], dtype=np.float64)


def sample_params(rng):
    """One random parameter vector [idx1, idx2, amp_p, amp_e, inc, phi]."""
    return np.array([
        rng.integers(0, Npix),
        rng.integers(0, Npix // 2),   # half-plane: removes the SBW reflection degeneracy
        rng.uniform(100.0, 5000.0),
        rng.uniform(10.0, 100.0),
        rng.uniform(0.0, np.pi / 2),
        rng.uniform(0.0, np.pi),
    ], dtype=np.float64)


def build_scene(params):
    """Scene image (Npix, Npix) = exozodi ring + planet pixel (star is zero)."""
    idx1, idx2, amp_p, amp_e, inc, phi = params
    alpha_rot = ALPHA * np.cos(phi) + BETA * np.sin(phi)
    beta_rot = 1 / np.cos(inc) * (-ALPHA * np.sin(phi) + BETA * np.cos(phi))
    r_rot2 = alpha_rot ** 2 + beta_rot ** 2
    im = amp_e * (r_rot2 / sig2) ** 2 * np.exp(-r_rot2 / sig2)
    im[int(idx1), int(idx2)] += amp_p
    return im


def transmission_scenery():
    """Precompute fixed SBW scenery, flattened to (NL*N, Npix*Npix)."""
    T = point_signal_SBW(t, ALPHA, BETA, baseline, lamda_ref=lambdas, Omega=Omega)
    return T.reshape(NL * N, Npix * Npix)


def signal_from_scene(im, T_flat):
    """Project a scene through the scenery -> signal image (NL, N)."""
    return (T_flat @ im.reshape(-1)).reshape(NL, N)


# --- Normalization ---
def norm_targets(y):
    return (y - LOWS) / (HIGHS - LOWS)


def denorm_targets(yn):
    return yn * (HIGHS - LOWS) + LOWS


# --- Model ---
PATCH = 4
H_PAD = NL + (-NL % PATCH)   # pad wavelength axis 31 -> 32 (divisible by PATCH)


class ViT(nn.Module):
    """Small Vision Transformer over the native (NL, N) signal image."""

    def __init__(self, dim=192, depth=6, heads=6, n_out=len(PARAMS)):
        super().__init__()
        self.proj = nn.Conv2d(IN_CH, dim, kernel_size=PATCH, stride=PATCH)  # patch embed
        self.gh, self.gw = H_PAD // PATCH, N // PATCH
        ntok = self.gh * self.gw
        self.cls = nn.Parameter(torch.zeros(1, 1, dim))
        self.pos = nn.Parameter(torch.zeros(1, ntok + 1, dim))
        nn.init.trunc_normal_(self.cls, std=0.02)
        nn.init.trunc_normal_(self.pos, std=0.02)
        layer = nn.TransformerEncoderLayer(dim, heads, dim * 4, batch_first=True,
                                           activation="gelu", norm_first=True)
        self.enc = nn.TransformerEncoder(layer, depth)
        self.norm = nn.LayerNorm(dim)
        self.head = nn.Linear(dim, n_out)

    def forward(self, x):                          # x: (B, 1, NL, N)
        x = F.pad(x, (0, 0, 0, H_PAD - NL))        # pad wavelength axis
        x = self.proj(x).flatten(2).transpose(1, 2)  # (B, ntok, dim)
        cls = self.cls.expand(x.size(0), -1, -1)
        x = torch.cat([cls, x], 1) + self.pos
        x = self.enc(x)
        return self.head(self.norm(x[:, 0]))


def build_model():
    return ViT()


IN_CH = 2  # channel 0: raw (absolute scale); channel 1: high-pass (structure)


def prep_image(sig, mean, std):
    """(NL, N) signal -> 2-channel (IN_CH, NL, N) tensor.

    ch0 keeps absolute brightness (-> Amplitude_exozodi); ch1 removes the large
    per-image pedestal so the faint planet ripple isn't swamped. Both use only
    the observable signal itself.
    """
    sig = np.asarray(sig, dtype=np.float64)
    ch0 = (sig - mean) / std
    ch1 = (sig - sig.mean()) / std
    return torch.as_tensor(np.stack([ch0, ch1]), dtype=torch.float32)


def demo():
    """Self-check: shapes and linearity of the projection."""
    T_flat = transmission_scenery()
    assert T_flat.shape == (NL * N, Npix * Npix), T_flat.shape
    rng = np.random.default_rng(0)
    p = sample_params(rng)
    sig = signal_from_scene(build_scene(p), T_flat)
    assert sig.shape == (NL, N), sig.shape
    # planet contributes linearly on top of exozodi
    p2 = p.copy(); p2[2] = 0.0  # amp_planet = 0
    sig_ez = signal_from_scene(build_scene(p2), T_flat)
    planet_only = T_flat[:, int(p[0]) * Npix + int(p[1])] * p[2]
    assert np.allclose(sig - sig_ez, planet_only.reshape(NL, N)), "linearity broken"
    # normalization round-trips
    assert np.allclose(denorm_targets(norm_targets(p)), p)
    x = prep_image(sig, sig.mean(), sig.std())
    assert x.shape == (IN_CH, NL, N), x.shape
    assert build_model()(x[None]).shape == (1, len(PARAMS))
    print("sbw.demo OK  sig", sig.shape, " T", T_flat.shape)


if __name__ == "__main__":
    demo()
