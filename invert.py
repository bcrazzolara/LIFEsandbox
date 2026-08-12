"""Analytic-ish inversion: recover the 6 SBW params from a signal image.

The forward map is *linear in the scene*:
    sig = amp_e * (T @ ring(inc, phi)) + amp_p * T[:, pixel]
so only (inc, phi) are nonlinear. Given them, the ring signal is fixed and the
problem is: find the single planet column T[:,k] and the two amplitudes that
best explain `sig`. That inner problem is an exact 2-basis least squares per
candidate pixel, solved in closed form and vectorised over all pixels. The
outer 2-D (inc, phi) search minimises the residual; at the true angles it hits
~0, so all six params come back to ~machine precision on noiseless data.

Known unrecoverables (structural, not solver faults):
  - pixel & its 180deg reflection share the same T column -> we search only the
    idx2<Npix/2 half-plane (the sampling convention).
  - at inc=0 the ring is a circle -> phi carries no information.
"""
import numpy as np
from scipy.optimize import minimize

from sbw import (transmission_scenery, sample_params, build_scene,
                 signal_from_scene, ALPHA, BETA, sig2, Npix, NAMES)

T = transmission_scenery()                 # (M, K) = (NL*N, Npix*Npix)
BB = (T * T).sum(0)                         # ||T[:,k]||^2, precomputed
# half-plane pixel mask matching the sampling convention (idx2 in [0, Npix/2))
_j = np.arange(Npix * Npix) % Npix
HALF = np.where(_j < Npix // 2)[0]


def ring_vec(inc, phi):
    ar = ALPHA * np.cos(phi) + BETA * np.sin(phi)
    br = (-ALPHA * np.sin(phi) + BETA * np.cos(phi)) / np.cos(inc)
    r2 = ar ** 2 + br ** 2
    return ((r2 / sig2) ** 2 * np.exp(-r2 / sig2)).reshape(-1)


def _best_pixel(inc, phi, y, yy):
    """Given angles, closed-form fit of (pixel, amp_e, amp_p); return residual."""
    a = T @ ring_vec(inc, phi)             # ring signal (M,)
    aa = a @ a
    ay = a @ y
    ab = a @ T                             # (K,)  a . T[:,k]
    by = y @ T                             # (K,)  y . T[:,k]
    det = aa * BB - ab ** 2                # 2x2 normal-eq determinant per pixel
    det = np.where(np.abs(det) < 1e-30, np.inf, det)
    ce = (BB * ay - ab * by) / det         # amp_e per pixel
    cp = (aa * by - ab * ay) / det         # amp_p per pixel
    resid = yy - (ce * ay + cp * by)       # ||y - proj||^2 per pixel
    resid[cp < 0] = np.inf                 # planet amplitude must be >= 0
    k = HALF[np.argmin(resid[HALF])]
    return resid[k], k, ce[k], cp[k]


def invert(sig):
    y = np.asarray(sig, dtype=np.float64).reshape(-1)
    yy = y @ y
    # coarse grid over (inc, phi) to seed the local search
    incs = np.linspace(0.02, np.pi / 2 - 0.02, 12)
    phis = np.linspace(0.0, np.pi, 12, endpoint=False)
    best = (np.inf, 0.0, 0.0)
    for inc in incs:
        for phi in phis:
            r = _best_pixel(inc, phi, y, yy)[0]
            if r < best[0]:
                best = (r, inc, phi)
    # refine angles (Nelder-Mead: inner argmin makes the objective non-smooth)
    res = minimize(lambda a: _best_pixel(a[0], a[1], y, yy)[0],
                   x0=[best[1], best[2]], method="Nelder-Mead",
                   options=dict(xatol=1e-6, fatol=1e-14, maxiter=2000))
    inc, phi = res.x
    inc = np.clip(inc, 0.0, np.pi / 2 - 1e-6)
    phi = phi % np.pi
    _, k, amp_e, amp_p = _best_pixel(inc, phi, y, yy)
    return np.array([k // Npix, k % Npix, amp_p, amp_e, inc, phi])


def demo(n=20, seed=0):
    rng = np.random.default_rng(seed)
    T_flat = T
    errs = []
    for _ in range(n):
        p = sample_params(rng)
        sig = signal_from_scene(build_scene(p), T_flat)
        est = invert(sig)  # noiseless
        # idx must be exact; amps to ~0.1%; angles skipped when inc~0 (phi free)
        assert est[0] == int(p[0]) and est[1] == int(p[1]), (p[:2], est[:2])
        assert abs(est[2] - p[2]) < 1.0, ("amp_p", p[2], est[2])
        assert abs(est[3] - p[3]) < 0.1, ("amp_e", p[3], est[3])
        assert abs(est[4] - p[4]) < 1e-3, ("inc", p[4], est[4])
        if p[4] > 0.1:                     # phi only meaningful away from inc=0
            dphi = abs((est[5] - p[5] + np.pi / 2) % np.pi - np.pi / 2)
            assert dphi < 1e-2, ("phi", p[5], est[5])
        errs.append(np.abs(est - p))
    print("invert.demo OK over", n, "samples")
    print("max |err| per param:",
          dict(zip(NAMES, np.max(errs, 0).round(6))))


def noise_test(n=30, seed=1, fracs=(0.0, 0.01, 0.05, 0.1, 0.2)):
    """Median per-param |err| vs Gaussian noise sigma = frac * sig.std()."""
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(n):
        p = sample_params(rng)
        samples.append((p, signal_from_scene(build_scene(p), T)))
    print(f"\nrobustness ({n} samples), median |err| per param:")
    print("  noise%  " + "  ".join(f"{x:>10}" for x in NAMES) + "   idx_hit%")
    for f in fracs:
        errs, hits = [], 0
        for p, sig in samples:
            noisy = sig + rng.normal(0, f * sig.std(), sig.shape)
            est = invert(noisy)
            hits += (est[0] == int(p[0]) and est[1] == int(p[1]))
            # angle error mod pi for phi
            e = np.abs(est - p)
            e[5] = abs((est[5] - p[5] + np.pi / 2) % np.pi - np.pi / 2)
            errs.append(e)
        med = np.median(errs, 0)
        print(f"  {f * 100:5.0f}%  " + "  ".join(f"{v:>10.4g}" for v in med)
              + f"   {100 * hits / n:5.0f}")


if __name__ == "__main__":
    demo()
    noise_test()
