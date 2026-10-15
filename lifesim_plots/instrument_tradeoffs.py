import contextlib
import gc
import importlib.util
import io
import json
import time
import warnings
from pathlib import Path
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pandas as pd
import torch
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm

from phringe.core.instrument import Instrument
from phringe.core.observation import Observation
from phringe.core.scene import Scene
from phringe.core.sources.local_zodi import LocalZodi
from phringe.lib.array_configuration import SBWConfiguration
from phringe.lib.beam_combiner import SingleBracewell
from phringe.main import PHRINGE
import phringe.io.sed_loader as sed_loader_module
from phringe.io.sed_loader import SEDLoader

from lifesimmc.core.pipeline import Pipeline
from lifesimmc.core.modules.base_module import BaseModule
from lifesimmc.core.modules.loading.setup_module import SetupModule
from lifesimmc.core.modules.generating.data_generation_module import DataGenerationModule
from lifesimmc.core.modules.generating.template_generation_module import TemplateGenerationModule
from lifesimmc.core.modules.processing.zca_whitening_module import ZCAWhiteningModule
from lifesimmc.core.modules.processing.ml_parameter_estimation_module import MLSEDEstimationModule
from lifesimmc.core.resources.data_resource import DataResource
from lifesimmc.core.resources.planet_resource import PlanetResource
from lifesimmc.core.resources.resource_collection import ResourceCollection
from lifesimmc.core.resources.transformation_resource import TransformationResource

import torch

RAD_TO_MAS = 180 / np.pi * 3600e3

def to_numpy(x):
    if torch.is_tensor(x):
        return x.detach().cpu().numpy()
    return np.asarray(x)

def to_float_scalar(x):
    if torch.is_tensor(x):
        return float(x.detach().cpu().item())
    return float(x)

def quiet(func, *args, **kwargs):
    # The lifesimmc modules print progress messages; keep the notebook output readable
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('ignore')
        return func(*args, **kwargs)



# --- Observation ---
T_OBS = 5 * 3600            # total integration time [s]
T_MOD = 10 * 3600           # modulation period [s] (one full array rotation); T_OBS / T_MOD = half a rotation
DIT = 2 * 60                # detector integration time [s] -> T_OBS / DIT = 60 time steps, as in the prototype

# --- Monte Carlo ---
N_REALISATIONS = 10          # noise realisations per (planet, setup) (~45 s each on CPU, much less on GPU)
BASE_SEED = 1000            # realisation i of every (planet, setup) uses seed BASE_SEED + i
FORCE_RECOMPUTE = False      # True: delete cached realisations and recompute all; False: reuse finished runs
MAX_POS_OFFSET_SIGMA = 5   # free fits further than this from +/- the true position count as failed

# --- Simulation ---
GRID_SIZE = 50              # spatial grid of the extended sources and the position-search templates
GPU_INDEX = 2            # e.g. 2 on the GPU machine; None -> CPU

# --- Output ---
# Per planet, next to its spectrum (as in the prototype): <spectrum folder>/setup_snr/realisations_<planet>.jsonl
# (one line per (setup, realisation); finished runs are reused) and the SNR-per-setup plot.
# All planets together (summary table, heatmap):
COMBINED_DIR = Path('plots/setup_snr')
COMBINED_DIR.mkdir(parents=True, exist_ok=True)




# =============== DEFINITION OF TARGETS =================

# planet name -> config file (relative to lifesim_plots), planet variable in it, optional replacement spectrum
TARGETS = {
    
    'PDS70 b': dict(config='setup/spectra/PDS70/pds_70.py', planet='pds_70_b'),
    'PDS70 c': dict(config='setup/spectra/PDS70/pds_70.py', planet='pds_70_c'),
    'HR 8799 d': dict(config='setup/spectra/HR8799/hr_8799.py', planet='hr_8799_d'),
    '51 Eri b': dict(config='setup/spectra/51Eri/_51_eri.py', planet='_51_eri_b'),
    'HR 8799 e': dict(config='setup/spectra/HR8799/hr_8799.py', planet='hr_8799_e'),
    'bet Pic b': dict(config='setup/spectra/BetaPic/bet_pic.py', planet='bet_pic_b'),
    'HR 8799 c': dict(config='setup/spectra/HR8799/hr_8799.py', planet='hr_8799_c'),
    'Kap And b': dict(config='setup/spectra/KapAnd/kap_and.py', planet='kap_and_b'),
    'HR 8799 b': dict(config='setup/spectra/HR8799/hr_8799.py', planet='hr_8799_b'),
    'DH Tau b': dict(config='setup/spectra/DHTau/dh_tau.py', planet='dh_tau_b'),
    'HD 106906 b': dict(config='setup/spectra/HD106906/hd_106906.py', planet='hd_106906_b'),
}

# =============== INSTRUMENT =================
# name -> D [m], B [m], T_mirror / T_instr / T_det [K], null depth, dark current [e-/s/pixel]




NULL_DEPTH = 1e-4
'''
SETUPS = {
    'Large cold':      dict(D=0.20,  B=2.0,  T_mirror=90, T_instr=80,  T_det=40, null_depth=NULL_DEPTH, dark_current=0.035),
    'Large warm':      dict(D=0.20,  B=2.0,  T_mirror=140, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=0.35),

    'Medium cold':     dict(D=0.15,  B=1.5,  T_mirror=90, T_instr=80,  T_det=40, null_depth=NULL_DEPTH, dark_current=0.035),
    'Medium warm':     dict(D=0.15,  B=1.5,  T_mirror=140, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=0.35),

    'Prototype cold':  dict(D=0.10,  B=1.0,  T_mirror=90, T_instr=80,  T_det=40, null_depth=NULL_DEPTH, dark_current=0.035),
    'Prototype warm':  dict(D=0.10,  B=1.0,  T_mirror=140, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=0.35),
    #'Prototype warm 5':  dict(D=0.10,  B=1.0,  T_mirror=145, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=1.0),

    'Small cold':      dict(D=0.075, B=0.75, T_mirror=90, T_instr=80,  T_det=40, null_depth=NULL_DEPTH, dark_current=0.035),
    'Small warm':      dict(D=0.075, B=0.75, T_mirror=140, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=0.35),
    #'Small warm 5':      dict(D=0.075, B=0.75, T_mirror=190, T_instr=150, T_det=80, null_depth=NULL_DEPTH, dark_current=1.0),

    'Short cold':      dict(D=0.05,  B=0.5,  T_mirror=90, T_instr=80,  T_det=40, null_depth=NULL_DEPTH, dark_current=0.035),
    'Short warm':      dict(D=0.05,  B=0.5,  T_mirror=140, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=0.35),
    #'Short hot 5':      dict(D=0.05,  B=0.5,  T_mirror=190, T_instr=150, T_det=80, null_depth=NULL_DEPTH, dark_current=1.0),
}
'''
SETUPS = {
    'Large cold':      dict(D=0.20,  B=2.0,  T_mirror=90, T_instr=80,  T_det=40, null_depth=NULL_DEPTH, dark_current=0.035),
    'Large warm':      dict(D=0.20,  B=2.0,  T_mirror=140, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=0.35),

    'Medium cold':     dict(D=0.15,  B=1.5,  T_mirror=90, T_instr=80,  T_det=40, null_depth=NULL_DEPTH, dark_current=0.035),
    'Medium warm':     dict(D=0.15,  B=1.5,  T_mirror=140, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=0.35),

    'Prototype cold':  dict(D=0.1,  B=1.0,  T_mirror=90, T_instr=80,  T_det=40, null_depth=NULL_DEPTH, dark_current=0.035),
    'Prototype warm':  dict(D=0.1,  B=1.0,  T_mirror=140, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=0.35),
    #'Prototype warm 5':  dict(D=0.10,  B=1.0,  T_mirror=145, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=1.0),

    'Small cold':      dict(D=0.075, B=0.75, T_mirror=90, T_instr=80,  T_det=40, null_depth=NULL_DEPTH, dark_current=0.035),
    'Small warm':      dict(D=0.075, B=0.75, T_mirror=140, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=0.35),
    #'Small warm 5':      dict(D=0.075, B=0.75, T_mirror=190, T_instr=150, T_det=80, null_depth=NULL_DEPTH, dark_current=1.0),

    'Short cold':      dict(D=0.05,  B=0.5,  T_mirror=90, T_instr=80,  T_det=40, null_depth=NULL_DEPTH, dark_current=0.035),
    'Short warm':      dict(D=0.05,  B=0.5,  T_mirror=140, T_instr=120, T_det=80, null_depth=NULL_DEPTH, dark_current=0.35),
    #'Short hot 5':      dict(D=0.05,  B=0.5,  T_mirror=190, T_instr=150, T_det=80, null_depth=NULL_DEPTH, dark_current=1.0),
}

def load_target(name):
    # Executes the config and returns its star, exozodi, disks and the planet (with `sed_file` if given)
    target = TARGETS[name]
    path = Path(target['config'])
    spec = importlib.util.spec_from_file_location(f'_target_{path.stem}', path)
    module = importlib.util.module_from_spec(spec)
    # A config may point at a missing spectrum (PLACEHOLDER.txt): load it without, then use `sed_file`
    original = sed_loader_module.SEDLoader
    sed_loader_module.SEDLoader = lambda **kw: original(**kw) if Path(kw['path_to_file']).exists() else None
    try:
        quiet(spec.loader.exec_module, module)
    finally:
        sed_loader_module.SEDLoader = original

    planet = getattr(module, target['planet'])
    if 'sed_file' in target:
        planet = planet.model_copy(update={'sed_loader': SEDLoader(
            path_to_file=target['sed_file'], sed_units=target['sed_units'], wavelength_units='um',
            observed_planet_radius=planet.radius, observed_host_star_distance=module.star.distance)})
    if planet.sed_loader is None:
        raise FileNotFoundError(f'{name}: the config has no spectrum file; set sed_file in TARGETS')
    out_dir = Path(planet.sed_loader.path_to_file).parent / 'setup_snr'
    out_dir.mkdir(exist_ok=True)
    slug = name.lower().replace(' ', '_')
    return dict(star=module.star, exozodi=getattr(module, 'exozodi', None), disks=list(getattr(module, 'disks', [])),
                planet=planet, results_file=out_dir / f'realisations_{slug}.jsonl',
                plot_file=out_dir / f'snr_per_setup_{slug}.png')


LOADED = {name: load_target(name) for name in TARGETS}
pd.DataFrame([dict(planet=name, star=t['star'].name, distance_pc=round(t['star'].distance / 3.0857e16, 1),
                   n_disks=len(t['disks']), spectrum=Path(t['planet'].sed_loader.path_to_file).name)
              for name, t in LOADED.items()]).set_index('planet')

# loading setups
pd.DataFrame(SETUPS).T.rename(columns={'D': 'D [m]', 'B': 'B [m]', 'T_mirror': 'T_mirror [K]',
                                       'T_instr': 'T_instr [K]', 'T_det': 'T_det [K]',
                                       'dark_current': 'dark current [e-/s/px]'})



def make_instrument(p):
    return Instrument(
        array_configuration_matrix=SBWConfiguration.acm,
        complex_amplitude_transfer_matrix=SingleBracewell.catm,
        kernels=SingleBracewell.kernels,
        nulling_baseline_min=0 * u.m,
        nulling_baseline_max=5 * u.m,
        throughput=0.12,
        wavelength_min=4 * u.um,
        wavelength_max=5 * u.um,
        wavelength_bands_boundaries=[],
        spectral_resolving_power=100,
        spectral_resolving_power_instrument=100,
        amplitude_perturbation=None,
        phase_perturbation=None,
        polarization_perturbation=None,
        null_depth=p['null_depth'],
        aperture_diameter=p['D'] * u.m,
        mirror_temperature=p['T_mirror'] * u.K,
        mirror_emissivity=0.06,
        mirror_throughput=0.94,
        instrument_temperature=p['T_instr'] * u.K,
        instrument_emissivity=1.0,
        quantum_efficiency=0.7,
        detector_temperature=p['T_det'] * u.K,
        detector_wavelength_min=0.6 * u.um,
        detector_wavelength_max=5.3 * u.um,
        pixel_size=18e-6 * u.m,
        pixels_per_wavelength_bin=2.2,
        dark_current=p['dark_current'],
    )

def make_observation(p):
    return Observation(
        solar_ecliptic_latitude=0 * u.deg,
        total_integration_time=T_OBS * u.s,
        detector_integration_time=DIT * u.s,
        modulation_period=T_MOD * u.s,
        nulling_baseline=p['B'] * u.m,
    )

def make_scene(target, with_planet=True):
    # Fresh copies of all sources: a source keeps a reference to the simulation it was last added to
    sources = [target['star'], target['exozodi'], *target['disks']] + ([target['planet']] if with_planet else [])
    scene = Scene()
    scene.add_source(LocalZodi())
    for source in sources:
        if source is not None:
            scene.add_source(source.model_copy())
    return scene



class PoissonWhiteningModule(BaseModule):
    # Mean subtraction per channel, then division by the per-sample Poisson noise sigma = sqrt(b_lambda + m_lambda,t),
    # with m the planet model counts (kernel x wavelength x time) and b the background level estimated from the data.

    def __init__(self, n_data_in, n_data_out, n_transformation_out, planet_counts):
        super().__init__()
        self.n_data_in = n_data_in
        self.n_data_out = n_data_out
        self.n_transformation_out = n_transformation_out
        self.planet_counts = np.clip(np.asarray(planet_counts, dtype=float), 0, None)

    def run(self, pipeline_resources):
        data = self.get_resource_from_name(self.n_data_in).get_data()
        nk, nl, nt = data.shape
        data_np = data.cpu().numpy().astype(float)
        background = data_np.mean(axis=2) - self.planet_counts.mean(axis=2)
        variance = np.clip(background, 1e-3 * data_np.mean(axis=2), None)[:, :, None] + self.planet_counts
        weight_np = 1 / np.sqrt(variance)
        weight_torch = torch.as_tensor(weight_np, dtype=data.dtype, device=data.device)

        def transformation(x):
            if isinstance(x, np.ndarray):
                return (x - x.mean(axis=2, keepdims=True)) * weight_np
            return (x - x.mean(dim=2, keepdim=True)) * weight_torch

        r_data_out = DataResource(self.n_data_out)
        r_data_out.set_data(transformation(data))
        return r_data_out, TransformationResource(name=self.n_transformation_out, transformation=transformation,
                                                  dof_lost=nk * nl)


class FitStartModule(BaseModule):
    # Wraps a fitted SED and position as a planet resource, to start MLSEDEstimationModule (n_planets_in) there.
    # The module only reads planet.spectral_energy_distribution[:, 0, 0] and planet.sky_coordinates[i, 0, 0, 0, 0].

    def __init__(self, n_planets_out, sed, pos_x, pos_y):
        super().__init__()
        self.n_planets_out = n_planets_out
        self.start = SimpleNamespace(
            spectral_energy_distribution=torch.as_tensor(np.asarray(sed, dtype=float))[:, None, None],
            sky_coordinates=torch.tensor([float(pos_x), float(pos_y)]).reshape(2, 1, 1, 1, 1))

    def run(self, pipeline_resources):
        r_planets_out = ResourceCollection[PlanetResource](name=self.n_planets_out)
        r_planets_out.collection.append(PlanetResource(name='', planet=self.start))
        return r_planets_out,


def prepare_setup(target, p):
    pipeline = Pipeline(seed=None, gpu_index=GPU_INDEX, grid_size=GRID_SIZE)   

    # setup
    module = SetupModule(n_setup_out='setup', n_planets_out='params_init', instrument=make_instrument(p), 
                        observation=make_observation(p), scene=make_scene(target))
    pipeline.add_module(module)
    quiet(pipeline.run)
    phringe = pipeline.get_resource('setup').phringe

    # generate templates
    module = TemplateGenerationModule(n_setup_in='setup', n_template_out='temp', fov=float(phringe.get_field_of_view().max()) / 2)
    pipeline.add_module(module)
    quiet(pipeline.run)

    planet = pipeline.get_resource('params_init').collection[0].planet
    sed_true = planet.spectral_energy_distribution[:, 0, 0].cpu().numpy()
    x_true = float(planet.sky_coordinates[0, 0, 0, 0, 0])
    y_true = float(planet.sky_coordinates[1, 0, 0, 0, 0])

    def model(sed, x, y):
        return quiet(phringe.get_model_counts, spectral_energy_distribution=sed, x_position=x, y_position=y,
                     kernels=True)

    def model_centred(sed, x, y):
        counts = model(sed, x, y)
        return counts - counts.mean(axis=2, keepdims=True)

    # compute jacobian of the data vector
    n_wl = len(sed_true)
    step = 1e-4 * float(phringe.get_wavelength_bin_centers().min()) / p['B']    # ~1e-4 of a fringe [rad]
    columns = [model_centred(np.eye(n_wl)[j], x_true, y_true) for j in range(n_wl)]
    columns.append((model_centred(sed_true, x_true + step, y_true)
                    - model_centred(sed_true, x_true - step, y_true)) / (2 * step))
    columns.append((model_centred(sed_true, x_true, y_true + step)
                    - model_centred(sed_true, x_true, y_true - step)) / (2 * step))
    jacobian = np.stack(columns, axis=-1)                          # (n_kernels, n_wl, n_time, n_wl + 2)

    # Background level from a planet-free simulation (independent noise realisation)
    reference = PHRINGE(seed=None, gpu_index=GPU_INDEX, grid_size=GRID_SIZE)
    quiet(reference.set, make_instrument(p))
    quiet(reference.set, make_observation(p))
    quiet(reference.set, make_scene(target, with_planet=False))
    background_counts = quiet(reference.get_counts, kernels=True).cpu().numpy()

    return dict(pipeline=pipeline, model=model, sed_true=sed_true, x_true=x_true, y_true=y_true, jacobian=jacobian,
                sep_mas=float(np.hypot(x_true, y_true) * RAD_TO_MAS),
                # channels in which the planet is inside the fibre field of view (outside: no signal at all)
                in_fov=np.abs(jacobian[..., :n_wl]).sum(axis=(0, 1, 2)) > 0,
                planet_counts=model(sed_true, x_true, y_true),   # (n_kernels, n_wl, n_time)
                background=background_counts.mean(axis=2),       # (n_kernels, n_wl)
                var_over_mean=float(np.median(background_counts.var(axis=2, ddof=1)
                                              / background_counts.mean(axis=2))))


def fov_status(prep):
    n_in = int(prep['in_fov'].sum())
    n_wl = len(prep['in_fov'])
    return 'ok' if n_in == n_wl else 'outside FOV' if n_in == 0 else f'partly outside FOV ({n_in}/{n_wl} ch)'


def fisher_snr(prep):
    # Expected SNR (known, marginalised) at the true parameters, noise = background + planet photon noise.
    # Only the channels with the planet inside the field of view carry information; none -> SNR 0.
    in_fov = prep['in_fov']
    if not in_fov.any():
        return dict(known=0.0, marg=0.0)
    sed, n = prep['sed_true'][in_fov], int(in_fov.sum())
    sigma = np.sqrt(prep['background'][:, :, None] + prep['planet_counts'])
    j = (prep['jacobian'] / sigma[..., None])[..., np.r_[in_fov, True, True]].reshape(-1, n + 2)
    fisher = j.T @ j
    cov = np.linalg.inv(fisher)
    return dict(known=float(np.sqrt(sed @ fisher[:n, :n] @ sed)),
                marg=float(np.sqrt(sed @ np.linalg.solve(cov[:n, :n], sed))))


# Bulk SNR for the SED
def fit_record(fit, sed_true, n_wl, known):
    # s_hat C^-1 s_hat (and s C^-1 s with the true SED) of one fit; NaN if lmfit returned no covariance
    cov_full = fit.cov
    if cov_full is None or not np.all(np.isfinite(cov_full)):
        return dict(q=np.nan, q_true=np.nan)
    cov = cov_full if known else cov_full[:n_wl, :n_wl]
    sed = np.asarray(fit.sed)
    return dict(q=float(sed @ np.linalg.solve(cov, sed)), q_true=float(sed_true @ np.linalg.solve(cov, sed_true)))


# Accuracy of the position in units of it position covariance
def position_offset(fit, x_true, y_true):
    
    if fit.cov is None or not np.all(np.isfinite(fit.cov)):
        return np.nan
    cov_pos = fit.cov[-2:, -2:]

    fitted = np.array([fit.pos_x, fit.pos_y])
    true = np.array([x_true, y_true])

    return float(min(np.sqrt(d @ np.linalg.solve(cov_pos, d)) for d in (fitted - true, fitted + true)))


def run_realisation(prep, seed):
    pipeline, sed_true, model = prep['pipeline'], prep['sed_true'], prep['model']
    n_wl = len(sed_true)
    torch.manual_seed(seed)
    np.random.seed(seed)

    # --- Pass 1: the prototype ---
    pipeline.add_module(DataGenerationModule(n_setup_in='setup', n_data_out='data'))
    pipeline.add_module(ZCAWhiteningModule(n_setup_in='setup', n_data_in='data', n_template_in='temp',
                                           n_data_out='data_zca', n_template_out='temp_zca',
                                           n_transformation_out='zca', diagonal_only=True, subtract_mean=True))
    pipeline.add_module(MLSEDEstimationModule(n_setup_in='setup', n_data_in='data_zca', n_template_in='temp_zca',
                                              n_transformation_in='zca', n_planets_in='params_init',
                                              n_planets_out='p1_known', fix_position=True))
    pipeline.add_module(MLSEDEstimationModule(n_setup_in='setup', n_data_in='data_zca', n_template_in='temp_zca',
                                              n_transformation_in='zca', n_planets_out='p1_marg'))  # grid search
    quiet(pipeline.run)
    p1_known = pipeline.get_resource('p1_known').collection[0]
    p1_marg = pipeline.get_resource('p1_marg').collection[0]

    record = dict(seed=seed, n_wl=n_wl,
                  pos_offset_sigma=position_offset(p1_marg, prep['x_true'], prep['y_true']),
                  pos_x_mas=float(p1_marg.pos_x * RAD_TO_MAS), pos_y_mas=float(p1_marg.pos_y * RAD_TO_MAS))
    for case, fit in (('known', p1_known), ('marg', p1_marg)):
        r = fit_record(fit, sed_true, n_wl, case == 'known')
        record.update({f'proto_q_{case}': r['q'], f'proto_q_true_{case}': r['q_true']})

    # --- Pass 2: Poisson noise model including the planet, from the pass-1 fit of the same case ---
    # known: position fixed at the truth again (params_init); marginalised: started at the pass-1 best fit
    pipeline.add_module(FitStartModule(n_planets_out='p1_marg_start', sed=p1_marg.sed, pos_x=p1_marg.pos_x,
                                       pos_y=p1_marg.pos_y))
    passes = [('known', model(p1_known.sed, prep['x_true'], prep['y_true']), 'params_init'),
              ('marg', model(p1_marg.sed, p1_marg.pos_x, p1_marg.pos_y), 'p1_marg_start')]
    for case, planet_counts, start in passes:
        pipeline.add_module(PoissonWhiteningModule(n_data_in='data', n_data_out=f'data_poisson_{case}',
                                                   n_transformation_out=f'poisson_{case}',
                                                   planet_counts=planet_counts))
        # The template is only read for its shape (no grid search when a start planet is given)
        pipeline.add_module(MLSEDEstimationModule(n_setup_in='setup', n_data_in=f'data_poisson_{case}',
                                                  n_template_in='temp_zca', n_transformation_in=f'poisson_{case}',
                                                  n_planets_in=start, n_planets_out=f'p2_{case}',
                                                  fix_position=(case == 'known')))
    quiet(pipeline.run)
    for case in ('known', 'marg'):
        fit = pipeline.get_resource(f'p2_{case}').collection[0]
        r = fit_record(fit, sed_true, n_wl, case == 'known')
        record.update({f'q_{case}': r['q'], f'q_true_{case}': r['q_true']})
        if case == 'marg':
            record['pos_offset_sigma_p2'] = position_offset(fit, prep['x_true'], prep['y_true'])
    return record



def run_key(target, setup, seed):
    return json.dumps(dict(target=target, setup=SETUPS[setup], seed=seed, t_obs=T_OBS, t_mod=T_MOD, dit=DIT,
                           grid_size=GRID_SIZE, version=4), sort_keys=True)


done, PREP = {}, {}
for target_name, target in LOADED.items():
    if FORCE_RECOMPUTE:
        target['results_file'].unlink(missing_ok=True)   # start from an empty file
    elif target['results_file'].exists():
        for line in target['results_file'].read_text().splitlines():
            record = json.loads(line)
            done[record['key']] = record

    for setup_name, params in SETUPS.items():
        # The setup is prepared even if all realisations are on file: it gives the Fisher reference
        prep = quiet(prepare_setup, target, params)
        status = fov_status(prep)
        todo = [BASE_SEED + i for i in range(N_REALISATIONS)
                if run_key(target_name, setup_name, BASE_SEED + i) not in done] if status == 'ok' else []
        print(f"{target_name} / {setup_name}: separation {prep['sep_mas']:.0f} mas, {status}, "
              f"var/mean = {prep['var_over_mean']:.2f}, {len(todo)} realisations to run")
        start = time.time()
        for i, seed in enumerate(todo, 1):
            record = run_realisation(prep, seed)
            record.update(key=run_key(target_name, setup_name, seed), target=target_name, setup=setup_name)
            done[record['key']] = record
            with target['results_file'].open('a') as f:
                f.write(json.dumps(record) + '\n')
            elapsed = time.time() - start
            print(f'\r  [{i}/{len(todo)}] {elapsed / i:.0f} s per realisation, '
                  f'~{elapsed / i * (len(todo) - i) / 60:.0f} min left', end='')
        if todo:
            print()

        # Keep only what the summary needs; the pipeline (templates, PHRINGE caches, data) is several 100 MB
        PREP[target_name, setup_name] = dict(sed_true=prep['sed_true'], sep_mas=prep['sep_mas'], in_fov=prep['in_fov'],
                                             status=status, fisher=fisher_snr(prep))
        del prep
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()



def snr_from_q(q, n_wl):
    # bias-corrected SNR from the realisations of s_hat C^-1 s_hat, its standard error, and the per-realisation scatter
    q = np.asarray(q, dtype=float)
    q = q[np.isfinite(q)]
    if len(q) == 0:
        return np.nan, np.nan, np.nan
    snr = np.sqrt(max(q.mean() - n_wl, 0))
    err = q.std(ddof=1) / np.sqrt(len(q)) / (2 * snr) if snr > 0 and len(q) > 1 else np.nan
    return snr, err, np.sqrt(np.clip(q - n_wl, 0, None)).std(ddof=1) if len(q) > 1 else np.nan


rows, per_realisation, n_failed = [], [], {}
for (target_name, setup_name), prep in PREP.items():
    status = prep['status']
    fisher = prep['fisher']
    n_wl = len(prep['sed_true'])
    records = [done[run_key(target_name, setup_name, BASE_SEED + i)] for i in range(N_REALISATIONS)] \
        if status == 'ok' else []
    failed = [sum(not r[key] <= MAX_POS_OFFSET_SIGMA for r in records)
              for key in ('pos_offset_sigma', 'pos_offset_sigma_p2')]
    if any(failed):
        n_failed[target_name, setup_name] = failed
    for case in ('known', 'marg'):
        # A fit that did not locate the planet (marginalised: wrong peak, judged on the fit of each pass) or returned
        # no covariance is not a detection. Leaving it out would only keep the realisations where the noise happened
        # to look like the planet and bias the SNR up; instead it enters with q = n_wl, i.e. SNR 0 after bias correction.
        located = lambda r, key, q_key: (case == 'known' or r[key] <= MAX_POS_OFFSET_SIGMA) and np.isfinite(r[q_key])
        used = [r for r in records if located(r, 'pos_offset_sigma_p2', f'q_{case}')]
        q = [r[f'q_{case}'] if located(r, 'pos_offset_sigma_p2', f'q_{case}') else n_wl for r in records]
        q_proto = [r[f'proto_q_{case}'] if located(r, 'pos_offset_sigma', f'proto_q_{case}') else n_wl
                   for r in records]
        # the uncorrected diagnostics stay averages over the located fits only
        q_used = [r[f'q_{case}'] for r in used]
        q_proto_used = [r[f'proto_q_{case}'] for r in records if located(r, 'pos_offset_sigma', f'proto_q_{case}')]
        snr, err, scatter = snr_from_q(q, n_wl)
        rows.append(dict(
            target=target_name, setup=setup_name, case=case, status=status, n_used=len(used),
            SNR_fisher=fisher[case],
            proto_reported=np.mean(np.sqrt(q_proto_used)) if q_proto_used else np.nan,
            proto_corrected=snr_from_q(q_proto, n_wl)[0],
            SNR_uncorrected=np.mean(np.sqrt(q_used)) if q_used else np.nan,
            SNR=snr, SNR_err=err, SNR_scatter=scatter,
            SNR_true_sed=np.sqrt(np.nanmean([r[f'q_true_{case}'] for r in used])) if used else np.nan,
        ))
        per_realisation += [dict(target=target_name, setup=setup_name, case=case, snr=np.sqrt(max(v - n_wl, 0)))
                            for v in q if np.isfinite(v)]

summary = pd.DataFrame(rows)
for column in ('proto_reported', 'SNR'):
    summary[f'{column}/fisher'] = summary[column] / summary['SNR_fisher']
summary.to_csv(COMBINED_DIR / 'snr_summary.csv', index=False)
print('Free fits on a wrong peak (pass 1, pass 2), left out of the marginalised average of their pass:',
      n_failed or 'none')
summary.set_index(['target', 'setup', 'case']).round(3)




# --- What to show (see the table above) ---
CASE = 'marg'               # 'marg' (marginalised position) or 'known' (known position)
ESTIMATE = 'SNR'            # 'SNR', 'SNR_uncorrected', 'proto_corrected', 'proto_reported' or 'SNR_fisher'
SNR_DETECTION = 7           # a planet counts as detected if the shown SNR after T_OBS reaches this
SIGMA_SNR_DETECTION = 1     # how uncertain the SNR is allowed to be in SIGMA of SNR

CASE_LABELS = {'marg': 'marginalised position', 'known': 'known position'}
ESTIMATE_LABELS = {'SNR': 'pass 2, bias-corrected', 'SNR_uncorrected': 'pass 2, not bias-corrected',
                   'proto_corrected': 'pass 1 (prototype), bias-corrected',
                   'proto_reported': 'pass 1 (prototype) as reported', 'SNR_fisher': 'Fisher expectation'}

plt.rcParams.update({'savefig.dpi': 300})
names = list(SETUPS)
planets = sorted(TARGETS, key=lambda t: PREP[t, names[0]]['sep_mas'])    # smallest angular separation on the left
SETUP_LABELS = {n: f"{n}\nD = {p['D']:g} m, B = {p['B']:g} m\n"
                   f"{p['T_mirror']:g}/{p['T_instr']:g}/{p['T_det']:g} K\n"
                   f"null {p['null_depth']:.0e}, dark {p['dark_current']:g}"
                for n, p in SETUPS.items()}

rows = summary[summary.case == CASE]
snr_table = rows.pivot(index='setup', columns='target', values=ESTIMATE).loc[names, planets]
# per-observation scatter, only computed for the result (pass 2, bias-corrected)
err_table = rows.pivot(index='setup', columns='target', values='SNR_scatter').loc[names, planets] \
    if ESTIMATE == 'SNR' else snr_table * np.nan
status_table = rows.pivot(index='setup', columns='target', values='status').loc[names, planets]
sep_mas = {t: PREP[t, names[0]]['sep_mas'] for t in planets}
# realisations in which the single observation reaches the detection threshold (failed fits count as SNR 0)
realisations = pd.DataFrame(per_realisation)
realisations = realisations[realisations.case == CASE]
n_above = (realisations.assign(above=realisations.snr >= SNR_DETECTION)
           .pivot_table(index='setup', columns='target', values='above', aggfunc='sum')
           .reindex(index=names, columns=planets))

values = snr_table.where(status_table == 'ok').astype(float).values
sigma_snr = err_table.where(status_table == 'ok').astype(float).values


detected = np.nan_to_num(values - SIGMA_SNR_DETECTION * sigma_snr, nan=0) >= SNR_DETECTION      # outside the FOV counts as not detected
n_detected = detected.sum(axis=1)

fig = plt.figure(figsize=( len(planets) + 4, 0.85 * len(names) + 2.6), layout='constrained')
grid = fig.add_gridspec(2, 2, width_ratios=[0.5, len(planets)], height_ratios=[len(names), 0.5])
ax_count, ax, cax = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, 1])

# --- SNR per planet and setup ---
finite = values[np.isfinite(values) & (values > 0)]
# log scale from SNR 1; lower values (e.g. 0 after the noise-bias correction) get the lowest colour
#norm = LogNorm(vmin=1, vmax=max(finite.max(), 10) if finite.size else 10, clip=True)
norm = LogNorm(vmin=1, vmax=50, clip=True)
cmap = plt.get_cmap('Blues').copy()
cmap.set_bad('0.85')                                   # cells without a value: planet (partly) outside the FOV
image = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, norm=norm, aspect='auto')
for (row, col), value in np.ndenumerate(values):
    if np.isfinite(value):
        err = err_table.values[row, col]
        text = f'{value:.0f}' if value >= 10 else f'{value:.1f}'
        if np.isfinite(err):
            text += f'\n'  + rf'$\pm$ {err:.1f}'
        text += f'\n{n_above.values[row, col]:.0f}/{N_REALISATIONS} ≥ {SNR_DETECTION:g}'
        ax.text(col, row, text, ha='center', va='center', fontsize=9,
                color='white' if norm(value) > 0.6 else 'black')
    else:
        ax.text(col, row, status_table.values[row, col].replace(' (', '\n('), ha='center', va='center',
                fontsize=7.5, color='0.35')
ax.set_xticks(range(len(planets)), [f'{t}\n({sep_mas[t]:.0f} mas)' for t in planets], fontsize=9)
ax.set_yticks([])
ax.xaxis.tick_top()

# --- Detections per setup ---
ax_count.imshow(n_detected[:, None], vmin=-0.5, vmax = len(planets) + 0.5, cmap='magma',  aspect='auto') # vmin=-0.5, vmax=len(planets) + 0.5,
for row, n in enumerate(n_detected):
    ax_count.text(0, row, f'{n} ', ha='center', va='center', fontsize=10, fontweight='bold',
                  color='white' if n / len(planets) < 0.5 else 'black')
ax_count.set_xticks([0], [f'detected\nSNR ≥ {SNR_DETECTION:g}'], fontsize=9)
ax_count.xaxis.tick_top()
ax_count.set_yticks(range(len(names)), [SETUP_LABELS[n] for n in names], fontsize=7)


for a, n_cols in ((ax, len(planets)), (ax_count, 1)):
    a.tick_params(length=0)
    for spine in a.spines.values():
        spine.set_visible(False)
    a.set_xticks(np.arange(n_cols + 1) - 0.5, minor=True)
    a.set_yticks(np.arange(len(names) + 1) - 0.5, minor=True)
    a.grid(which='minor', color='white', lw=2)
    a.tick_params(which='minor', length=0)

# --- Colour bar below, with the detection threshold ---
colorbar = fig.colorbar(image, cax=cax, orientation='horizontal')
ticks = [t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000) if t <= norm.vmax]
colorbar.set_ticks(ticks, labels=[f'{t:g}' for t in ticks])
colorbar.set_label(f'bulk SNR after {T_OBS / 3600:g} h ({CASE_LABELS[CASE]}); red line: detection threshold',
                   fontsize=9)
cax.axvline(SNR_DETECTION, color='red', lw=4)

fig.suptitle(f'SNR, {CASE_LABELS[CASE]}, {ESTIMATE_LABELS[ESTIMATE]} ({N_REALISATIONS} realisations each; '
             f'fits that miss the planet count as SNR 0)\n'
             rf'cell: SNR $\pm$ scatter of a single observation, number of observations with SNR ≥ {SNR_DETECTION:g}',
             fontsize=10)
fig.savefig(COMBINED_DIR / f'snr_heatmap_{CASE}_{ESTIMATE}.png')#, transparent=True)
plt.show()