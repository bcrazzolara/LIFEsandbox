# LIFEsandbox

Research sandbox for the **LIFE** (Large Interferometer For Exoplanets) mission, using a **single Bracewell (SBW)** nulling interferometer. Owner: bcrazzolara (ETH). GitHub: `bcrazzolara/LIFEsandbox`, branch `master`. Not a package: scripts, notebooks and data folders, no tests, no build, no requirements file.

There are two unrelated parts. **`lifesim_plots/` is the active work.** The root-level SBW code is a legacy side project.

## Environment

- Conda env **`life`**: `/scratch/user/bcrazzolara/miniconda3/envs/life` (Python 3.12, torch 2.14, numpy 1.26, astropy, lmfit, pandas, scipy, matplotlib). The default `python` outside the env does **not** have `phringe`/`lifesimmc`; use `/scratch/user/bcrazzolara/miniconda3/envs/life/bin/python`.
- Machine `rainbow`: 7x RTX 2080 Ti (GPU 0-6), shared. GPU index is **hardcoded** in code (`GPU_INDEX = 3` in `instrument_tradeoffs.py`, `gpu_index=2` in notebooks / `prototype_pipeline.py`); `None` = CPU. Check `nvidia-smi` for a free GPU before choosing one.
- Two sibling repos (your forks, installed **editable** into `life`) are the real dependencies. They live outside this repo; **read them, do not modify them unless asked**:
  - `/scratch/user/bcrazzolara/PHRINGE` (`phringe` 1.8.2): the forward simulator. `phringe.core.{instrument,observation,scene}`, `phringe.core.sources.{star,planet,exozodi,local_zodi,debris_disk}`, `phringe.lib.{array_configuration,beam_combiner}` (`SBWConfiguration.acm`, `SingleBracewell.catm/.kernels`), `phringe.io.sed_loader.SEDLoader`, `phringe.main.PHRINGE` (`.set(...)`, `.get_counts(kernels=True)`, `.get_model_counts(...)`, `.get_wavelength_bin_centers()`, `.get_field_of_view()`). Recent work there: static null depth, debris disk, thermal/dark-current noise fluxes.
  - `/scratch/user/bcrazzolara/LIFEsimMC` (`lifesimmc` 1.3.5): data-analysis pipeline on top of PHRINGE. `Pipeline(seed, gpu_index, grid_size)` runs named **modules** that read/write named **resources** in a dict (`pipeline.get_resource(name)`). Modules in `lifesimmc/core/modules/{loading,generating,processing}`.
- Heavy runs: **ask the user before running simulations or anything on a GPU.** Light checks (syntax, imports, `python -c`, small demos) are fine. Do not launch the full trade-off sweep on your own.
- Notebooks are normally run by the user in VS Code with the `life` kernel.

## Layout

```
LIFEsandbox/
├── lifesim_plots/                 # ACTIVE: PHRINGE + LIFEsimMC studies
│   ├── instrument_tradeoffs.py    # MAIN study (.ipynb is the same code as notebook)
│   ├── prototype_pipeline.py/.ipynb   # single-target prototype run + plots
│   ├── detector_budget.ipynb      # single-target (bet Pic b) detector/noise budget
│   ├── target_sky_map.ipynb       # Hammer sky map of the target stars
│   ├── full_sbw_pipeline.ipynb, SBW_lifesimmc_example_advanced.ipynb  # older "Sun Twin" demos
│   ├── colormap.py                # custom blue/red colormaps (pos_map, neg_map, full_map)
│   ├── setup/
│   │   ├── prototype_setup.py     # `inst`, `obs`: the baseline prototype Instrument/Observation
│   │   ├── config/dh_tau.py       # stray old-style config (spectra/ is the real location)
│   │   └── spectra/<Target>/      # one folder per host star, see "Targets"
│   └── plots/                     # outputs: setup_snr/{snr_summary.csv, heatmaps, detections_grid}, target_sky_map.png, ...
├── sbw.py, gen_data.py, train.py, invert.py, app.py, transmission_maps.py, SBW.ipynb   # LEGACY toy NN project
└── plot_notebooks/                # LEGACY analytic SBW plots (NCC, transmission maps, signal GIFs)
```

`lifesim_plots/` scripts use paths relative to `lifesim_plots/` (`from setup.prototype_setup import inst, obs`, `setup/spectra/...`, `plots/setup_snr`). **Run them with `lifesim_plots/` as the working directory**, otherwise imports and output paths break. Notebooks do the same (kernel cwd = `lifesim_plots/`).

## The physics / measurement in one paragraph

A single Bracewell interferometer (2 apertures, baseline B, null-depth `null_depth`) rotates for `T_OBS` (5 h, half of `T_MOD` = 10 h rotation period), observing 4-5 um at R=100 (about 22 wavelength channels) in `DIT` = 120 s steps. A planet off-axis produces a modulated signal in (wavelength x time) per kernel. The question studied: **for which directly imaged planets does a given instrument (aperture D, nulling baseline B, thermal state) reach bulk SNR >= 7 on the planet spectrum?** Noise = photon noise of star leakage, local zodi, exozodi, debris disks, mirror/instrument thermal emission, detector dark current.

## lifesim_plots: how the analysis works

All of this is in `instrument_tradeoffs.py` (the most important file, about 700 lines, top-to-bottom script / notebook):

1. **Config constants** (top): `T_OBS, T_MOD, DIT, N_REALISATIONS, BASE_SEED, FORCE_RECOMPUTE, MAX_POS_OFFSET_SIGMA, GRID_SIZE, GPU_INDEX`.
2. **`TARGETS`** dict: name -> `dict(config=<spectra/.../x.py>, planet=<variable name in that file>)`; optional `sed_file/sed_units` overrides a missing spectrum. `load_target()` execs the config file via importlib (with a patched `SEDLoader` so a missing `PLACEHOLDER.txt` does not crash) and returns star, exozodi, disks, planet and the results/plot paths.
3. **Setups**: `GRID_D x GRID_B x GRID_THERMAL` ('cold'/'warm': mirror/instrument/detector temperature, dark current 0.035 / 0.35 e-/s/px, `NULL_DEPTH=1e-4`) -> `GRID_SETUPS`. `SETUPS` is the subset (D=0.2 m, all B) used for the heatmap. `ALL_SETUPS = {**SETUPS, **GRID_SETUPS}`. The old named setups (Large/Medium/Prototype/Small/Short cold/warm) sit in a dead triple-quoted string. `make_instrument(p)` / `make_observation(p)` / `make_scene(target)` build PHRINGE objects (instrument fixed: 4-5 um, R=100, throughput 0.12, QE 0.7, NIRSpec-like detector).
4. **`prepare_setup`** (once per setup x target): runs `SetupModule` (+ `TemplateGenerationModule` once per setup, reused across planets), builds the true planet SED/position, a finite-difference **Jacobian** of the noise-free model w.r.t. (SED per channel, x, y), and a planet-free background simulation. `in_fov` flags channels where the planet is inside the field of view; `fov_status` -> 'ok' / 'partly outside FOV' / 'outside FOV'. `fisher_snr` gives the Fisher-matrix SNR reference (known / marginalised position).
5. **`run_realisation(prep, seed)`**: one noise realisation, two passes.
   - Pass 1 ("prototype"): `DataGenerationModule` -> `ZCAWhiteningModule(diagonal_only, subtract_mean)` -> `MLSEDEstimationModule` twice (known position `fix_position=True`; marginalised = grid-search position).
   - Pass 2: custom `PoissonWhiteningModule` (noise sigma^2 = background + planet counts from the pass-1 fit) and `FitStartModule` (injects a fit as start planet) -> `MLSEDEstimationModule` again.
   - Result record: `q = s_hat^T C^-1 s_hat` (and with true SED), position offset in sigma, for known/marg and pass 1/2.
6. **Caching**: results are appended as JSON lines to `setup/spectra/<Target>/setup_snr/realisations_<slug>.jsonl`. The record key is `json.dumps(dict(target, setup params, seed, t_obs, t_mod, dit, grid_size, version=4))`. Finished (target, setup, seed) combos are reused, and identical setup parameters share results across names. **If you change anything that alters the physics/analysis, bump `version` in `run_key`** (or set `FORCE_RECOMPUTE=True`, which deletes the jsonl files, so ask first). Changing instrument constants inside `make_instrument` (not in the key) silently reuses stale results.
7. **Summary / plots**: `snr_from_q` = bias-corrected SNR `sqrt(max(mean(q) - n_wl, 0))`. Failed fits (position further than `MAX_POS_OFFSET_SIGMA` from +/- truth, or no covariance) count as `q = n_wl` (SNR 0), not dropped. A planet counts as detected if `SNR - SIGMA_SNR_DETECTION*scatter >= SNR_DETECTION (7)`; outside FOV = not detected. Writes `plots/setup_snr/snr_summary.csv`, `snr_heatmap_<CASE>_<ESTIMATE>.png`, `detections_grid_*.png`. `CASE` is 'marg' or 'known'; `ESTIMATE` is 'SNR', 'SNR_uncorrected', 'proto_corrected', 'proto_reported' or 'SNR_fisher'. Note: `vmax = len(planets) - 5` in the heatmap colour scales is a manual tweak, not a bug.

Cost: about 45 s per realisation on CPU (much less on GPU). The full grid is 50 setups x 17 targets x `N_REALISATIONS`, hence the caching.

### `prototype_pipeline.py` / `.ipynb`
Runs the single baseline `inst, obs` from `setup/prototype_setup.py` (D=0.25 m, B=0.75 m, 120/80 K, 5 h) on one target (currently bet Pic b, switch by editing the import/`planet = ...` lines near the top). Steps: scene (star + LocalZodi + exozodi + disks + planet) -> SetupModule -> DataGenerationModule + TemplateGenerationModule -> ZCA whitening -> `MLSEDEstimationModule` (free position and `fix_position=True`) -> metrics (r/sigma_r, 2pi/sigma_phi, position accuracy, flux SNR known/marg, sqrt(chi2/n) goodness of fit) -> plots saved **into the target's spectrum folder** (`data.png`, `cost_function.png`, `flux_estimate.png`, `position_estimate.png`, `snr_summary.txt`). Note the +/-(x, y) position ambiguity of a single Bracewell: the plots pick whichever mirror image of the truth is closer. `detector_budget.ipynb` is the same flow for bet Pic b with detector-noise analysis.

## Targets (`lifesim_plots/setup/spectra/<Target>/`)

Each folder holds a config `.py` plus spectrum `.txt` files and a `setup_snr/` results folder. Config file convention (copy an existing one, e.g. `GJ504/gj_504.py` or `BetaPic/bet_pic.py`):

- `star = Star(name, distance, mass, radius, temperature, right_ascension, declination)` (TIC catalogue values), `exozodi = Exozodi(level, inclination, raan)`, optional `disks = [DebrisDisk(...)]`, one `Planet(...)` per planet with `SEDLoader(path_to_file, sed_units, wavelength_units, [observed_planet_radius, observed_host_star_distance])` and orbital elements (NASA Exoplanet Archive + a chosen `true_anomaly` for a 2030-07 epoch, comments record the separation in mas and the whereistheplanet check).
- Variable names matter: `TARGETS` refers to the planet by Python variable name (e.g. `_51_eri_b`, `_2mass_j12073346_3932539_b`: leading underscore because names cannot start with a digit).
- SED units: `W/m2/um/sr` (surface brightness) needs no radius/distance; `W/m2/um` (flux at the observer) **requires** `observed_planet_radius` and `observed_host_star_distance`.
- Many configs were auto-generated by a catalogue script (`config/catalogue.py`, not in this repo) with `# review this` comments: placeholders for exozodi level (3.0 = HOSTS median), inclination, raan, argument of periapsis. Several still point to `spectra/PLACEHOLDER.txt` (no spectrum), such as DH Tau, which then fails in `load_target` unless `sed_file` is given.
- Active targets in `TARGETS` (17): AB Pic b, PDS70 b/c, HR 8799 b/c/d/e, 51 Eri b, bet Pic b, Kap And b, DH Tau b, HD 106906 b, HD 4113 C, HD 95086 b, GJ 504 b, HIP 65426 b, 2MASS J1207 b. `VHS1256/` has only a spectrum and is not yet wired in. `HD206893/` has a config but is not in `TARGETS`.
- bet Pic has a resolved edge-on disk: three `DebrisDisk` belts with placeholder optical depths (marked in comments).

## Legacy: root SBW toy project (do not touch unless asked)

Self-contained NumPy/PyTorch experiment, no PHRINGE. Idealised SBW with an 80x80 px scene (exozodi ring + one planet pixel) projected through a precomputed transmission "scenery" (`transmission_maps.point_signal_SBW`) to a (31 wavelengths x 108 time steps) signal.
- `sbw.py`: constants, `sample_params`, `build_scene`, `transmission_scenery`, `ViT` (patch 4, 2 input channels), `demo()` self-check.
- `gen_data.py` -> `data/{train,test,stats}.npz` (10k/500 samples); `train.py` -> `model.pt`; `app.py` Flask debug UI (its title says "ConvNext" but the model is a ViT); `invert.py` closed-form least-squares inversion (recovers 6 params exactly on noiseless data; `demo()` and `noise_test()`).
- The 6 parameters are idx1, idx2, planet amplitude, exozodi amplitude, inc, phi. idx2 is restricted to a half-plane because SBW cannot distinguish a point from its 180-degree reflection. This is the same +/-(x, y) ambiguity seen in the PHRINGE analysis.
- `plot_notebooks/` has analytic derivations (NCC curvatures, Bessel-function position uncertainty), animated GIFs and publication-style plots. `plot_notebooks/transmission_maps.py` is identical to the root one; the two `colormap.py` copies differ only in slice indices and which one is pos/neg.
- `data/` and `model.pt` are not in the repo; regenerate with `gen_data.py` / `train.py`.

## Working conventions

- **Notebooks are large (0.2-8 MB) with embedded outputs.** Prefer editing the `.py` equivalent where one exists (`instrument_tradeoffs.py`, `prototype_pipeline.py`) and say if the notebook needs the same change; edit `.ipynb` only with NotebookEdit/JSON-aware tools, never by text replace, and never dump whole notebooks into context (extract cell sources with a small script).
- The `.py` and `.ipynb` versions of the same analysis drift apart (the notebooks were last run at different times). Check both before assuming they match.
- Code style in the repo: dense, heavily commented scripts, explanatory comments in plain English that say *why*; matplotlib figures saved with `transparent=True`, dpi 300. Follow that style; keep PHRINGE/LIFEsimMC calls wrapped in `quiet(...)` where the existing code does (the modules print a lot).
- Units: PHRINGE uses SI internally (radians for sky positions, `RAD_TO_MAS`); pass astropy quantities or unit strings (`'19.44 pc'`) to Star/Planet.
- Tracked in git: the result `.jsonl` files and `__pycache__` are committed (noisy diffs). The working tree usually has many uncommitted modified spectra/results files; **never `git add -A`, never commit or push unless asked**, and do not revert files you did not change.
- Pointers for context outside this repo: PHRINGE and LIFEsimMC sources (read-only) answer almost every "what does this module do / what are the parameters" question faster than guessing.
