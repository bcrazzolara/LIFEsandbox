

from phringe.core.instrument import Instrument
from phringe.core.observation import Observation

from phringe.lib.array_configuration import SBWConfiguration
from phringe.lib.beam_combiner import SingleBracewell
from phringe.core.perturbations.power_law_psd_perturbation import PowerLawPSDPerturbation

import astropy.units as u


# Use the predefined ideal LIFE baseline instrument, i.e. without any instrumental noise
inst = inst = Instrument(
    # Interferometer 
    array_configuration_matrix=SBWConfiguration.acm,
    complex_amplitude_transfer_matrix=SingleBracewell.catm,
    kernels=SingleBracewell.kernels,
    nulling_baseline_min=0. * u.m,
    nulling_baseline_max=5 * u.m,
    throughput=0.12,

    # Science Wavelengths
    wavelength_min=4 * u.um,
    wavelength_max=5 * u.um,
    wavelength_bands_boundaries=[],  # Alternatively e.g. ['9 um']
    spectral_resolving_power=100,       # simulation resolving power
    spectral_resolving_power_instrument=100,# "real" instrument spectral power, needed if simulation spectral_resolving_power is reduced to save comp. power

    # Instrumental Instability Noise
    amplitude_perturbation=None,
    phase_perturbation=None,
    polarization_perturbation=None,
    #amplitude_perturbation=PowerLawPSDPerturbation(coefficient=1, rms=0.000001 * u.percent),
    #phase_perturbation=PowerLawPSDPerturbation(coefficient=1, rms=0.000015 * u.nm, chromatic=True),
    #polarization_perturbation=PowerLawPSDPerturbation(coefficient=1, rms=0.000001 * u.rad),

    # Static on-axis null depth I_null/I_bright (reduced fringe visibility); None = perfect null
    null_depth=1e-3,

    # Mirror 
    aperture_diameter=0.25 * u.m,
    mirror_temperature=150 * u.K,  
    mirror_emissivity=0.06,
    mirror_throughput=0.94,

    # Pre-fiber/instrument environment: emissivity here is a genuine configurable parameter
    # (LIFEsim's reference default is 1.0, i.e. treated as a perfect emitter).
    instrument_temperature=120 * u.K, 
    instrument_emissivity=1.0,

    # Detector thermal background: pixel_size/detector wavelength range are LIFEsim's reference
    # Based on detector properties of JWST/NIRSpec, PRISM/CLEAR
    # https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-operations/nirspec-fs-operations/nirspec-fs-wavelength-ranges-and-gaps#gsc.tab=0
    # https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-detectors/nirspec-detector-performance#gsc.tab=0
    # https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-detectors#gsc.tab=0
    quantum_efficiency=0.7,
    detector_temperature=70 * u.K,          # NIRSpec SCAs at 42.8K
    detector_wavelength_min=0.6 * u.um,     # NIRSpec
    detector_wavelength_max = 5.3 * u.um,   # NIRSpec
    pixel_size=18e-6 * u.m,                 # NIRSpec
    pixels_per_wavelength_bin=2.2,          # assumption
    dark_current=1,                         # e-/px/s
)


# User the predefined observation for the LIFE baseline design
obs = Observation(
    solar_ecliptic_latitude=0 * u.deg,  # Alternatively: '0 deg' or 0
    total_integration_time=10*3600* u.s,  # Alternatively: '1 d' or 86400
    detector_integration_time='60 s',  # Alternatively: '600 s' or 600
    modulation_period=10*3600* u.s,  # Alternatively: '1 d' or 86400
    nulling_baseline=0.5
)
