import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 438702139), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='HIP 65426',
    distance='108.875 pc',
    mass='2.299999952316284 Msun',
    radius='1.735 Rsun',
    temperature='9097.0 K',
    right_ascension='13.410027 hourangle',
    declination='-51.504457 deg',
)


exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: ATMO 2020 (Phillips et al. 2020) model, log(g)=5.0, absolutely scaled (factor 1.516e-21) to a real measured magnitude -- see atmo_spectra.py for the log(g)/temperature-range caveats.
sed_loader_hip_65426_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'hip_65426_b.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='16.8135 Rearth',
    observed_host_star_distance='108.875 pc',
)
hip_65426_b = Planet(
    sed_loader=sed_loader_hip_65426_b,
    name='HIP 65426 b',
    propagate_orbit=False,
    mass='2860.47 Mearth',
    radius='16.8135 Rearth',
    temperature='1500.0 K',
    semi_major_axis='92.0 au',
    eccentricity=0.0,  # not available from the catalogue -- review this
    inclination='0 deg',  # not available from the catalogue -- review this
    raan='0 deg',  # not available from the catalogue -- review this
    argument_of_periapsis='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)
