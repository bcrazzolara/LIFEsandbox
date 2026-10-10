import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 102076870), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='2MASS J12073346-3932539',
    distance='64.308 pc',
    mass='0.18799999356269836 Msun',
    radius='0.219 Rsun',
    temperature='2825.0 K',
    right_ascension='12.125963 hourangle',
    declination='-39.548338 deg',
)


exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: ATMO 2020 (Phillips et al. 2020) model, log(g)=5.0, absolutely scaled (factor 4.233e-21) to a real measured magnitude -- see atmo_spectra.py for the log(g)/temperature-range caveats.
sed_loader__2mass_j12073346_3932539_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'twa27b_Wm2um.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='15.69257298 Rearth',
    observed_host_star_distance='64.308 pc',
)
_2mass_j12073346_3932539_b = Planet(
    sed_loader=sed_loader__2mass_j12073346_3932539_b,
    name='2MASS J12073346-3932539 b',
    propagate_orbit=False,
    mass='1652.7077164 Mearth',
    radius='15.69257298 Rearth',
    temperature='1163.0 K',
    semi_major_axis='46.0 au',
    eccentricity=0.0,  # not available from the catalogue -- review this
    inclination='0 deg',  # not available from the catalogue -- review this
    raan='0 deg',  # not available from the catalogue -- review this
    argument_of_periapsis='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)
