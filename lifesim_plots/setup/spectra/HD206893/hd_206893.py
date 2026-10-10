import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 155782207), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='HD 206893',
    distance='40.7583 pc',
    mass='1.3600000143051147 Msun',
    radius='1.311 Rsun',
    temperature='6554.3 K',
    right_ascension='21.756085 hourangle',
    declination='-12.783350 deg',
)


exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# TODO: point path_to_file at this planet's actual spectrum -- no catalogue provides one
sed_loader_hd_206893_c = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'PLACEHOLDER.txt'),
    sed_units='W/sr/m2/um',
    wavelength_units='um',
)
hd_206893_c = Planet(
    sed_loader=sed_loader_hd_206893_c,
    name='HD 206893 c',
    propagate_orbit=False,
    mass='3527.8953177 Mearth',
    radius='16.36511182 Rearth',
    temperature='1181.9 K',
    semi_major_axis='3.74 au',
    eccentricity=0.283,
    inclination='142.5 deg',
    argument_of_periapsis='27.0 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# TODO: point path_to_file at this planet's actual spectrum -- no catalogue provides one
sed_loader_hd_206893_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'PLACEHOLDER.txt'),
    sed_units='W/sr/m2/um',
    wavelength_units='um',
)
hd_206893_b = Planet(
    sed_loader=sed_loader_hd_206893_b,
    name='HD 206893 b',
    propagate_orbit=False,
    mass='6197.6539365 Mearth',
    radius='14.01122587 Rearth',
    temperature='1429.2 K',
    semi_major_axis='10.75 au',
    eccentricity=0.069,
    inclination='142.0 deg',
    argument_of_periapsis='102.0 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)
