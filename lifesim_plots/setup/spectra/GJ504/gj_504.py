import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 397587084), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='GJ 504',
    distance='17.5299 pc',
    mass='1.1200000047683716 Msun',
    radius='1.352 Rsun',
    temperature='6029.2 K',
    right_ascension='13.279587 hourangle',
    declination='9.424156 deg',
)


exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# TODO: point path_to_file at this planet's actual spectrum -- no catalogue provides one
sed_loader_gj_504_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'PLACEHOLDER.txt'),
    sed_units='W/sr/m2/um',
    wavelength_units='um',
)
gj_504_b = Planet(
    sed_loader=sed_loader_gj_504_b,
    name='GJ 504 b',
    propagate_orbit=False,
    mass='1271.3 Mearth',
    radius='12.10569916 Rearth',
    temperature='510.0 K',
    semi_major_axis='43.5 au',
    eccentricity=0.0,  # not available from the catalogue -- review this
    inclination='0 deg',  # not available from the catalogue -- review this
    raan='0 deg',  # not available from the catalogue -- review this
    argument_of_periapsis='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)
