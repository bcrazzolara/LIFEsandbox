import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 298810162), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='51 Eri',
    distance='29.7575 pc',
    mass='1.649999976158142 Msun',
    radius='1.487 Rsun',
    temperature='7295.0 K',
    right_ascension='4.626704 hourangle',
    declination='-2.473549 deg',
)


exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# TODO: point path_to_file at this planet's actual spectrum -- no catalogue provides one
sed_loader__51_eri_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'PLACEHOLDER.txt'),
    sed_units='W/sr/m2/um',
    wavelength_units='um',
)
_51_eri_b = Planet(
    sed_loader=sed_loader__51_eri_b,
    name='51 Eri b',
    propagate_orbit=False,
    mass='3464.3296363 Mearth',
    temperature='807.0 K',
    semi_major_axis='10.4 au',
    eccentricity=0.57,
    inclination='141.0 deg',
    argument_of_periapsis='82.0 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
    radius=None,  # TODO: REQUIRED, not available from the catalogue -- you must supply this
)
