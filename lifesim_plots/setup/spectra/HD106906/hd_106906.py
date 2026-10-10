import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 390472975), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='HD 106906',
    distance='103.029 pc',
    mass='1.3200000524520874 Msun',
    radius='2.124 Rsun',
    temperature='6473.8 K',
    right_ascension='12.298109 hourangle',
    declination='-55.975526 deg',
)


# WARNING: this system has a known/resolved circumstellar disk -- the default exozodi level below
# is almost certainly wrong for it. Look up a real disk brightness estimate before using this.
exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# TODO: point path_to_file at this planet's actual spectrum -- no catalogue provides one
sed_loader_hd_106906_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'PLACEHOLDER.txt'),
    sed_units='W/sr/m2/um',
    wavelength_units='um',
)
hd_106906_b = Planet(
    sed_loader=sed_loader_hd_106906_b,
    name='HD 106906 b',
    propagate_orbit=False,
    mass='3496.0 Mearth',
    radius='12.4 Rearth',
    temperature='1800.0 K',
    semi_major_axis='650.0 au',
    eccentricity=0.0,  # not available from the catalogue -- review this
    inclination='0 deg',  # not available from the catalogue -- review this
    raan='0 deg',  # not available from the catalogue -- review this
    argument_of_periapsis='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)
