import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 399637637), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='HD 95086',
    distance='86.2279 pc',
    mass='1.7699999809265137 Msun',
    radius='1.487 Rsun',
    temperature='7595.0 K',
    right_ascension='10.950839 hourangle',
    declination='-68.667347 deg',
)


# WARNING: this system has a known/resolved circumstellar disk -- the default exozodi level below
# is almost certainly wrong for it. Look up a real disk brightness estimate before using this.
exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: ATMO 2020 (Phillips et al. 2020) model, log(g)=5.0, absolutely scaled (factor 2.234e-21) to a real measured magnitude -- see atmo_spectra.py for the log(g)/temperature-range caveats.
sed_loader_hd_95086_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'hd_95086_b.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='12.9 Rearth',
    observed_host_star_distance='86.2279 pc',
)
hd_95086_b = Planet(
    sed_loader=sed_loader_hd_95086_b,
    name='HD 95086 b',
    propagate_orbit=False,
    mass='1589.0 Mearth',
    radius='12.9 Rearth',
    temperature='1000.0 K',
    semi_major_axis='55.7 au',
    eccentricity=0.0,  # not available from the catalogue -- review this
    inclination='0 deg',  # not available from the catalogue -- review this
    raan='0 deg',  # not available from the catalogue -- review this
    argument_of_periapsis='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)
