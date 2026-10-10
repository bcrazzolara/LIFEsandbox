import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 176945734), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='kap And',
    distance='50.0177 pc',
    mass='2.819999933242798 Msun',
    radius='2.312 Rsun',
    temperature='10839.0 K',
    right_ascension='23.673474 hourangle',
    declination='44.333933 deg',
)


exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: ATMO 2020 (Phillips et al. 2020) model, log(g)=5.0, absolutely scaled (factor 5.129e-21) to a real measured magnitude -- see atmo_spectra.py for the log(g)/temperature-range caveats.
sed_loader_kap_and_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'kap_and_b.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='15.91675259 Rearth',
    observed_host_star_distance='50.0177 pc',
)
kap_and_b = Planet(
    sed_loader=sed_loader_kap_and_b,
    name='kap And b',
    propagate_orbit=False,
    mass='5498.4314411 Mearth',
    radius='15.91675259 Rearth',
    temperature='1791.0 K',
    semi_major_axis='55.0 au',
    eccentricity=0.0,  # not available from the catalogue -- review this
    inclination='0 deg',  # not available from the catalogue -- review this
    raan='0 deg',  # not available from the catalogue -- review this
    argument_of_periapsis='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)
