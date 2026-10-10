import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 179413040), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='PDS 70',
    distance='113.064 pc',
    mass='0.6399999856948853 Msun',
    radius='1.186 Rsun',
    temperature='4117.0 K',
    right_ascension='14.136154 hourangle',
    declination='-41.397937 deg',
)


# WARNING: this system has a known/resolved circumstellar disk -- the default exozodi level below
# is almost certainly wrong for it. Look up a real disk brightness estimate before using this.
exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: ATMO 2020 (Phillips et al. 2020) model, log(g)=5.0, absolutely scaled (factor 3.419e-21) to a real measured magnitude -- see atmo_spectra.py for the log(g)/temperature-range caveats.
sed_loader_pds_70_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'pds_70_b.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='21.96960217 Rearth',
    observed_host_star_distance='113.064 pc',
)
pds_70_b = Planet(
    sed_loader=sed_loader_pds_70_b,
    name='PDS 70 b',
    propagate_orbit=False,
    mass='1017.0509024 Mearth',
    radius='21.96960217 Rearth',
    temperature='1392.0 K',
    semi_major_axis='20.8 au',
    eccentricity=0.17,
    inclination='131.0 deg',
    argument_of_periapsis='161.0 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: ATMO 2020 (Phillips et al. 2020) model, log(g)=5.0, absolutely scaled (factor 1.166e-20) to a real measured magnitude -- see atmo_spectra.py for the log(g)/temperature-range caveats.
sed_loader_pds_70_c = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'pds_70_c.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='22.19378179 Rearth',
    observed_host_star_distance='113.064 pc',
)
pds_70_c = Planet(
    sed_loader=sed_loader_pds_70_c,
    name='PDS 70 c',
    propagate_orbit=False,
    mass='2383.7130525 Mearth',
    radius='22.19378179 Rearth',
    temperature='1054.0 K',
    semi_major_axis='34.3 au',
    eccentricity=0.037,
    inclination='130.5 deg',
    argument_of_periapsis='53.0 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)
