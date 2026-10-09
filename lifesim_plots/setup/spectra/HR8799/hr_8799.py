import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 245368902), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='HR 8799',
    distance='41.2441 pc',
    mass='1.6100000143051147 Msun',
    radius='1.493 Rsun',
    temperature='7204.6 K',
    right_ascension='23.124643 hourangle',
    declination='21.134253 deg',
)


# WARNING: this system has a known/resolved circumstellar disk -- the default exozodi level below
# is almost certainly wrong for it. Look up a real disk brightness estimate before using this.
exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)



# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: ATMO 2020 (Phillips et al. 2020) model, log(g)=5.0, absolutely scaled (factor 2.757e-21) to a real measured magnitude -- see atmo_spectra.py for the log(g)/temperature-range caveats.
sed_loader_hr_8799_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'hr_8799_b.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='13.0 Rearth',
    observed_host_star_distance='41.2441 pc',
)
hr_8799_b = Planet(
    sed_loader=sed_loader_hr_8799_b,
    name='HR 8799 b',
    propagate_orbit=False,
    mass='2000.0 Mearth',
    radius='13.0 Rearth',
    temperature='1200.0 K',
    semi_major_axis='68.0 au',
    inclination='30.27 deg',
    eccentricity=0.0,  # not available from the catalogue -- review this
    raan='0 deg',  # not available from the catalogue -- review this
    argument_of_periapsis='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: ATMO 2020 (Phillips et al. 2020) model, log(g)=5.0, absolutely scaled (factor 6.434e-21) to a real measured magnitude -- see atmo_spectra.py for the log(g)/temperature-range caveats.
sed_loader_hr_8799_c = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'hr_8799_c.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='13.0 Rearth',
    observed_host_star_distance='41.2441 pc',
)
hr_8799_c = Planet(
    sed_loader=sed_loader_hr_8799_c,
    name='HR 8799 c',
    propagate_orbit=False,
    mass='3000.0 Mearth',
    radius='13.0 Rearth',
    temperature='1200.0 K',
    semi_major_axis='38.0 au',
    eccentricity=0.5,
    inclination='29.43 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    argument_of_periapsis='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)



# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: ATMO 2020 (Phillips et al. 2020) model, log(g)=5.0, absolutely scaled (factor 5.535e-21) to a real measured magnitude -- see atmo_spectra.py for the log(g)/temperature-range caveats.
sed_loader_hr_8799_d = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'hr_8799_d.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='13.0 Rearth',
    observed_host_star_distance='41.2441 pc',
)
hr_8799_d = Planet(
    sed_loader=sed_loader_hr_8799_d,
    name='HR 8799 d',
    propagate_orbit=False,
    mass='3000.0 Mearth',
    radius='13.0 Rearth',
    temperature='1300.0 K',
    semi_major_axis='24.0 au',
    eccentricity=0.6,
    inclination='38.63 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    argument_of_periapsis='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: ATMO 2020 (Phillips et al. 2020) model, log(g)=5.0, absolutely scaled (factor 8.234e-21) to a real measured magnitude -- see atmo_spectra.py for the log(g)/temperature-range caveats.
sed_loader_hr_8799_e = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'hr_8799_e.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='13.11453 Rearth',
    observed_host_star_distance='41.2441 pc',
)
hr_8799_e = Planet(
    sed_loader=sed_loader_hr_8799_e,
    name='HR 8799 e',
    propagate_orbit=False,
    mass='3178.3 Mearth',
    radius='13.11453 Rearth',
    temperature='1150.0 K',
    semi_major_axis='16.4 au',
    eccentricity=0.15,
    inclination='25.0 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    argument_of_periapsis='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)
