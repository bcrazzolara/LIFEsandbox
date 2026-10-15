import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 118586897), queried via config/catalogue.py on 2026-10-15.
star = Star(
    name='HD 4113',
    distance='41.8727 pc',
    mass='1.0 Msun',
    radius='1.085 Rsun',
    temperature='5655.4 K',
    right_ascension='0.720165 hourangle',
    declination='-37.982633 deg',
)


exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet spectrum: no published spectrum available for download (Cheetham et al. 2018 SPHERE YJH spectrum, GRAVITY+ K band
# not public), so ATMO 2020 (Phillips et al. 2020) model, Teff=600 K, log(g)=5.0, absolutely scaled (factor 3.611e-21) to the
# SPHERE/IRDIS H2 magnitude 19.7 (Cheetham et al. 2018, Table 4). Teff 500-600 K, log g 4.5-5 from their atmospheric fits
# (COND photometric temperatures are 580-690 K); log g = 5.0 is the only gravity in the local ATMO grid -- see atmo_spectra.py.
sed_loader_hd_4113_c = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'hd_4113_c.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='15.7 Rearth',
    observed_host_star_distance='41.8727 pc',
)

# Planet parameters: Cheetham et al. 2018 (https://arxiv.org/abs/1712.05217): dynamical mass 66 Mjup, radius 1.4-1.5 Rjup
# (1.4 used), 535 mas projected separation = 22.4 au at 41.87 pc. Semi-major axis is set to the projected separation and the
# orbit to circular face-on as a simplification (the real orbit is eccentric, e~0.44 Cheetham 2018 / ~0.65 and a~50 au
# Rickman et al. 2022 as reported secondhand) -- review this. Not in the NASA Exoplanet Archive.
hd_4113_c = Planet(
    sed_loader=sed_loader_hd_4113_c,
    name='HD 4113 C',
    propagate_orbit=False,
    mass='20976 Mearth',  # 66 Mjup
    radius='15.7 Rearth',  # 1.4 Rjup
    temperature='600.0 K',
    semi_major_axis='22.4 au',  # projected separation -- review this
    eccentricity=0.0,  # simplification -- review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
    argument_of_periapsis='0 deg',  # unconstrained -- review this
    true_anomaly='0 deg',  # unconstrained -- review this
)
