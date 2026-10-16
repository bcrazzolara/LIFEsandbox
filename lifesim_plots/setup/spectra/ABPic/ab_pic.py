import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi


# Star parameters: TESS Input Catalog (TIC 260351540), queried via config/catalogue.py on 2026-10-16.
star = Star(
    name='AB Pic',
    distance='50.0475 pc',
    mass='0.8399999737739563 Msun',
    radius='1.019 Rsun',
    temperature='5027.0 K',
    right_ascension='6.320254 hourangle',
    declination='-58.054313 deg',
)


exozodi = Exozodi(
    level=3.0,  # not target-specific -- HOSTS survey median, review this
    inclination='0 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)


# Planet spectrum: no downloadable spectrum found (Bonnefoy et al. 2010 SINFONI, Palma-Bifani et al. 2023 are not public as
# files), so ATMO 2020 (Phillips et al. 2020) model, Teff=1700 K, log(g)=5.0, absolutely scaled (factor 7.574e-21) to the
# NACO Ks magnitude 14.14 (Chauvin et al. 2005). Teff matches the Exo-REM fit of Palma-Bifani et al. 2023
# (https://arxiv.org/abs/2211.01474; Bonnefoy et al. 2010 find 2000 K). Their log g is 4.5 +- 0.3, the local ATMO grid only has
# log g = 5.0 -- see atmo_spectra.py. The scale factor implies R ~ 1.9 Rjup at 50.05 pc, consistent with the radius below.
sed_loader_ab_pic_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'ab_pic_b.txt'),
    sed_units='W/m2/um',
    wavelength_units='um',
    observed_planet_radius='19.06 Rearth',
    observed_host_star_distance='50.0475 pc',
)

# Planet parameters: mass 13.5 Mjup (4290.5 Mearth, NASA Exoplanet Archive pscomppars / Chauvin et al. 2005 evolutionary mass).
# The archive's radius (12.3 Rearth) and temperature (15.8 K, an insolation Teq) are not the self-luminous values and were
# replaced: radius 1.7 Rjup (Exo-REM K-band fit, Palma-Bifani et al. 2023; evolutionary models give 1.53 Rjup), Teff 1700 K.
# Semi-major axis 190 au, seen edge-on, from the orbit fit of Palma-Bifani et al. 2023 (archive: 260 au projected separation).
# Eccentricity, raan, argument of periapsis and true anomaly are not constrained -- review this.
ab_pic_b = Planet(
    sed_loader=sed_loader_ab_pic_b,
    name='AB Pic b',
    propagate_orbit=False,
    mass='4290.5 Mearth',
    radius='19.06 Rearth',  # 1.7 Rjup
    temperature='1700.0 K',
    semi_major_axis='190 au',
    eccentricity=0.0,  # not available -- review this
    inclination='90 deg',  # edge-on, Palma-Bifani et al. 2023
    raan='0 deg',  # not available -- review this
    argument_of_periapsis='0 deg',  # not available -- review this
    true_anomaly='0 deg',  # not available -- review this
)
