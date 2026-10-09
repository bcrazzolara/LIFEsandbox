import os as _os

from phringe.core.sources.planet import Planet
from phringe.io.sed_loader import SEDLoader
from phringe.core.sources.star import Star
from phringe.core.sources.exozodi import Exozodi
from phringe.core.sources.debris_disk import DebrisDisk


# Star parameters: TESS Input Catalog (TIC 270577175), queried via config/catalogue.py on 2026-09-14.
star = Star(
    name='bet Pic',
    distance='19.44 pc',
    mass='1.95 Msun',  
    radius='1.544 Rsun',
    temperature='8038.7 K',
    right_ascension='5.788080 hourangle',
    declination='-51.066511 deg',
)


# WARNING: this system has a known/resolved circumstellar disk -- the default exozodi level below
# is almost certainly wrong for it. Look up a real disk brightness estimate before using this.
exozodi = Exozodi(
    level=0.0,  # not target-specific -- HOSTS survey median, review this
    inclination='88.993 deg',  # unconstrained -- review this
    raan='0 deg',  # unconstrained -- review this
)

# Inner warm dust belts of the beta Pic debris disk. Belt positions 6.4, 16 and 30 au from mid-IR (10 um)
# spectroscopic imaging, Okamoto et al. 2004 (Nature 431, 660). Radial widths are assumed (not constrained there).
# Vertical aspect ratio 0.05 ~ thin dynamical population of beta Pic, 0.046 in the ARKS double-Gaussian fit
# (https://www.aanda.org/articles/aa/full_html/2026/07/aa59709-26/aa59709-26.html); the disk is ~1.5x thicker in
# the mid-IR than at mm. Inclination/raan match bet Pic b (coplanar; the inner disk is warped by a few deg towards
# the planet's orbit). Optical depths are PLACEHOLDERS -- calibrate against the observed IR excess before use.
# At 4-5 um only the 6.4 au belt contributes noticeably (~190 K blackbody grains; tau=1e-4 gives ~1e-4 of the
# stellar flux, ~10% of bet Pic b); the 16/30 au belts are >1e4x fainter there. temperature_factor=1 assumes
# blackbody grains; the sub-micron silicates seen by Okamoto et al. are hotter (superheated), so this is a lower limit.
disks = [
    DebrisDisk(
        name=f'bet Pic belt {r0} au',
        radial_profile='gaussian',
        reference_radius=f'{r0} au',
        radial_width=f'{width} au',            # assumed
        inner_radius=f'{r0 - 2 * width} au',   # truncate at +/-2 sigma to avoid an unphysical hot inner tail
        outer_radius=f'{r0 + 2 * width} au',
        optical_depth=1e-4,                     # placeholder -- review this
        inclination='0.993 deg',               # = bet Pic b
        raan='0 deg',                           # = bet Pic b (placeholder there too)
        scale_height=0.05,                      # edge-on: must integrate along the line of sight
        temperature_factor=1.0,
    )
    for r0, width in [(6.4, 2.0), (16, 3.0), (30, 5.0)]
]

# Planet parameters: Teff, log g, [M/H], C/O and radius from Houlle et al. 2025 (https://arxiv.org/pdf/2508.18366),
# Table 3, GRAVITY+MATISSE (fixed) model. Semi-major axis/eccentricity/inclination/argument of periapsis from
# NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# Planet spectrum: Exo-REM (cloud_R500 grid), quadrilinear-interpolated to Teff=1529 K, log g=3.83, [M/H]=0.30,
# C/O=0.539 -- native per-steradian intensity, no observed_planet_radius/observed_host_star_distance needed.
sed_loader_bet_pic_b = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)),  'betapicb_exorem_interpolated_spectrum.dat'),
    sed_units='W/m2/um/sr',
    wavelength_units='um',
)


bet_pic_b = Planet(
    sed_loader=sed_loader_bet_pic_b,
    name='bet Pic b',
    propagate_orbit=False,
    mass='2765.11 Mearth',
    radius='23.37 Rearth',  # 2.085 Rjup, Houlle et al. 2025, consistent with the Exo-REM SED above
    temperature='1529.0 K',             # https://arxiv.org/pdf/2508.18366
    semi_major_axis='10.07 au',
    eccentricity=0.105,
    inclination='88.993 deg',
    argument_of_periapsis='23.7 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)



'''
# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Radial Velocity, queried via config/catalogue.py on 2026-09-14.
# TODO: point path_to_file at this planet's actual spectrum -- no catalogue provides one
sed_loader_bet_pic_c = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'PLACEHOLDER.txt'),
    sed_units='W/sr/m2/um',
    wavelength_units='um',
)
bet_pic_c = Planet(
    sed_loader=sed_loader_bet_pic_c,
    name='bet Pic c',
    propagate_orbit=False,
    mass='2765.1071409 Mearth',
    radius='12.5 Rearth',
    temperature='1250.0 K',
    semi_major_axis='2.7 au',
    eccentricity=0.307,
    inclination='88.93 deg',
    argument_of_periapsis='-114.9 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)


# Planet parameters: NASA Exoplanet Archive (pscomppars), discovered via Imaging, queried via config/catalogue.py on 2026-09-14.
# TODO: point path_to_file at this planet's actual spectrum -- no catalogue provides one
sed_loader_bet_pic_d = SEDLoader(
    path_to_file=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'spectra', 'PLACEHOLDER.txt'),
    sed_units='W/sr/m2/um',
    wavelength_units='um',
)
bet_pic_d = Planet(
    sed_loader=sed_loader_bet_pic_d,
    name='bet Pic d',
    propagate_orbit=False,
    mass='762.7881768 Mearth',
    radius='14.12331568 Rearth',
    temperature='600.0 K',
    semi_major_axis='26.0 au',
    eccentricity=0.44,
    inclination='89.0 deg',
    argument_of_periapsis='191.0 deg',
    raan='0 deg',  # not available from the catalogue -- review this
    true_anomaly='0 deg',  # not available from the catalogue -- review this
)
'''
