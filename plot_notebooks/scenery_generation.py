import numpy as np

def generate_scenery(sep_planet_as, phase_planet, flux_planet, inc_exozodi, phi_exozodi, r0_exozodi, sig_exozodi, flux_exozodi, ALPHA, BETA):

    IM = np.zeros_like(ALPHA)
    Npix = ALPHA.shape[0]
    fov_as = np.max(ALPHA)


    # apply rotation
    ALPHA_rot = ALPHA * np.cos(np.deg2rad(phi_exozodi)) + BETA * np.sin(np.deg2rad(phi_exozodi))
    BETA_rot = 1/np.cos(np.deg2rad(inc_exozodi)) * (-ALPHA * np.sin(np.deg2rad(phi_exozodi)) + BETA * np.cos(np.deg2rad(phi_exozodi)))

    R_rot2 = ALPHA_rot**2 + BETA_rot**2
    IM = flux_exozodi * np.exp(-0.5 * ((np.sqrt(R_rot2) - r0_exozodi) / sig_exozodi)**2) / np.cos(np.deg2rad(inc_exozodi))


    # add star 
    IM[int(1/2*Npix), int(1/2*Npix)] += 0

    # add planet
    idx_planet = int(1/2*Npix + sep_planet_as/fov_as * 1/2*Npix * np.cos(np.deg2rad(phase_planet)))
    idy_planet = int(1/2*Npix + sep_planet_as/fov_as * 1/2*Npix * np.sin(np.deg2rad(phase_planet)))
    IM[idx_planet, idy_planet] += flux_planet

    


    return IM

