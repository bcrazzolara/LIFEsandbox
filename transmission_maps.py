import numpy as np


# Vk is the complex field amplitude for a telescope sitting at the position pos, in the focal plane of 
# the interferometer. The focal plane is perpendicular to the distance vector of the telescope to the source.
def Vk(pos = (0, 0), lambda_ref = 10e-6, alpha = np.linspace(-0.1, 0.1, 1000), beta = np.linspace(-0.1, 0.1, 1000)):
    xk, yk = pos
    ALPHA, BETA = np.meshgrid(alpha * (np.pi/180) / 3600, beta * (np.pi/180) / 3600, indexing='ij') # [rad]
    fringe_pattern = np.exp(2*np.pi*1j*(xk*ALPHA + yk*BETA)/lambda_ref)
    return fringe_pattern 

# =================== T_(lambda, phi) (alpha, beta) =====================================
def transmission_map_SBW(baseline, rot_angle, lambda_ref = 10e-6, alpha = np.linspace(-0.1, 0.1, 1000), beta = np.linspace(-0.1, 0.1, 1000)):
    '''
    baseline: distance between the two telescopes
    rot_angle: rotation angle of the interferometer in radians
    lambda_ref: reference wavelength in meters
    alpha: array of alpha angles [as]
    beta: array of beta angles [as]
    
    Returns the transmission map of the single bracewell interferometer
    '''
    xk = baseline/2 * np.cos(rot_angle) 
    yk = baseline/2 * np.sin(rot_angle) 
    V1 = Vk(pos=(xk, yk), lambda_ref=lambda_ref, alpha=alpha, beta=beta)
    V2 = Vk(pos=(-xk, -yk), lambda_ref=lambda_ref, alpha=alpha, beta=beta)
    D1 = 1/2*(V1 + np.exp(1j*np.pi)*V2)

    TSBW = np.abs(D1)**2
    return TSBW


def transmission_map_DBW(baseline, rot_angle, q = 6, lambda_ref = 10e-6, alpha = np.linspace(-0.1, 0.1, 1000), beta = np.linspace(-0.1, 0.1, 1000)):
    '''
    baseline: distance between the two telescopes
    rot_angle: rotation angle of the interferometer in radians
    q: ratio of the imaging baseline to the nulling baseline
    lambda_ref: reference wavelength in meters
    alpha: array of alpha angles [as]
    beta: array of beta angles [as]

    Returns the transmission map of the double bracewell interferometer
    '''
    xk = baseline/2
    yk = q*baseline/2
    V1 = Vk(pos=(xk*np.cos(rot_angle) - yk*np.sin(rot_angle), yk*np.cos(rot_angle) + xk*np.sin(rot_angle)), lambda_ref=lambda_ref, alpha=alpha, beta=beta)
    V2 = Vk(pos=(-xk*np.cos(rot_angle) - yk*np.sin(rot_angle), yk*np.cos(rot_angle) - xk*np.sin(rot_angle)), lambda_ref=lambda_ref, alpha=alpha, beta=beta)
    V3 = Vk(pos=(xk*np.cos(rot_angle) + yk*np.sin(rot_angle), -yk*np.cos(rot_angle) + xk*np.sin(rot_angle)), lambda_ref=lambda_ref, alpha=alpha, beta=beta)
    V4 = Vk(pos=(-xk*np.cos(rot_angle) + yk*np.sin(rot_angle), -yk*np.cos(rot_angle) - xk*np.sin(rot_angle)), lambda_ref=lambda_ref, alpha=alpha, beta=beta)


    D1 = 1/np.sqrt(4)*(V1 + np.exp(1j*np.pi)*V2)
    D2 = 1/np.sqrt(4)*(V3 + np.exp(1j*np.pi)*V4)

    O1 = 1/np.sqrt(4)*(D1 + np.exp(1j*np.pi/2)*D2)  
    O2 = 1/np.sqrt(4)*(D1 - np.exp(1j*np.pi/2)*D2)

    TDBW = np.abs(O1)**2 - np.abs(O2)**2
    return TDBW


# =================== S_(alpha, beta) (lambda, t) =====================================
def point_signal_SBW(t, alpha, beta, baseline, lamda_ref = 10e-6, Omega = 1e-4):
    '''
    t: time in seconds
    alpha: array of alpha angles [as]
    beta: array of beta angles [as]
    baseline: distance between the two telescopes
    lamda_ref: reference wavelength in meters
    
    Returns the time signal of a point source at (alpha, beta) as a function of time and wavelength for a single bracewell
    '''
    alpha_rad = alpha * (np.pi/180) / 3600
    beta_rad = beta * (np.pi/180) / 3600
    alpha_t_rad = np.multiply.outer(np.cos(Omega * t) * baseline, alpha_rad) - np.multiply.outer(np.sin(Omega * t) * baseline, beta_rad)
    return np.sin(np.pi  * np.multiply.outer(1/lamda_ref, alpha_t_rad)) ** 2


def point_signal_DBW(t, alpha, beta, baseline, lamda_ref = 10e-6, Omega = 1e-4, q = 6):
    '''
    t: time in seconds
    alpha: array of alpha angles [as]
    beta: array of beta angles [as]
    baseline: distance between the two telescopes
    lamda_ref: reference wavelength in meters
    q: ratio of the imaging baseline to the nulling baseline
    
    Returns the time signal of a point source at (alpha, beta) as a function of time and wavelength for a double bracewell
    '''
    alpha_rad = alpha * (np.pi/180) / 3600
    beta_rad = beta * (np.pi/180) / 3600
    alpha_t_rad = np.multiply.outer(np.cos(Omega * t) * baseline, alpha_rad) - np.multiply.outer(np.sin(Omega * t) * baseline, beta_rad)
    beta_t_rad = np.multiply.outer(np.sin(Omega * t) * baseline, alpha_rad) + np.multiply.outer(np.cos(Omega * t) * baseline, beta_rad)
    return np.sin(np.pi * np.multiply.outer(1/lamda_ref, alpha_t_rad))**2 * np.sin(2 * np.pi * q  * np.multiply.outer(1/lamda_ref, beta_t_rad))
