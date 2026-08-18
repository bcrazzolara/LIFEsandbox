# Image with planet - averaged over time 
import numpy as np

def cross_corr_coeff(signal, template, centered = True):
    S = signal.copy()
    T = template.copy() 
    if centered:
        S -= np.mean(signal, axis = 1)[:, None] # average along time axis
        T -= np.mean(template, axis = 1)[:,None,:,:] # average along time axis

    return np.sum(T * S[:, :, None, None], axis = 1) / np.linalg.norm(T, axis = 1) / np.linalg.norm(S, axis = 1)[:, None, None]

def cross_corr_coeff_merged(signal, template, centered = True):
    S = signal.copy()
    T = template.copy() 
    if centered:
        S   -=np.mean(signal, axis = 1)[:, None] # average along time axis
        T -= np.mean(template, axis = 1)[:,None,:,:] # average along time axis

    return np.sum(T * S[:, :, None, None], axis = (0,1)) / np.linalg.norm(T, axis = (0,1)) / np.linalg.norm(S)