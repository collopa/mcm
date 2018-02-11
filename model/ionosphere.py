###########################
# Packages
###########################
from scipy import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

##############################
# Dependencies
#################################
import constants as c


############################
# Functions
#############################
def N_z(z, angle_inc, sigma):
    """
        Gives charge number of electrons per m^3
    """
    n0 = 1
    # sigma = 4.032e18
    kappa = c.gas_const*c.T_tropo/(c.p0 * c.grav)
    z0 = kappa*np.log(sigma*kappa*c.p0) 
    inside = 0.5*(1 - (z - z0)/kappa - 1./np.cos(angle_inc) * np.exp(-(z - z0)/kappa))  
    return n0 * np.exp(inside)

def refIndex(freq, height, angle_inc, sigma):
    """
        Refractive index for an isotropic medium without electron damping.
        Based on a function of atmospheric height. Frequency of RF.
    """
    N = N_z(height, angle_inc, sigma)
    ang_freq = 2*np.pi*freq
    X = N * c.electron_charge**2/(c.free_perm * e_mass * ang_freq**2)
    n2 = 1 - X
    
    return np.sqrt(n2)

