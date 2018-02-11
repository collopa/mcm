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
def N_z(z, angle_inc, recomb_coeff):
    """
        Gives charge number of electrons per m^3
    """
    n0 = 1
    kappa = c.gas_const*c.T_tropo/(c.p0 * c.grav)
    return kappa
