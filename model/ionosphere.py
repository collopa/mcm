#########################
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
def chargeDenHeight(z, angle_inc, recomb_coeff):
    """
        Gives charge number of electrons per m^3
    """
    n0 = 1
    kappa = c.R*c.T_tropo/(c.p0 * c.g)
    return kappa