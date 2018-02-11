#!/usr/bin/python
# Authors: cadams@g.hmc.edu, rbarclay@g.hmc.edu, cbecker@g.hmc.edu

##### Default Parameters #####
##### Temperatures #####
T_tropo =  288.  # [Kelvin]
T_stato = 270.   # [Kelvin]
T_meso = 183.    # [Kelvin]
T_thermo = 1773. # [Kelvin]
T_iono = 350.    # [Kelvin]

##### Ionosphere #####
gas_const = 8.3144598 # [Joule Mole^-1 Kelvin^-1]
grav = 9.80665 # [Kilogram Meter^-3]
p0 = 1.225 # [Kilogram Meter^-3]
air_dielec_const = 1.000 # [Coulomb^2 Second^2 Meter^-3 Kilogram^-1]
MW_air = 0.02895 # [Kiloram Mole^-1]
electron_charge = 1.6e-19 #[Coulomb]
mag_field = 50e-6 #[Tesla]

##### Radio #####
iono_height = 300e3 # [Meter]
TxRx_height = 100. # [Meter]

#---------------------------------------------------------------#

def get_default():
    return {
        'T_thermo': T_thermo,
        'T_meso': T_meso,
        'T_tropo': T_tropo,
        'T_strato': T_strato,
        'T_iono': T_iono, 
        'gas_const': gas_const,
        'grav': grav,
        'p0': p0
        'air_dielec_const': air_dielec_const,
        'MW_air': MW_air,
        'electron_charge': electron_charge,
        'mag_field': mag_field,
        'iono_height': iono_height,
        'TxRx_height': TxRx_height
    }
#---------------------------------------------------------------#
if __name__ == '__main__':

        constants = get_default()
