#!/usr/bin/python
# Authors: cadams@g.hmc.edu, rbarclay@g.hmc.edu, cbecker@g.hmc.edu

##### Default Parameters #####
T_tropo =  288.  # [Kelvin]
T_stato = 270.   # [Kelvin]
T_meso = 183.    # [Kelvin]
T_thermo = 1773. # [Kelvin]
T_iono = 350.    # [Kelvin]

gas_const = 8.3144598    # J/mol*Kelvin
grav = 9.80665      # m/s^2
p0 = 1.225       # kg/m^3
air_dielec_const = 1.000 # [Coulomb^2 Second^2 Meter^-3 Kilogram^-1]
MW_air = 0.02895           # [Kilorams Mole^-1]

#---------------------------------------------------------------#

def get_default():
    return {
            'T_therm': T_therm,
            'T_meso': T_meso,
            'T_tropo': T_tropo,
            'T_strato': T_strato ,
            'T_iono': T_iono, 
            'air_dielec_const': air_dielec_const,
            'MW_air': MW_air,
            'gas_const': gas_const,
            'R': gas_const,
            'grav': grav,
            'p0': p0 
            }
}

#---------------------------------------------------------------#
if __name__ == '__main__':

        constants = get_default()
