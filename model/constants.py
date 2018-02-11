#!/usr/bin/python
# Authors: cadams@g.hmc.edu, rbarclay@g.hmc.edu, cbecker@g.hmc.edu

##### Default Parameters #####
T_tropo =  288 # [Kelvin]
T_stato = 270 # [Kelvin]
T_meso = 183 # [Kelvin]
T_thermo = 1773 # [Kelvin]
T_iono = 350 # [Kelvin]

air_dielec_const = 1.000 # [Coulomb^2 Second^2 Meter^-3 Kilogram^-1]
MW_air = 28.95 # [Grams Mole^-1]

#---------------------------------------------------------------#

def get_default():
    return {
        'temperature':{
            'T_therm': T_therm,
            'T_meso': T_meso,
            'T_tropo': T_tropo,
            'T_strato': T_strato ,
            'T_iono': T_iono
            },
        'atmosphere':{
            'air_dielec_const': air_dielec_const,
            'MW_air': MW_air
            }
}

#---------------------------------------------------------------#
if __name__ == '__main__':

        constants = get_default()
