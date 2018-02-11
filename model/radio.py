#!/usr/bin/python
# Author: cbecker@g.hmc.edu

#-----------------------------------------------------------------------#
'''
TODO:
learn about fft and ifft
fresnel(ocean_index_of_refr, iono_index_of_refr, angle_of_incidence, f0_w)
dispersion(f0_w_post_Fresnel)
'''
#-----------------------------------------------------------------------#


'''
Inputs: discrete signal vector as a function of time
        indices of refraction of ionosphere
        Angle of incidience (+ or - from parallel to Earth's surface)
        Distance to receiver
        Severity of Weather on a scale of 0-1, with 1 being the most intense
Outputs: discrete signal vector as function of time
         SNR at each N value         

Where f(t) = \sum_{n=0}^{N} f_{n}(t), this file:
-- decomposes the input signal into its Fourier components, f_{0_{n}}(t) --> f_{0_{n}}(w)
-- for each component, f0n(w), it computes fn(w) according to Fresnel's equations and dispersive equations
-- uses the inverse Fourier transform to recompose the output signal, f_{n}(w) --> f_{n}(t)

Notes:
-- Fresnel's equations require angle of incidence and indices of refraction of the two media.
-- The dispersive equations require the conductivity of air.
-- If angle of incidence is less than zero Earth's surface), odd N reflections are off of the ocean/land
-- If angle of incidence is greater than zero, even N reflections are off of the ionosphere
-- This model assumes Gaussian noise which accumulates over time, 
   so a random, normal noise vector is added to the signal vector every so often. 
-- The frequency of white noise addition is dependent on weather conditions. (Weather conditions on a 5 point scale)
'''

import os
import numpy as numpy
import math as math
import constants as constants
import ionosphere as iono_dict 

def get_num_refl(angle_of_incidence, distance):
    d = distance
    z = TxRx_height
    h = iono_height
    a = angle_of_incidence
    N = d * math.tan(a) / (h - z)
    return N

def calc_noise(weather_severity, N, array_length):
    base_stdev = array_length/100.0
    wc = 1/weather_severity # weather contribution
    Nc = math.sqrt(N) # Reflection number contribution
    stdev = base_stdev * wc * Nc
    noise = numpy.random.normal(0, stdev, array_length)
    return noise

def fresnel(ocean_index_of_refr, angle_of_incidence, f0_w):
    a = angle_of_incidence
    print constants.electron_charge


    
    
#-----------------------------------------------------------------------#

if __name__ == '__main__':

    print constants.electron_charge

    # Read in file arguments
    args = sys.argv[1:]
    ocean_index_of_refr = args[1]
    angle_of_incidence = args[2]
    distance =  args[3] # distance on Earth's surface to be travelled
    weather_severity = args[4]
    input_file = args[5] # input signal

    # Get input data from .csv
    # Decompose the input signal, f0_t, to it's Fourier components f0_w
    [time, f0_t] = numpy.genfromtxt(input_file, delimiter = ',')
    f0_w = fft(f0_t)

    # Number of reflections
    N = get_num_refl(angle_of_incidence, distance)

    # Calculate noise vector to be added to output signal
    noise = calc_noise(weather_severity, N, array_length)
    
    # For each f0_w, compute the output after Fresnel and dispersion
    f0_w_post_Fresnel = fresnel(ocean_index_of_refr, iono_index_of_refr, angle_of_incidence, f0_w, constants)
    f0_w_post_dispersion = dispersion(f0_w_post_Fresnel, constants)
    f_w = f0_w_post_dispersion

    # Compose output signal and add in accumulated noise
    # Export output data to .csv
    f_t_clean = ifft(f_w)
    f_t = f_t_clean + noise
    output_file = numpy.savetxt('radio_output.csv', [time, f_t], delimiter = ',')

    return output_file

    

    
