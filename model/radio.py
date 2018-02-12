#!/usr/bin/python
# Author: cbecker@g.hmc.edu

#-----------------------------------------------------------------------#
'''
TODO:
learn about fft and ifft
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
-- The frequency of white noise addition is dependent on weather conditions. (Weather conditions on a 0-1 scale)
'''

import os
import numpy as numpy
import math as math
import constants as constants
import ionosphere_output as iono_dict 

def decompose(time, f0_t):
    max_time = max(time)
    ns = len(time) # number of samples
    L = ns + 1
    Fs = max_time/ns # sampling frequency
    min_freq = -ns/2
    max_freq = ns/2 - ns/L
    step_size = ns/L
    k = min_freq:step_size:max_freq
    f0_k = numpy.fft(f0_t)
    return [k, f0_k]  

def get_num_refl(angle_of_incidence, distance):
    d = distance
    z = constants.TxRx_height
    h = constants.iono_height
    a = angle_of_incidence
    N = math.floor(d * math.tan(a) / (h - z))
    return N

def calc_noise(weather_severity, N, array_length):
    base_stdev = array_length/100.0
    wc = 1/weather_severity # weather contribution
    Nc = math.sqrt(N) # Reflection number contribution
    stdev = base_stdev * wc * Nc
    noise = numpy.random.normal(0, stdev, array_length)
    return noise

def fresnel(ocean_index_of_refr, angle_of_incidence, [k, f0_k], N):    

    # Constants
    theta_I = angle_of_incidence
    theta_R = theta_I # By the law of reflection
    n1 = (constants.air_dielec_const)**2
    n2_ocean = ocean_index_of_refr

    # Determine multipliers, which determine the amplification of Fresnel effects
    if N % 2 == 0:
        iono_multiplier = N/2
        ocean_multiplier = N/2
    else:
        if theta_I > 0:
            iono_multiplier = (N - 1) / 2 + 1
            ocean_multiplier = (N - 1) / 2
        else:
            iono_multiplier = (N - 1) / 2 
            ocean_multiplier = (N - 1) / 2 + 1
    refl_tup = ((n2_iono, iono_multiplier), (n2_ocean, ocean_multiplier))
        
    # Compute Fresnel-modified signal
    f0_k_post_Fresnel = numpy.empty())
    for f0n_k in f0_k:
        mag_sig = abs(f0n_k)
        exp_sig = f0n_k/mag_sig
        fn_k = 1
        for tup in refl_tup:
            n2 = tup[0]
            m = tup[1]
            theta_T = numpy.arcsin(n1 * numpy.sin(theta_I) / n2) # By Snell's law
            beta = n2/n1
            alpha = numpy.cos(theta_T) / numpy.cos(theta_I)
            fn_k = m * (alpha - beta) / (alpha + beta) * f0n_k
        numpy.append(f0_k_post_Fresnel, fn_k)

    # Return the Fresnel-modified signal
    return f0_k_post_Fresnel

def dispersion([k, f0_k_post_Fresnel]):
    sigma = constants.air_elec_cond
    f0_k_post_dispersion = numpy.empty(f0_k_post_Fresnel)
    omega = 2 * math.pi * k # convert to angular frequency
    for n in range(len(k)):
        wn = w[n]
        f0n_k = f0_k_post_Fresnel[n]
        kappa_n = wn * math.sqrt(constants.epsilon_perm * constants.mu_perm / 2) * math.sqrt(1 + math.sqrt(1 + (constants.air_elec_cond/ comstants.epsilon_perm / wn)**2)))
        new_f0_k = f0n_k * math.exp(-kappa)
    return f0_k_post_dispersion

    
    
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
    # Decompose the input signal, f0_t, to it's Fourier components f0_k
    [time, f0_t] = numpy.genfromtxt(input_file, delimiter = ',')
    [k, f0_k] = decompose(time, f0_t)
    numpy.savetxt('FT_input_signal.csv', [k, f0_k], delimiter = ',')

    # Number of reflections
    N = get_num_refl(angle_of_incidence, distance)

    # Calculate noise vector to be added to output signal
    noise = calc_noise(weather_severity, N, array_length)
    
    # For each f0_k, compute the output after Fresnel and dispersion
    f0_k_post_Fresnel = fresnel(ocean_index_of_refr, iono_index_of_refr, angle_of_incidence, [k, f0_k] , N)
    f0_k_post_dispersion = dispersion([k, f0_k_post_Fresnel])
    f_k = f0_k_post_dispersion

    # Compose output signal and add in accumulated noise
    # Export output data to .csv
    f_t_clean = ifft([time, f_k])
    f_t = f_t_clean + noise
    output_file = numpy.savetxt('radio_output.csv', [time, f_t], delimiter = ',')
    
    

    
