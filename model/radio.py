#!/usr/bin/python
# Author: cbecker@g.hmc.edu

#-----------------------------------------------------------------------#
'''
TODO:
import iono and ocean indices correctly, change MAIN!
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
import sys
import numpy as numpy
import math as math
import constants as constants
from operator import add
factor = 1.0e4
c = 3e8 # speed of light

def get_z(t, f0_t):
    z = numpy.empty(len(t))
    for i in range(len(t)):
        zi = t[i] - math.arcsin(f0_t[i])
        numpy.append(z, zi)
    return z

def decompose(t, f0_t):
    max_time = max(t)
    ns = float(len(t)) # number of samples
    Fs = max_time/ns # sampling frequency
    min_freq = -ns/2
    max_freq = ns/2 
    step_size = 1
    k0 = range(int(min_freq*factor), int(max_freq*factor), int(step_size*factor))
    k = map((lambda x: x/factor), k0)
    f0_k = numpy.fft.fft(f0_t)
    return [k, f0_k]  

def get_num_refl(angle_of_incidence, distance):
    d = distance
    z = constants.TxRx_height
    h = constants.iono_height
    a = angle_of_incidence
    N = math.floor(d * math.tan(a) / (h - z))
    return N

def calc_noise(weather_severity, N, array_length, f_t_real):
    sum_f_t = 0.0
    for num in f_t_real: sum_f_t += num
    avg = sum_f_t/array_length
    base_stdev = array_length/100.0
    wc = weather_severity + 0.1 # weather contribution
    Nc = math.sqrt(N) # Reflection number contribution
    stdev = base_stdev * wc * Nc
    noise = numpy.random.normal(avg, stdev, array_length)
    return noise

def add_fresnel(iono_index_of_refr, ocean_index_of_refr, angle_of_incidence, k, f0_k, N):    

    # Constants
    theta_I = angle_of_incidence
    theta_R = theta_I # By the law of reflection
    n1 = (constants.air_dielec_const)**2
    n2_iono = iono_index_of_refr
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
    f0_k_post_Fresnel = numpy.ones(len(k))
    for n in range(len(k)):
        fn_k = f0_k[n]
        for tup in refl_tup:
            n2 = tup[0]
            m = tup[1]
            theta_T = numpy.arcsin(n1 * numpy.sin(theta_I) / n2) # By Snell's law
            # By Fresnel's equations:
            beta = n2/n1
            alpha = numpy.cos(theta_T) / numpy.cos(theta_I)
            fn_k = (alpha - beta) / (alpha + beta) * fn_k # *m
            f0_k_post_Fresnel[n] = fn_k

    # Return the Fresnel-modified signal
    return f0_k_post_Fresnel

def add_dispersion(z, t, k, f0_t, f0_k_post_Fresnel):
    sigma = constants.air_elec_cond
    f0_k_post_dispersion = numpy.ones(len(z))
    omega = numpy.ones(len(k))
    kappa = numpy.ones(len(k))
    for i in range(len(k)):
        omega[i] = 2 * math.pi * k[i] # convert to angular frequency
    for n in range(len(k)):
        zn = z[n]
        wn = omega[n]
        if wn != 0:
            kappa_n = wn * math.sqrt(constants.epsilon_perm * constants.mu_perm / 2) \
                         * math.sqrt(1 \
                         + math.sqrt(1 \
                         + (constants.air_elec_cond/ constants.epsilon_perm / wn)**2))
            fn_k = math.exp(-kappa_n * zn) * f0_k_post_Fresnel[n]
        else:
            kappa_n = 0
            fn_k = f0_k_post_Fresnel[n]
        f0_k_post_dispersion[n] = fn_k
        kappa[n] = kappa_n
    return [kappa, f0_k_post_dispersion]
    
#-----------------------------------------------------------------------#

if __name__ == '__main__':

    # Read in file arguments
    args = sys.argv[1:]
    iono_index_of_refr = float(args[0])
    ocean_index_of_refr = float(args[1])
    angle_of_incidence = math.radians(float(args[2])) # in degrees, we convert to radians
    distance = float(args[3])*1000 # distance on Earth's surface to be traveled in kilometers, we convert to meters
    weather_severity = float(args[4])
    input_file = args[5] # input signal

    # Get input data from .csv
    # Decompose the input signal, f0_t, to it's Fourier components f0_k
    [z, t, f0_t] = numpy.genfromtxt(input_file, delimiter = ',')
    [k, f0_k] = decompose(t, f0_t)
    k_real = map((lambda x: numpy.real(x)), k)
    f0_k_real = map((lambda x: numpy.real(x)), f0_k)
    numpy.savetxt('FT_input_signal.csv', [k_real, f0_k_real], delimiter = ',', fmt='%1.3f')

    # Number of reflections
    N = get_num_refl(angle_of_incidence, distance)
    
    # For each f0_k, compute the output after Fresnel and dispersion
    f0_k_post_Fresnel = add_fresnel(iono_index_of_refr, ocean_index_of_refr, angle_of_incidence, k, f0_k, N)
    [kappa, f0_k_post_dispersion] = add_dispersion(z, t, k, f0_t, f0_k_post_Fresnel)
    f_k = f0_k_post_dispersion
    numpy.savetxt('FT_output_signal.csv', [k_real, f_k, kappa], delimiter = ',', fmt='%1.4e')
    
    # Compose output signal and add in accumulated noise
    # Export output data to .csv
    f_t_real = numpy.real(numpy.fft.ifft(f_k))
   
    # Calculate noise vector to be added to output signal
    noise = calc_noise(weather_severity, N, len(t), f_t_real)
    f_t = map(add, f_t_real, noise)

    numpy.savetxt('radio_output.csv', (z, t, f_t), delimiter = ',', fmt='%1.3f')
    
    

    
