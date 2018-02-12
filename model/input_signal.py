#!/usr/bin/python
# Author: cbecker@g.hmc.edu

import os
import sys
import numpy
import math
factor = 1.0e12
c = 3e8 # speed of light

def compute_signal(component_params, space, time, ns):
    z0 = range(0, int(space * factor), int(round(space * factor / ns)))
    z = map((lambda x: x/factor), z0)
    t0 = range(0, int(time * factor), int(round(time * factor / ns)))
    t = map((lambda x: x/factor), t0)
    signal = numpy.ones(len(z))
    # Get time from space
    # the amount of time over which we see the signal propagate
    for i in range(len(t)):
        sig_i = 0
        for j in range(1, len(component_params), 5):
            A = component_params[j]
            w = component_params[j+1]
            d = component_params[j+2]
            B = component_params[j+3]
            sig_i += A * math.sin(w*t[i] - w/c * z[i] + d) + B
        signal[i] = sig_i
    return [z, t, signal]

        
#-----------------------------------------------------------------------#

if __name__ == '__main__':

    # Read in file arguments
    args = sys.argv[1:]
    space = float(args[0])*1000 # the distance over which we see the signal propagate, we convert to meters
    time = float(args[1]) # time over which we plot the signal
    ns = float(args[2]) # number of samples
    component_params_file = args[3]

    # Get input data from .csv
    component_params = numpy.genfromtxt(component_params_file, delimiter = ',')[1:]
    new_comp_params = []
    for params in component_params:
        for param in params:
            new_comp_params.append(param)
    component_params = new_comp_params

    # Compute signal from components
    [z, t, signal] = compute_signal(component_params,  space, time, ns)
    
    # Export output to .csv
    signal_file = numpy.savetxt('input_signal.csv', (z, t, signal), delimiter = ',', fmt='%1.3f')
