#!/usr/bin/python
# Author: cbecker@g.hmc.edu

import os
import sys
import numpy
factor = 1.0e4

def compute_signal(component_params, time, space, ns):
    t = 1/factor * range(0, int(time * factor), round(time * factor / ns))
    z = 1/factor * range(0, int(space * factor), round(space * factor / ns))
    signal = 0
    for i in range(0, len(component_params), 4):
        A = component_params[i]
        w = component_params[i+1]
        d = component_params[i+2]
        B = component_params[i+3]
        phi_i = A * math.sin(w * t[i] - z[i] + d) + B
        signal += phi_i
    return [z, t, signal]

        
#-----------------------------------------------------------------------#

if __name__ == '__main__':

    # Read in file arguments
    args = sys.argv[1:]
    time = float(args[0]) # the amount of time over which we see the signal propagate
    space = float(args[1]) # the distance over which we see the signal propagate
    ns = args[2] # number of samples
    component_params_file = args[3]

    # Get input data from .csv
    component_params = numpy.genfromtxt(component_params_file, delimiter = ',')
    numpy.flatten(component_params)

    # Compute signal from components
    [z, t, signal] = compute_signal(component_params, time, space, ns)

    # Export output to .csv
    signal_file = numpy.savetxt('input_signal.csv', [z, t, signal], delimiter = ',')
