#!/bin/bash                                                                                            
# Author: cbecker@g.hmc.edu

NS=100 # number of samples
DISTANCE=500 # distance from LA to Tokyo, 8820 km
TIME=30 # time over which we want to observe the signal
ANGLE=30 # angle of incidence (from parallel to Earth) in degrees
WEATHER=0 # weather severity on a scale of 0-1 where 1 is most severe

python input_signal.py $DISTANCE $TIME $NS 'component_params.csv'
python radio.py '1.5' '1.5' $ANGLE $DISTANCE $WEATHER 'input_signal.csv'
