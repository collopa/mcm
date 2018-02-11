import os
import clawpack.geoclaw.topotools as topotools
import matplotlib.pyplot as plt
#from numpy import genfromtxt
import numpy
import pandas
import csv
import math
from functools import reduce
product = reduce((lambda x, y: x * y), [1, 2, 3, 4])
#from netCDF4 import Dataset

df = pandas.read_csv("OSMC_flattened_eff7_eca9_6785.csv", parse_dates=[0], header=0,delimiter=",")


def getColData(column):

	avgs4 = []
	nDate = df['time'][0][0:10]
	dayData = []
	numMeas = 0

for i in range(df['time'].count()):

	np1Date = df['time'][i][0:10]

	if nDate == np1Date:
		dayTempData.append(df['ztmp'][i])
		daySalData.append(df['zsal'][i])
		numMeas +=1 #increment number of measurements in each day
			
	else: #new day- divide into 4 chunks and average each
		secondStop = math.floor(numMeas/2)
		firstStop = math.floor(secondStop/2)
		thirdStop = math.floor(numMeas/2) + secondStop

		l1 = list(filter(lambda x: not math.isnan(x), daySalData[0:firstStop]))
		l2 = list(filter(lambda x: not math.isnan(x), daySalData[0:firstStop]))
		l3 = list(filter(lambda x: not math.isnan(x), daySalData[0:firstStop]))
		l4 = list(filter(lambda x: not math.isnan(x), daySalData[0:firstStop]))

		if len(l1) != 0:
			avg1Sal = sum(l1)/len(l1)
		else:
			avg1Sal = 'nan'

		if len(l1) != 0:
			avg2Sal = sum(l2)/len(l2)
		else:
			avg2Sal = 'nan'

		if len(l1) != 0:
			avg3Sal = sum(l3)/len(l3)
		else:
			avg3Sal = 'nan'

		if len(l1) != 0:
			avg4Sal = sum(l4)/len(l4)
		else:
			avg4Sal = 'nan'

		salToAppend = [avg1Sal, avg2Sal, avg3Sal, avg4Sal]
		if i < 100:
			print(salToAppend)
		sal4.extend(salToAppend)

		tempToAppend = [sum(dayTempData[0:firstStop])/(firstStop-1), sum(dayTempData[firstStop:secondStop])/(secondStop - firstStop-1), 
		sum(dayTempData[secondStop:thirdStop])/(thirdStop - secondStop - 1), sum(dayTempData[thirdStop:])/(numMeas - thirdStop - 1)]
		temp4.extend(tempToAppend)

		del dayTempData[:]
		del daySalData[:]
		numMeas = 0
		nDate = np1Date


print(sal4[0:50])

