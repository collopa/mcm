import os
import clawpack.geoclaw.topotools as topotools
import matplotlib.pyplot as plt
import pandas
import csv
import math
from functools import reduce
from datetime import datetime
from datetime import timedelta
import time
import numpy as np
from numpy import pi, r_, array, savetxt
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit

df = pandas.read_csv("OSMC_flattened_6977_e409_1aff.csv", parse_dates=[0], header=0,delimiter=",") #data for 32.3N - 144.6E
#df = pandas.read_csv("OSMC_flattened_eff7_eca9_6785.csv", parse_dates=[0], header=0,delimiter=",")

def getColData(column):

	avgList = []
	times = []
	nDate = df['time'][0][0:10]
	dayData = []
	numMeas = 0

	for i in range(df['time'].count()):

		np1Date = df['time'][i][0:10]

		if nDate == np1Date:
			dayData.append(df[column][i])
			numMeas +=1 #increment number of measurements in each day
				
		else: #new day- divide into 4 chunks and average each
			secondStop = math.floor(numMeas/2)
			firstStop = math.floor(secondStop/2)
			thirdStop = math.floor(numMeas/2) + secondStop

			l1 = list(filter(lambda x: not math.isnan(x), dayData[0:firstStop]))
			l2 = list(filter(lambda x: not math.isnan(x), dayData[firstStop:secondStop]))
			l3 = list(filter(lambda x: not math.isnan(x), dayData[secondStop:thirdStop]))
			l4 = list(filter(lambda x: not math.isnan(x), dayData[thirdStop:]))

			if len(l1) != 0:
				avg1 = sum(l1)/len(l1)
			else:
				avg1 = 'nan'

			if len(l2) != 0:
				avg2 = sum(l2)/len(l2)
			else:
				avg2 = 'nan'

			if len(l3) != 0:
				avg3 = sum(l3)/len(l3)
			else:
				avg3 = 'nan'

			if len(l4) != 0:
				avg4 = sum(l4)/len(l4)
			else:
				avg4 = 'nan'

			toAppend = [avg1, avg2, avg3, avg4]

			avgList.extend(toAppend)
			times.extend(get4Times(nDate))

			del dayData[:]
			numMeas = 0
			nDate = np1Date
	return(avgList, times)

def get4Times(day):
	times4 = []
	for i in range(4):
		date = datetime.strptime(day+"T01:00:00Z", "%Y-%m-%dT%H:%M:%SZ") + timedelta(hours = (2 + 6*i))
		times4.append(date.timestamp())
	return times4


'''def getTimes(startDate, length):
	datetimes = []
	currentDatetime = datetime.strptime(startDate, "%Y-%m-%dT%H:%M:%SZ")
	for i in range(length):
		datetimes.append(currentDatetime)
		currentDatetime += timedelta(hours = 3)
	return datetimes
	'''

def refractiveIndexList(salList, tempList):
	rIndices = []
	l = 10e10 #nm
	if len(salList) == len(tempList):
		for i in range(len(salList)):
			t = float(tempList[i])
			s = float(salList[i])
			if (not np.isnan(t)) and (not np.isnan(s)):
				rIndices.append(1.3247 - 2.5e-6*t**2 + s*(2e-4 - 8e-7*t) + 3300*l**(-2) - 3.2e7*l**(-4))
			else:
				rIndices.append(np.nan)
	return rIndices

def func(t, a, b, c, d, e, f):
	return (a*math.cos(b*t + c) + d*math.sin(e*t + f))

def clean(x,y):
	newX = []
	newY = []
	for i in range(len(x)):
		if not np.isnan(y[i]):
			newX.append(x[i])
			newY.append(y[i])
	return newX, newY


def main():
	startDate = df['time'][0]

	salData = getColData('zsal')[0]
	tempData = getColData('ztmp')[0]
	times = getColData('zsal')[1]
	rIndices = refractiveIndexList(salData, tempData)

	print(times[0])

	# Fit the first set
	# fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/p[1]*x+p[2]) + p[3]*x # Target function
	# errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
	# p0 = [-15., 0.8, 0., -1.] # Initial guess for the parameters
	# p1, success = optimize.leastsq(errfunc, p0[:], args=(Tx, tX))

	# time = np.linspace(Tx.min(), Tx.max(), 100)
	# plt.plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-")

	# plt.subplot(311)
	# plt.scatter(times, tempData)

	# plt.subplot(312)
	# plt.scatter(times, salData)

	#plt.subplot(313)
	#valid = ~(np.isnan(np.array(times)) | np.isnan(np.array(rIndices)))
	#print(times[valid])
	#x, y = clean(times, rIndices)
	#popt, pcov = curve_fit(func, x, y)
	#plt.plot(times, func(rIndices, *popt))

	with open('newData.csv','w') as f:
		writer = csv.writer(f, delimiter=",")
		writer.writerows(zip(times, rIndices))

	spl = UnivariateSpline(times, rIndices)
	xs = np.linspace(1478163600, 1508349600, 100000)
	plt.plot(xs, spl(xs), 'g', lw=3)
	plt.plot(array(times), array(rIndices))
	plt.ylim(1.34,1.345)
	plt.show()

if __name__ == "__main__":
	main()


