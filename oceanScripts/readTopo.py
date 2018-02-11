import os
import clawpack.geoclaw.topotools as topotools
import matplotlib.pyplot as plt
#from numpy import genfromtxt
import numpy
import pandas
import csv
import math
from functools import reduce
from datetime import datetime
from datetime import timedelta
product = reduce((lambda x, y: x * y), [1, 2, 3, 4])
#from netCDF4 import Dataset

df = pandas.read_csv("OSMC_flattened_eff7_eca9_6785.csv", parse_dates=[0], header=0,delimiter=",") #data for 32.3N - 144.6E


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
		times4.append(datetime.strptime(day+"T01:00:00Z", "%Y-%m-%dT%H:%M:%SZ"))
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
	l = 500 #nm
	if len(salList) == len(tempList):
		for i in range(len(salList)):
			t = float(tempList[i])
			s = float(salList[i])
			if (not numpy.isnan(t)) and (not numpy.isnan(s)):
				rIndices.append(1.3247 - 2.5e-6*t**2 + s*(2e-4 - 8e-7*t) + 3300*l**(-2) - 3.2e7*l**(-4))
			else:
				rIndices.append(numpy.nan)
	return rIndices



def main():
	startDate = df['time'][0]

	salData = getColData('zsal')[0]
	tempData = getColData('ztmp')[0]
	times = getColData('zsal')[1]

	# plt.subplot(311)
	# plt.scatter(times, tempData)

	# plt.subplot(312)
	# plt.scatter(times, salData)

	# plt.subplot(313)
	plt.scatter(times, refractiveIndexList(salData, tempData))
	plt.show()

if __name__ == "__main__":
	main()


