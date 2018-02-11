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

df = pandas.read_csv("OSMC_flattened_eff7_eca9_6785.csv", parse_dates=[0], header=0,delimiter=",")


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




def main():
	startDate = df['time'][0]

	plt.subplot(211)
	plt.scatter(getColData('ztmp')[1], getColData('ztmp')[0])

	plt.subplot(212)
	plt.scatter(getColData('zsal')[1], getColData('zsal')[0])
	plt.show()

if __name__ == "__main__":
	main()


