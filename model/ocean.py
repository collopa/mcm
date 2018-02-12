import math

def freqToLs(freq):
	return 3.0e8/freq

def getPhi(A, ls, lr):
	return 2*math.sqrt(A**2 + ls**2/16)/lr

def getR(phi):
	return (math.cos(phi/2))**2

def getN(R):
	return (1-math.sqrt(R))/(1+math.sqrt(R))


#gets pseudo index of refraction as a function of 
#amplitude (A), ocean wavelength (ls), and radiowave length (lr)
def ocean_index_of_refr(amp_ocean, lambda_ocean, freq_radio): 
	return getN(getR(getPhi(amp_ocean,lambda_ocean,freqToLs(freq_radio))))


