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
def output(A, ls, freq): 
	return getN(getR(getPhi(A,ls,freqToLs(freq))))

def main():
	print(output(10,15,10e6))
	print(freqToLs(3000))
	print(getPhi(10, 10, 1e5))
	print(getR(getPhi(10,10,freqToLs(10e6))))

if __name__ == "__main__":
	main()
