# Launch the script: python mStar.py csvFile.csv inputs.txt
# Inputs to the script:
# 1- csv file (without header) and with the following columns
# 	|_time (s), 
# 	|_position of the piston (position(t0) = 0, position(end) = stroke) (m),
# 	|_pressure (Pa)
# 2- text file specifying the following (in this order)
# 	|_species name of the mixture (Name, in IUPAC form or common form or a synonym registered in PubChem, or CAS number)
# 	|_mass fractions for each species in the mixture
# 	|_time of beginning of compression (s), if None is specified, then it is detected automatically
# 	|_time of end of compression, if None is specied, this is by default the time of maximum pressure
# 	|_volume (m3) of the main chamber at TDC (without crevices but including other dead volumes)
# 	|_volume of the crevices (not including other dead voumes) (m3)
# 	|_radius of the combustion chamber (m)
# 	|_temperature at the beginning of compression (K)
# 	|_Stroke (m), if the stroke is obtained from the position vector, then indicate 0

# Output: The ratio between the volume of the crevices and the volume needed to absorb
# the boundary layer. The ratio should be bigger than 1 for the crevices to be effective

# Details of the implementation are provided in the following article:
# N. Bourgeois, H. Jeanmart, G. Winckelmans, O. Lamberts, F. Contino
# How to ensure the interpretability of experimental data in Rapid Compression Machines? A method to validate piston crevice designs
# Combustion and Flame, 2018

# Contributors to this script:
# N. Bourgeois, M. Pochet, F. Contino
#

import numpy as np
try:
    import pandas as pd
except ImportError:
    print("- Error: Module pandas not found.\n")
    print(" Please install it using your favourite library manager.\n")
    print(" Alternatively, you can run command: python -m pip install --upgrade pandas\n")

import scipy.integrate as integrate
from scipy.signal import savgol_filter
import math
try: 
    from thermo.chemical import Mixture
except ImportError:
    print("- Error: Module thermo not found.\n")
    print(" Please install it using your favourite library manager.\n")
    print(" Alternatively, you can run command: python -m pip install --upgrade thermo\n")

import sys
import warnings
import matplotlib.pyplot as plt

def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'

warnings.formatwarning = custom_formatwarning

# Default parameters
pmin = 10000 # Minimum pressure to check if in Pa
Tmin = 250 # Minimum temperature to check if in K
A = 1.2 # Scaling constant for vortex boundary layer, defined by Bourgeois et al.
thresholdInitialTime = 0.1
smoothFactor = 41


# Function to read the time, position and pressure from csv file
def readCSV(filename, t0, tEnd, stroke):
	df = pd.read_csv(filename, sep=None,header=None,usecols=[0,1,2], engine='python', names=['time', 'position', 'pressure'])
	
	pressure = df.get('pressure').values
	if np.min(pressure) < pmin:
		raise ValueError("\n\n The pressure is unlikely to be in Pa,\n if yes set pmin in your input file to a value below the minimum pressure.")

	position = df.get('position').values
	if position[0] > position[-1]:
		raise ValueError("\n\n The position vector should increase. The end position is smaller than the start.")		

	position -= position[0] # Make sure it starts at 0

	time = df.get('time').values
	if tEnd==None:
		print("\n\n The compression end time is determined automatically...")
		tEnd = time[np.argmax(pressure)]
		print("\n tEnd = %f" % (tEnd))
	if t0==None:
		print("\n\n The compression start time is determined automatically...\n")
		[position, time] = smooth(position, smoothFactor, time)
		speed = np.gradient(position, time)
		[speed, time] = smooth(speed, smoothFactor, time)
		accel = np.gradient(speed, time)
		[accel, time] = smooth(accel, smoothFactor, time)
		t0 = time[np.argmax(accel>thresholdInitialTime*np.max(accel))]
		print("\n t0 = %f" % (t0))
		

	if t0 > np.max(time) or t0 < np.min(time):
		raise ValueError("\n\n The initial time t0 is not within the specified time series.")
	criteria = (df['time'] >= t0) & (df['time'] <= tEnd)
	time = df.loc[criteria, 'time'].values  - t0
	
	pressure = df.loc[criteria, 'pressure'].values

	position = df.loc[criteria, 'position'].values
	position -= position[0]
	if stroke != (position[-1] - position[0]):
		warnings.warn("\n\n The difference of position from begin to end does not correspond to the stroke.\n The position has been scaled so that it corresponds.")
		if stroke > 0:
			position = position/(position[-1] - position[0])*stroke

	#plt.plot(time, position)

	return [time, position, pressure]

def smooth(y, box_pts, x):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    totalLength = len(y) - box_pts + 1
    endPoint = int((box_pts+1)/2 - 1) + totalLength
    x = x[int((box_pts+1)/2 - 1):endPoint]
    return [y_smooth, x]

# Function to read the file of inputs
def readInput(filename):
	with open(filename) as f:
		inputs = f.read().splitlines()

	species = inputs[0].split('#')[0].split(',')
	ws = [eval(i) for i in inputs[1].split(',')]
	t0 = eval(inputs[2])
	tEnd = eval(inputs[3])
	VClearance = eval(inputs[4])
	VCrev = eval(inputs[5])
	R = eval(inputs[6])
	T0 = eval(inputs[7])
	stroke = eval(inputs[8])
	return [species, ws, t0, tEnd, VClearance, VCrev, R, T0, stroke]


def main(csvFile,inputFile):
	[species, ws, t0, tEnd, VClearance, VCrev, R, T0, stroke] = readInput(inputFile)
	mixture = Mixture(species, ws=ws)
	RStar = mixture.R_specific # Ideal gas constant [J/(kg K)]

	try:
	    [time, position, pressure] = readCSV(csvFile, t0, tEnd, stroke)
	except ValueError as err:
	    print(err.message)

	p0 = pressure[0] # Initial pressure
	pTDC = pressure[-1] # Pressure at top dead center
	V0 = VClearance + VCrev + stroke*math.pi*R**2
	VTDC = VClearance + VCrev #V chamber without crevices
	H = VClearance/(math.pi*R**2) + stroke # Clearance height at BDC
	nPoly = np.log(pTDC/p0)/np.log((V0)/VTDC)
	VChamber = V0 - VCrev - position*math.pi*R**2
	TChamber = T0*(VChamber[0]/VChamber)**(nPoly - 1.0)


	nu = np.zeros(len(TChamber))
	rhoChamber = np.zeros(len(TChamber))
	for Ti in range(len(TChamber)):
		mixture.calculate(T=TChamber[Ti],P=pressure[Ti])
		rhoChamber[Ti] = pressure[Ti]/(RStar*TChamber[Ti])
		nu[Ti] = mixture.mug/rhoChamber[Ti]

	vPiston = np.gradient(position, time)

	mDotStar = np.zeros(len(time))
	for ti in range(2,len(time)):
		nut = integrate.simps(nu[0:ti], time[0:ti])
		deltaOmega = A*np.sqrt(4.0*nut)*np.sqrt((H - position[ti])/H)
		mDotStar[ti] = rhoChamber[ti]*math.pi*(R**2 - (R - deltaOmega)**2)*vPiston[ti]

	mStar = integrate.simps(mDotStar, time)
	mChamberTDC = rhoChamber[-1]*VChamber[-1]
	mTot = p0*V0/(RStar*T0)
	mCrevTDC = mTot - mChamberTDC
	TCrevTDC = pTDC*VCrev/(mCrevTDC*RStar)
	rhoCrevTDC = pTDC/(RStar*TCrevTDC)
	VStar = mStar/(rhoCrevTDC - rhoChamber[0])

	print("\n The ratio of the volume of the crevices to the volume needed to absorb the vortex is :")
	print(" %f" % (VCrev/VStar))
	if (VCrev/VStar) > 1.0:
		print("\n Your crevices seem to be big enough.")


main("csvFile.csv","inputs.txt")
