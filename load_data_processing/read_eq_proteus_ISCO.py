###
###
### daniel vogler
### davogler@ethz.ch
###
### run in folder with respective load data as:
### python read_eq_proteus_ISCO.py ./
###
###


import csv
import collections

import numpy as np
import math

import matplotlib
import matplotlib.pyplot as pl
import matplotlib.colors as mcolors
import matplotlib.cm as cm

import sys
import glob, os
import warnings
import re


warnings.simplefilter("error")

#######

####### definitions

expRadius = 0.061

##################

print "System input variables"
print str(sys.argv)
print 

# read in filepath from command line
filepath =  sys.argv[1]
# figure path
figurePath = "./"
#loadCurve =  sys.argv[2]

### Read in experimental data #########################

print "###########################################"
print
print "Read experimental data"
print filepath
print 
print "###########################################"
print "Processing mechanical experimental data"
print


# define measurement locations on sample
mechTestFiles = []

# find all files with given string
os.chdir(filepath)
variableName = "ethz_vd"
searchString = str(variableName+"*mechData*.txt")
print "\t Searching test in path %s \n" %filepath
#
for file in glob.glob(searchString):
	mechTestFiles.append(file)
	print "\t\t Found test file %s " %file
print "\n\t\t --> %d test files found for %s \n" %(len(mechTestFiles), filepath)
mechTestFiles = sorted(mechTestFiles)

# search string
expSearchString = 'Index'
expIndicesGlobal = [[] for _ in range(len(mechTestFiles))]
#
# initialize variables
listExpIndex = [[] for _ in range(len(mechTestFiles))]
listExpElement = [[] for _ in range(len(mechTestFiles))]
listExpMarker = [[] for _ in range(len(mechTestFiles))]
listExpTime = [[] for _ in range(len(mechTestFiles))]
listExpForceDiff = [[] for _ in range(len(mechTestFiles))]
listExpMachineDisplacement = [[] for _ in range(len(mechTestFiles))]
listExpSensorDisplacement = [[] for _ in range(len(mechTestFiles))]
listExpSigmaZ = [[] for _ in range(len(mechTestFiles))]
# initialize variables
expIndex = [[] for _ in range(len(mechTestFiles))]
expElement = [[] for _ in range(len(mechTestFiles))]
expMarker = [[] for _ in range(len(mechTestFiles))]
expTime = [[] for _ in range(len(mechTestFiles))]
expForceDiff = [[] for _ in range(len(mechTestFiles))]
expMachineDisplacement = [[] for _ in range(len(mechTestFiles))]
expSensorDisplacement = [[] for _ in range(len(mechTestFiles))]
expSigmaZ = [[] for _ in range(len(mechTestFiles))]
#
maxLoadIdx = []
maxSigmaZIdx = []
# 
for loadFileCounter in range(len(mechTestFiles)):
	print "\tOpening file %s" %mechTestFiles[loadFileCounter]
	fileToLoad = str(filepath+mechTestFiles[loadFileCounter])
	with open(fileToLoad) as file:
	        expFile = csv.reader(file, delimiter=';', skipinitialspace=True)
	        for line in expFile:
	                expData = list(expFile)
	        #
			# find search string occurences in list
	        expIndices = [l for l, s in enumerate(expData) if str(expSearchString) in s]
	expIndicesGlobal[loadFileCounter] = expIndices
	# seperate data sheets
	for l in range(expIndices[0]+2, expIndices[1]-3):
	    listExpIndex[loadFileCounter].append( float(expData[l][0]) )
	    listExpElement[loadFileCounter].append( float(expData[l][1]) )
	    listExpMarker[loadFileCounter].append( float(expData[l][2]) )
	    listExpTime[loadFileCounter].append( expData[l][3]) 
	    listExpForceDiff[loadFileCounter].append( float(expData[l][4]) )
	    listExpMachineDisplacement[loadFileCounter].append( float(expData[l][5]) )
	    listExpSensorDisplacement[loadFileCounter].append( float(expData[l][6]) )
	    listExpSigmaZ[loadFileCounter].append( float(expData[l][4])/(math.pi*expRadius**2) )


# convert to array and print max results
for loadFileCounter in range(len(mechTestFiles)):
	expIndex[loadFileCounter] = np.array(listExpIndex[loadFileCounter]).astype(np.float)
	expForceDiff[loadFileCounter] = np.array(listExpForceDiff[loadFileCounter]).astype(np.float)
	expMachineDisplacement[loadFileCounter] = np.array(listExpMachineDisplacement[loadFileCounter]).astype(np.float)
	expSensorDisplacement[loadFileCounter] = np.array(listExpSensorDisplacement[loadFileCounter]).astype(np.float)
	expSigmaZ[loadFileCounter] = expForceDiff[loadFileCounter]/(math.pi*expRadius**2)/1e3
	# find maximum load and stress indexes
	maxLoadIdx.append( expForceDiff[loadFileCounter].argmax(axis=0) )
	maxSigmaZIdx.append( expSigmaZ[loadFileCounter].argmax(axis=0) )
	#
	print "Maximum force difference during experiment:"
	print "%s" %mechTestFiles[loadFileCounter]
	print "%3.2f [kN]" %expForceDiff[loadFileCounter][maxLoadIdx[loadFileCounter]]
	print "%3.2f [MPa]" %expSigmaZ[loadFileCounter][maxSigmaZIdx[loadFileCounter]]
	print 

# formula for sigma tensile
# sigma_tensile = ( 2*PrimaryFailureLoad) / (pi * diameter * thickness) in [MPa]



print "###########################################"
print
print "Read experimental data"
print filepath
print 
print "###########################################"
print "Processing flow experimental data"
print

# define measurement locations on sample
flowTestFiles = []

# find all files with given string
os.chdir(filepath)
variableName = "ethz_vd"
searchString = str(variableName+"*flowData*.txt")
print "\t Searching test in path %s \n" %filepath
#
for file in glob.glob(searchString):
	flowTestFiles.append(file)
	print "\t\t Found test file %s " %file
print "\n\t\t --> %d test files found for %s \n" %(len(flowTestFiles), filepath)
flowTestFiles = sorted(flowTestFiles)

# search string
expSearchString = 'Index'
expIndicesGlobal = [[] for _ in range(len(flowTestFiles))]
#
listExpTime = [[] for _ in range(len(flowTestFiles))]
listExpPressure = [[] for _ in range(len(flowTestFiles))]
listExpFlowrate = [[] for _ in range(len(flowTestFiles))]
listExpVolume = [[] for _ in range(len(flowTestFiles))]
#
expTime = [[] for _ in range(len(flowTestFiles))]
expPressure = [[] for _ in range(len(flowTestFiles))]
expFlowrate = [[] for _ in range(len(flowTestFiles))]
expVolume = [[] for _ in range(len(flowTestFiles))]
expTimeHMS = [[0 for x in range(4)] for _ in range(len(flowTestFiles))]
#
maxFlowRateIdx = []
maxPressureIdx = []
#
expFileLines = []
# 
for loadFileCounter in range(len(flowTestFiles)):
	print "\tOpening file %s" %flowTestFiles[loadFileCounter]
	fileToLoad = str(filepath+flowTestFiles[loadFileCounter])
	expFile = []
	expFileLines = []
	with open(fileToLoad,'r') as f:
	    next(f)
	    expFile=csv.reader(f,delimiter='\t')
	    for lines in expFile:
	    	expFileLines.append(lines)
    # read out shape
	expFileLinesShape = np.shape(expFileLines)
	# seperate data sheets
	for l in range( expFileLinesShape[0] ):
		listExpTime[loadFileCounter].append( float(expFileLines[l][0]) )
		listExpPressure[loadFileCounter].append( float(expFileLines[l][1]) )
		listExpFlowrate[loadFileCounter].append( float(expFileLines[l][2]) )
		listExpVolume[loadFileCounter].append( float(expFileLines[l][3]) )



print "\n\n####################################################################################"

# convert to array and print max results
for loadFileCounter in range(len(flowTestFiles)):
	expTime[loadFileCounter] = ( np.array(listExpTime[loadFileCounter]).astype(np.float) - np.array(listExpTime[loadFileCounter][0]).astype(np.float) ) *1e4#/10
	expPressure[loadFileCounter] = np.array(listExpPressure[loadFileCounter]).astype(np.float)*1e2
	expFlowrate[loadFileCounter] = np.array(listExpFlowrate[loadFileCounter]).astype(np.float)
	expVolume[loadFileCounter] = np.array(listExpVolume[loadFileCounter]).astype(np.float)
	# find maximum load and stress indexes
	maxFlowRateIdx.append( expFlowrate[loadFileCounter].argmax(axis=0) )
	maxPressureIdx.append( expPressure[loadFileCounter].argmax(axis=0) )
	#
	print "Maximum flow rate and pressure during experiment:"
	print "%s" %flowTestFiles[loadFileCounter]
	print "%3.2f [ml/s]" %expFlowrate[loadFileCounter][maxFlowRateIdx[loadFileCounter]]
	print "%3.0f [kPa]" %expPressure[loadFileCounter][maxPressureIdx[loadFileCounter]]
	print 



print "##############################################################"
print "Plotting Figures"
print

# plots
# settings
markerSize = 10.0
lineStyle = 'none'
legendLocation = "upper left"
Color = ['b', 'r', 'm', 'g']
conversionFactor = 1000
fontSize = 15

# save path
substring =  filepath.find("ethz_vd", 8)
saveString = filepath[substring:substring+8]
saveString = str(figurePath+saveString+"_")
print "saving figures with save string %s" %saveString


normalize = mcolors.Normalize( vmin=0.0, vmax=len(mechTestFiles) )
colormap = cm.hot#cm.gnuplot

# plot - exp results vs mean of simulation
pl.figure()

for loadFileCounter in range(len(mechTestFiles)):
	pl.plot( expSensorDisplacement[loadFileCounter] - expSensorDisplacement[loadFileCounter][0] , expSigmaZ[loadFileCounter], color=colormap(normalize(loadFileCounter)), marker='o', label=loadFileCounter+1, markersize=markerSize, linestyle="-")
	#pl.plot([], [], color=colormap(normalize(k)), label=mechTestFiles[loadFileCounter] )
	#pl.plot( expSensorDisplacement[0:maxSigmaZIdx], expSigmaZ[0:maxSigmaZIdx], c='r', marker='o', label='Experiment', markersize=markerSize, linestyle=lineStyle)

pl.ylabel('SigmaZ [MPa]', fontsize = fontSize)
pl.xlabel('Displacement [mm]', fontsize = fontSize)
pl.legend(loc=legendLocation, numpoints = 1,title="Cycle [-]")
pl.grid(b=True, which='major', color='lightgrey', linestyle='-')
pl.tick_params(axis='both', which='major', labelsize=fontSize)
pl.savefig( str(saveString+"displacementSensor_vs_sigmaZ.pdf"), bbox_inches='tight' )
pl.savefig( str(saveString+"displacementSensor_vs_sigmaZ.png"), bbox_inches='tight' )



# plot - fluid pressure vs sigma z of simulation
pl.figure()

for loadFileCounter in range(len(mechTestFiles)):
	pl.plot( expSigmaZ[loadFileCounter] , expPressure[loadFileCounter], color=colormap(normalize(loadFileCounter)), marker='o', label=loadFileCounter+1, markersize=markerSize, linestyle="-")

pl.xlabel('Sigma Z [MPa]', fontsize = fontSize)
pl.ylabel('Fluid Pressure [kPa]', fontsize = fontSize)
pl.legend(loc=legendLocation, numpoints = 1,title="Cycle [-]")
pl.grid(b=True, which='major', color='lightgrey', linestyle='-')
pl.tick_params(axis='both', which='major', labelsize=fontSize)
pl.savefig( str(saveString+"fluidPressure_vs_sigmaZ.png"), bbox_inches='tight' )
pl.savefig( str(saveString+"fluidPressure_vs_sigmaZ.pdf"), bbox_inches='tight' )



# plot - fluid pressure vs sigma z of simulation
pl.figure()

for loadFileCounter in range(len(mechTestFiles)):
	pl.plot( expSensorDisplacement[loadFileCounter] - expSensorDisplacement[loadFileCounter][0] , expPressure[loadFileCounter], color=colormap(normalize(loadFileCounter)), marker='o', label=loadFileCounter+1, markersize=markerSize, linestyle="-")

pl.xlabel('Displacement [mm]', fontsize = fontSize)
pl.ylabel('Fluid Pressure [kPa]', fontsize = fontSize)
pl.legend(loc=legendLocation, numpoints = 1,title="Cycle [-]")
pl.grid(b=True, which='major', color='lightgrey', linestyle='-')
pl.tick_params(axis='both', which='major', labelsize=fontSize)
pl.savefig( str(saveString+"fluidPressure_vs_fractureClosure.png"), bbox_inches='tight' )
pl.savefig( str(saveString+"fluidPressure_vs_fractureClosure.pdf"), bbox_inches='tight' )


pl.show()
exit()
