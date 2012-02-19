#!/usr/bin/env python
#
# compNeuroPython.py
# Created by Nicholas Cain on 9/14/11.

import os
import numpy as np
import pbsTools as pt
import pylab as pl

#-------------------------------------------------------------------------------

class NTL(object):

	def __init__(self, nameTimeListIn=None):

		# Null creator:
		if nameTimeListIn == None:
			data = []
		
		self.data = nameTimeListIn
		
#-------------------------------------------------------------------------------

class dirList(object):

	def __init__(self, dir="./"):
	
		
		self.fileList = os.listdir(dir)
		self.root = os.path.abspath(dir)
		
	def __iter__(self):
		return self.forward()

	def forward(self):
		current_item = 0
		while (current_item < len(self.fileList)):
			currfile = self.fileList[current_item]
			current_item += 1
			yield currfile

#-------------------------------------------------------------------------------

def getDataFromFile(fileName):

	# Load data from file:
	f1 = open(fileName)
	data = getDataFromStream(f1)
	f1.close()
	
	# Parse data:
	prefix, suffix = fileName.rsplit('.',1)
	spikeList = formatList(data, inputFormat=suffix)
	
	return spikeList

#-------------------------------------------------------------------------------

def getDataFromStream(inputPipe):

	# Read and parse data:
	data = [line for line in inputPipe.read().split('\n')][0:-1]
	inputPipe.close()
	
	return data
	
#-------------------------------------------------------------------------------

def formatList(dataList, inputFormat = 'ntf'):

	if inputFormat == 'ntf':
		
		# "ntf" stands for "Neuron-Time Format"
		delimiter = '\t'
		
	elif inputFormat == 'bsd':
		
		# "bsd" stands for "Brian Spike Data"
		delimiter = ','
		
	return NTL([[line.split(delimiter)[0].strip(),float(line.split(delimiter)[1].strip())] for line in dataList])
	
#-------------------------------------------------------------------------------
	
def printSpike(neuronTimeList, printFormat = 'ntf'):

	if printFormat == 'ntf':
		
		# "ntf" stands for "Neuron-Time Format"
		printString = neuronTimeList[0] + '\t' + str(neuronTimeList[1])
		
	elif printFormat == 'bsd':
		
		# "bsd" stands for "Brian Spike Data"
		printString = neuronTimeList[0] + ', ' + str(neuronTimeList[1])
		
	return printString
	
#-------------------------------------------------------------------------------

def getNamesFromList(spikeList):

	# Create a name list:
	names, times = zip(*spikeList)

	return list(set(names))
	
#-------------------------------------------------------------------------------

def firingRate(inputData, 
			   t=None,
			   window=50, 
			   nameList='All', 
			   plotsOn=False, 
			   dt=1,
			   xlabel='t',
			   ylabel='Firing Rate'):
		
	if isinstance(inputData, NTL):
		spikeList = inputData
	else:
		spikeList = NTL(inputData)
	
	# Import numpy tools:
	from numpy import arange, zeros, histogram, floor

	# Create the data lists, and restrict to valid names/times:
	if nameList == 'All':
		names, times = zip(*spikeList.data)
		nameList = list(set(names))
	popSize = len(nameList)
	validTimesList = []
	for neuronName, spikeTime in spikeList.data:
		if neuronName in nameList:
			validTimesList.append(spikeTime)
	
	# Create plotting axes:
	tMin = validTimesList[0]
	tMax = validTimesList[-1]
	if t==None:
		t = arange(tMin,tMax,dt)
	y = zeros(len(t))
	
	# Create histogram:
	events, edges = histogram(validTimesList, bins=t)
	maxInd = len(events)-1
	
	# Create firing rate function:
	def getFR(t):
	
		if t < tMin + window/2:
			lInd = 0
		else:
			lInd = floor((t-window/2)/dt)
		if t > tMax - window/2:
			rInd = maxInd
		else:
			rInd = floor((t+window/2)/dt)
			
		NSpikes = sum(events[lInd:rInd])

		return float(NSpikes)/(rInd-lInd)/popSize/dt/.001
		
	# Compute FR at each point in time domain:
	for i in range(len(t)):
		y[i] = getFR(t[i])
		
	if plotsOn:
		
		# Import plotting package:
		import pylab as pl
		
		# Make plots:
		pl.plot(t,y)
		pl.xlabel(xlabel)
		pl.ylabel(ylabel)
		
	return t,y
	
#-------------------------------------------------------------------------------

def firstCrossingSpikes(spikeList1, spikeList2, thetaList,
			   window=50, 
			   nameList='All', 
			   dt=1):
	
	# Get Firing Rate Data:
	ytVals1 = firingRate(spikeList1, window=window, nameList=nameList, plotsOn=False,
						dt=dt, xlabel='t', ylabel='Firing Rate')
	ytVals2 = firingRate(spikeList2, t=ytVals1[0], window=window, nameList=nameList, plotsOn=False,
						dt=dt, xlabel='t', ylabel='Firing Rate')
	
	return firstCrossing(ytVals1, ytVals2, thetaList)
	
#-------------------------------------------------------------------------------

def firstCrossing(ytVals1, ytVals2, thetaList):
			   
	if not isinstance(thetaList,list):
		thetaList = [thetaList]
	thetaList.sort()
	
	# Get Firing Rate Data:
	tVals, yVals1 = ytVals1
	tVals, yVals2 = ytVals2
	
	# Walk the time indices, and report first crossing:
	thetaInd = 0
	resultList = [-1]*len(thetaList)
	crossTimeList = [float("inf")]*len(thetaList)
	for i in range(len(tVals)):
		if yVals1[i] > thetaList[thetaInd] or yVals2[i] > thetaList[thetaInd]:
			if yVals1[i] > yVals2[i]:
				resultList[thetaInd]=1
			else:
				resultList[thetaInd]=0
			crossTimeList[thetaInd]=tVals[i]
			if thetaInd < len(thetaList)-1:
				thetaInd += 1
			else:
				break
				

	
	return crossTimeList, resultList
	
			
#-------------------------------------------------------------------------------

#def loadDir(dir='./'):
#	
#	# Walk the load directory:
#	run = []
#	UUIDHash = {}
#	for root, dirs, files in os.walk(dir):
#		for f in files:
#		
#			fullFileName = os.path.join(os.path.abspath(dir),f)
#		
#			# Extract data from file name:
#			fileNamePrefix,fileNameSuffix = f.rsplit('.',1)
#			whichPool, UUIDAndSettings = fileNamePrefix.split('_',1)
#			UUIDAndSettingsList = UUIDAndSettings.split('_')
#			UUID = UUIDAndSettingsList[0]
#			settingsList = UUIDAndSettingsList[1:]
#			
#			# Load data:
#			spikeData = getDataFromFile(fullFileName)
#			
#			# Add
#			if not UUID in UUIDHash.keys():
#				UUIDHash[UUID]=len(run)
#				run.append({'UUID':UUID, 
#							'Coh':settingsList[0],
#							'TOn':settingsList[1],
#							'TOff':settingsList[2],
#							'TMax':settingsList[3],
#							'inputCorr':settingsList[4]})
#				run[-1]['poolData']={whichPool:{'spikes':spikeData.data,
#												'FR':firingRate(spikeData)}}
#			else:
#				thisUUIDInd = UUIDHash[UUID]
#				run[thisUUIDInd]['poolData'][whichPool] = {'spikes':spikeData.data,
#														   'FR':firingRate(spikeData)}
#
#	return run

#-------------------------------------------------------------------------------


	
#-------------------------------------------------------------------------------

#def psychoChrono(data, saveResults=1):
#
#	# Import packages:
#	import numpy as np
#	
#	# Get Coherence Values:
#	CVals = []
#	for trial in data:
#		if not float(trial['Coh']) in CVals:
#			CVals.append(float(trial['Coh']))
#			
#	# Get theta values:
#	thetaVals = data[0]['theta']
#	
#	# Sort and make a hash table:
#	CVals.sort()
#	CValsHash = {}
#	for i in range(len(CVals)):
#		CValsHash[CVals[i]]=i
#		
#	# Compute number of samples:
#	N = len(data)/len(CVals)
#			
#	# Organize each trial by Coh type:
#	FCRTDict = {}
#	for i in range(len(thetaVals)):
#		RTVals = np.zeros(len(CVals))
#		FCVals = np.zeros(len(CVals))
#		
#		for trial in data:
#			RTVals[CValsHash[float(trial['Coh'])]] += float(trial['RT'][i])/N
#			FCVals[CValsHash[float(trial['Coh'])]] += float(trial['FC'][i])/N
#		
#		FCRTDict[thetaVals[i]] = {'FC':FCVals,'RT':RTVals}
#		
#	if saveResults:
#		
#		import pbsTools as pt
#		pt.pickle(FCRTDict, "../psychoChronoAnalysis.dat")
#	
#	return FCRTDict

#-------------------------------------------------------------------------------

def splitFileName(fileName):

	fileNamePartList = fileName.rsplit('.',1)
	if len(fileNamePartList) == 1:
		fileNamePrefix = fileNamePartList[0]
		fileNameSuffix = -1
	else:
		fileNamePrefix, fileNameSuffix = fileName.rsplit('.',1)
	
	return fileNamePrefix, fileNameSuffix

#-------------------------------------------------------------------------------

def ntfToFRFile(fileName):

	# Parse name:
	fileNamePrefix, fileNameSuffix = splitFileName(fileName)
	if not fileNameSuffix == "ntf":
	
		print "Not spike file."
		return 0
	# Get list of files in directory:
	fileList = dirList("./").fileList
	if not ((fileNamePrefix + ".fr") in fileList):
		# Get data:
		data = getDataFromFile(fileName)
	
		# Make into FR function:
		t,y = firingRate(data)

		# Write to fr file:
		doubleListToFile(t, y, fileNamePrefix + ".fr")
		return
	else:
		print "File already in directory"
		return 0
		
#-------------------------------------------------------------------------------
	
def doubleListToFile(L1, L2, fileName):
	
	f = open(fileName, 'w')
	try:	
		for i in range(len(L1)):
			f.write(str(L1[i]) + '\t' + str(L2[i]) + '\n')
			
		f.close()
		return 0
	except:
		f.close()
		os.remove(f)
		return -1
		
#-------------------------------------------------------------------------------
		
def doubleListFromFile(fileName):
	
	f = open(fileName, 'r')
	L1 = []
	L2 = []
	try:	
		for line in f:
			if line == "\n":
				break
			else:
				a,b = line.split("\t")
				L1.append(a)
				L2.append(b)
			
		return L1, L2
	except:
		f.close()
		return -1
		
#-------------------------------------------------------------------------------
		
def tripleListToFile(L1, L2, L3, fileName):
	
	f = open(fileName, 'w')
	try:	
		for i in range(len(L1)):
			f.write(str(L1[i]) + '\t' + str(L2[i]) + '\t' + str(L3[i]) + '\n')
			
		f.close()
		return 0
	except:
		f.close()
		os.remove(f)
		return -1
	
#-------------------------------------------------------------------------------

def tripleListFromFile(fileName):
	
	f = open(fileName, 'r')
	L1 = []
	L2 = []
	L3 = []
	try:	
		for line in f:
			if line == "\n":
				break
			else:
				a,b,c = line.split("\t")
				L1.append(a)
				L2.append(b)
				L3.append(c)
			
		return L1, L2, L3
	except:
		f.close()
		return -1
		
#-------------------------------------------------------------------------------
	
def applyFxnToDir(fxn, dir="./", verbose=True):
	
	for f in dirList(dir):
		fxn(f)
		
		if verbose:
			print "File " + f + " done."
	return
	
#-------------------------------------------------------------------------------

def ntfToFRDir(dir="./", verbose=True):
	applyFxnToDir(ntfToFRFile, dir=dir, verbose=verbose)
	return

#-------------------------------------------------------------------------------

def getUUID(fileName):

	if isType(fileName):
		fileNamePrefix, fileNameSuffix = splitFileName(fileName)
		whichPool, UUIDAndSettings = fileNamePrefix.split('_',1)
		UUIDAndSettingsList = UUIDAndSettings.split('_')
	else:
		UUIDAndSettingsList = [-1]
	
	return UUIDAndSettingsList[0]

#-------------------------------------------------------------------------------

def getUUIDList(dir="./"):

	UUIDList = []
	for f in dirList(dir):
		currUUID = getUUID(f)
		if (not currUUID in UUIDList) and (not currUUID == -1):
			UUIDList.append(currUUID)
			
	return UUIDList

#-------------------------------------------------------------------------------

def getSettings(fileName):

	fileNamePrefix, fileNameSuffix = splitFileName(fileName)
	whichPool, UUIDAndSettings = fileNamePrefix.split('_',1)
	UUIDAndSettingsList = UUIDAndSettings.split('_')

	return UUIDAndSettingsList[1:]


#-------------------------------------------------------------------------------

def getUUIDListCoh(Coh, dir = "./"):
	
	UUIDList = getUUIDList(dir=dir)
	UUIDListNew = []
	for UUID in UUIDList:
		exampleFileName = findFileName([UUID],1)[0]
		settings = getSettings(exampleFileName)
		if settings[0] == Coh:
			UUIDListNew.append(UUID)
			
	return UUIDListNew
	
#-------------------------------------------------------------------------------
	
def findFileName(substringList, N=100000000):

	validFileNames = []
	for f in dirList("./"):
		if all(map(f.count, substringList)):
			validFileNames.append(f)
			if len(validFileNames) >= N:
				break
			
	return validFileNames

#-------------------------------------------------------------------------------

#def thresholdTestCoh(Coh, thetaList):
#	UUIDList = getUUIDListCoh(Coh)
#	RTList = []
#	FCList = []
#	for UUID in UUIDList:
#		RT, FC = thresholdTestUUID(UUID, thetaList)
#		RTList.append(RT)
#		FCList.append(FC)
#		
#	return RTList, FCList

#-------------------------------------------------------------------------------

def meanRTFCCoh(Coh, thetaList):
	
	RTListMean = np.zeros(len(thetaList))
	FCListMean = np.zeros(len(thetaList))
	
	RTList, FCList = thresholdTestCoh(Coh, thetaList)
	
	for RT in RTList:
		RT = [float(val) for val in RT]
		RTListMean += np.array(RT)
	for FC in FCList:
		FC = [float(val) for val in FC]
		FCListMean += np.array(FC)
		
	return RTListMean/float(len(RTList)), FCListMean/float(len(FCList))
	
#-------------------------------------------------------------------------------
	
def psychoChrono(thetaList, saveResults=1, verbose=1):

	CohListStr = ["0.0","3.2","6.4","12.8","25.6","51.2"]
	CohList = [float(Coh) for Coh in CohListStr]
	
	UUIDList = getUUIDList(dir="./")
	
	thetaCohDict = {}
	for theta in thetaList:
		for Coh in CohList:
			thetaCohDict[theta, Coh] = [0,0,0]
	for UUID in UUIDList:
		if verbose == 1:
			print "Processing UUID " + UUID
		f = findFileName(["thresholdTest",UUID], N=1)[0]
		fSetting = findFileName([".ntf",UUID], N=1)[0]
		currCoh = float(getSettings(fSetting)[0])
		thetaListIn, RTList, FCList = tripleListFromFile(f)
		counter = 0
		for theta in thetaListIn:
			thetaCohDict[float(theta), currCoh][0] += float(RTList[counter])
			thetaCohDict[float(theta), currCoh][1] += float(FCList[counter])
			thetaCohDict[float(theta), currCoh][2] += 1
			counter += 1
		
		
	psyChrDict = {}
	for theta in thetaList:
		psyChrDict[theta] = [CohList,[],[]]
	for Coh in CohList:
		for theta in thetaList:
			psyChrDict[theta][1].append(thetaCohDict[theta, Coh][0]/thetaCohDict[theta, Coh][2])
			psyChrDict[theta][2].append(thetaCohDict[theta, Coh][1]/thetaCohDict[theta, Coh][2])
			
	print "total sims: " + str(thetaCohDict[theta, Coh][2])

	pt.pickle(psyChrDict,"../psychoChronoAnalysis.dat")

	
	return psyChrDict
	
#-------------------------------------------------------------------------------

def thresholdTestDir(thetaList, verbose=1):
	
	UUIDList = getUUIDList(dir="./")
	
	for UUID in UUIDList:
		thresholdTestUUID(UUID, thetaList, verbose=verbose)
	return
		
		
#-------------------------------------------------------------------------------
	
def thresholdTestUUID(UUID, thetaList, verbose=1):
	saveString = "thresholdTest_" + UUID
	for theta in thetaList:
		saveString += "_" + str(theta)
	fileName = saveString + ".dat"
	
	if os.path.isfile(fileName):
		thetaList, RTList, FCList = tripleListFromFile(fileName)
		if verbose:
			print "UUID " + UUID + " loaded."
	else:
	
		if verbose:
			print "Testing UUID: " + UUID
        
		try:
			GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
			GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
		except IndexError:
			ntfToFRFile(findFileName([UUID, ".ntf", "GESel1"])[0])
			GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
                
			ntfToFRFile(findFileName([UUID, ".ntf", "GESel2"])[0])
			GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
    
						
		t1,y1 = doubleListFromFile(GESel1FileName)
		t2,y2 = doubleListFromFile(GESel2FileName)

		t1 = [float(val) for val in t1]
		t2 = [float(val) for val in t2]
		y1 = [float(val) for val in y1]
		y2 = [float(val) for val in y2]
		
		RTList, FCList = firstCrossing([t1,y1], [t2,y2], thetaList)
		
		tripleListToFile(thetaList, RTList, FCList, fileName)
		if verbose:
			print "UUID " + UUID + " tested and saved."
	return RTList, FCList

#-------------------------------------------------------------------------------

def getRoitmanData(subject='n'):

	CohVals = [0, 3.2, 6.4, 12.8, 25.6, 51.2]
	if subject == 'both':
		CohValsB, RTValsB, FCValsB, CohDictB = getRoitmanData(subject='b')
		CohValsN, RTValsN, FCValsN, CohDictN = getRoitmanData(subject='n')
		return CohValsB[0], [RTValsB[0], RTValsN[0]], [FCValsB[0], FCValsN[0]], [CohDictB[0], CohDictN[0]]
		
	elif subject == 'b':
		FCVals = [0.50462962962963, 0.615560640732265, 0.738532110091743, 0.93348623853211, 0.995412844036697, 1	]
		RTVals = [787.601851851852, 776.871853546911, 738.5, 669.220183486239, 559.967889908257, 464.413242009132]	
	elif subject == 'n':
		FCVals = [0.495741056218058, 0.661590524534687, 0.804753820033956, 0.947189097103918, 0.994915254237288, 1]		
		RTVals = [853.93867120954, 851.991539763113, 801.504244482173, 694.926746166951, 529.932203389831, 392.464406779661]

	else:
		print "unrecognized subject"
		return -1
	
	counter = 0
	CohDict = {}
	for CohVal in CohVals:
		CohDict[CohVal] = [RTVals[counter], FCVals[counter]]
		counter += 1
		
	return [CohVals], [RTVals], [FCVals], [CohDict]

#-------------------------------------------------------------------------------

def psychoChronoAnalysis(data, theta, plotStyle='-', dots=True, subject = 'n'):

	x = data[theta][0]
	RT = data[theta][1]
	FC = data[theta][2]
	
	if dots == True:
		CohVals, RTVals, FCVals, CohDict = getRoitmanData(subject=subject)
	
	pl.figure(1)
	pl.plot(x, RT, plotStyle)
	if dots == True:
		if subject == 'both':
			pl.plot(CohVals[0], RTVals[0],'bo')
			pl.plot(CohVals[1], RTVals[1],'ro')
		elif subject == 'b':
			pl.plot(CohVals[0], RTVals[0],'bo')
		elif subject == 'n':
			pl.plot(CohVals[0], RTVals[0],'ro')
	
	pl.figure(2)
	pl.plot(x, FC, plotStyle)
	if dots == True:
		if subject == 'both':
			pl.plot(CohVals[0], FCVals[0],'bo')
			pl.plot(CohVals[1], FCVals[1],'ro')
		elif subject == 'b':
			pl.plot(CohVals[0], FCVals[0],'bo')
		elif subject == 'n':
			pl.plot(CohVals[0], FCVals[0],'ro')
	
#-------------------------------------------------------------------------------
	
def speedAccTradeoff(data, C, plotStyle='-', dots=True, highlightTheta=None, subject = 'n'):

	thetaVals = data.keys()
	CVals = data[thetaVals[0]][0]
	CInd = CVals.index(C)
	RT = [0]*len(thetaVals)
	FC = [0]*len(thetaVals)
	
	if dots == True:
		CohVals, RTVals, FCVals, CohDict = getRoitmanData(subject=subject)
		
	if not highlightTheta == None:
		specialThetaRT = data[highlightTheta][1][CInd]
		specialThetaFC = data[highlightTheta][2][CInd]
	
	counter = 0
	for theta in thetaVals:
		RT[counter] = 350 + data[theta][1][CInd]
		FC[counter] = data[theta][2][CInd]
		counter += 1
	
	pl.figure(3)
	pl.plot(RT, FC, plotStyle)
	if dots == True:
		if subject == 'both':
			pl.plot(CohDict[0][C][0], CohDict[0][C][1],'bo')
			pl.plot(CohDict[1][C][0], CohDict[1][C][1],'ro')
		elif subject == 'b':
			pl.plot(CohDict[0][C][0], CohDict[0][C][1],'bo')
		elif subject == 'n':
			pl.plot(CohDict[0][C][0], CohDict[0][C][1],'ro')
	if not highlightTheta == None:
		pl.plot([specialThetaRT], [specialThetaFC],'b.')
	pl.ylim([.49,1.01])
	
#-------------------------------------------------------------------------------
		
def speedAccTradeoffFull(data, plotStyle='-', dots=True, highlightTheta=None, subject = 'n'):

	CohVals = [0, 3.2, 6.4, 12.8, 25.6, 51.2]
	
	for CohVal in CohVals:
		speedAccTradeoff(data, CohVal, plotStyle=plotStyle, dots=dots, highlightTheta=highlightTheta, subject=subject)
	
#-------------------------------------------------------------------------------
	
def reliability(thetaList, saveResults=1):

	
	UUIDList = getUUIDList(dir="./")
	
	thetaDict = {}
	for theta in thetaList:
		thetaDict[theta] = [0,0,0]
	for UUID in UUIDList:
		print "Processing UUID " + UUID
		f = findFileName(["thresholdTest",UUID], N=1)[0]
		thetaListIn, RTList, FCList = tripleListFromFile(f)
		counter = 0
		for theta in thetaListIn:
			thetaDict[float(theta)][0] += float(RTList[counter])
			thetaDict[float(theta)][1] += float(FCList[counter])
			thetaDict[float(theta)][2] += 1
			counter += 1
		
		
	for theta in thetaList:
		thetaDict[float(theta)][0] = thetaDict[theta][0]/thetaDict[theta][2]
		thetaDict[float(theta)][1] = thetaDict[theta][1]/thetaDict[theta][2]
		thetaDict[float(theta)].pop(2)
			
	pt.pickle(thetaDict,"../thetaAnalysis.dat")

	
	return thetaDict
	
#-------------------------------------------------------------------------------

def isType(filename, types=['ntf','dat','fr']):

	if isinstance(types,str):
		types = [types]
	
	fileNamePrefix, fileNameSuffix = splitFileName(filename)
	
	if fileNameSuffix in types:
		return True
	else:
		return False


#-------------------------------------------------------------------------------

def stripUUIDAndSettings(filename):

	if isType(filename):

		prefix = filename.split("_")
		if not len(prefix) == 1:
			suffix = filename.rsplit(".")[-1]
		
			newFileName = prefix[0] + '.' +  suffix
			
			if not os.path.exists(newFileName):
				os.rename(filename, newFileName)
	else:
		print "Disallowed extension; skipping..."
		
	return 0
	
#-------------------------------------------------------------------------------

def stripUUIDAndSettingsDir(dir="./", verbose=True):

	applyFxnToDir(stripUUIDAndSettings, dir=dir, verbose=verbose)
	return	

#-------------------------------------------------------------------------------

def plotFR(filename,plotstyle = '-', figure=1):

	if isType(filename,'fr'):
		t,y = doubleListFromFile(filename)
		pl.figure(figure)
		pl.plot(t,y,plotstyle)
	
	return 0

#-------------------------------------------------------------------------------

def FDReproduceTestUUID(UUID):

	fileName = "FDReproduceTest_" + UUID + ".dat"
	
	if os.path.isfile(fileName):
		aCorrect, bCorrect = doubleListFromFile(fileName)
		aCorrect = aCorrect[0]
		bCorrect = bCorrect[0]

		print "UUID " + UUID + " loaded."
	else:
	
		GESel1aFileName = findFileName([UUID, ".fr", "GESel1a"])[0]
		GESel2aFileName = findFileName([UUID, ".fr", "GESel2a"])[0]
		GESel1bFileName = findFileName([UUID, ".fr", "GESel1b"])[0]
		GESel2bFileName = findFileName([UUID, ".fr", "GESel2b"])[0]
						
		t1a,y1a = doubleListFromFile(GESel1aFileName)
		t2a,y2a = doubleListFromFile(GESel2aFileName)
		t1b,y1b = doubleListFromFile(GESel1bFileName)
		t2b,y2b = doubleListFromFile(GESel2bFileName)

		y1aEnd = float(y1a[-1])
		y2aEnd = float(y2a[-1])
		y1bEnd = float(y1b[-1])
		y2bEnd = float(y2b[-1])
		
		
		if y1aEnd > y2aEnd:
			aCorrect = 1
		else:
			aCorrect = 0
			
		if y1bEnd > y2bEnd:
			bCorrect = 1
		else:
			bCorrect = 0
		
		doubleListToFile([aCorrect], [bCorrect], fileName)
		
		print "UUID " + UUID + " tested and saved."
	return aCorrect, bCorrect


#-------------------------------------------------------------------------------

def FDReproduceTestDir():
	
	UUIDList = getUUIDList(dir="./")
	
	for UUID in UUIDList:
		FDReproduceTestUUID(UUID)
	return
		
#-------------------------------------------------------------------------------
		
def FDReproduceAnalysis():

	UUIDList = getUUIDList(dir="./")
	
	resultList = []
	for UUID in UUIDList:
		print "Processing UUID " + UUID
		f = findFileName(["FDReproduceTest",UUID], N=1)[0]
		aCorrect, bCorrect = doubleListFromFile(f)
		resultList.append([float(aCorrect[0].strip()), float(bCorrect[0].strip())])

	resultData = np.zeros((2,len(resultList)))
	for i in range(len(resultList)):
		resultData[0,i] = resultList[i][0]
		resultData[1,i] = resultList[i][1]
	
	a,b = zip(*resultList)
	doubleListToFile(a, b, 'FDReproduceAnalysis.dat')
		
#-------------------------------------------------------------------------------
	
def singleListToFile(L1, fileName):
	
	f = open(fileName, 'w')
	try:	
		for i in range(len(L1)):
			f.write(str(L1[i]) + '\n')
			
		f.close()
		return 0
	except:
		f.close()
		os.remove(f)
		return -1
		
#-------------------------------------------------------------------------------
		
def singleListFromFile(fileName):
	
	f = open(fileName, 'r')
	L1 = []
	try:	
		for line in f:
			if line == "\n":
				break
			else:
				L1.append(line.strip())

		return L1
	except:
		f.close()
		return -1

#-------------------------------------------------------------------------------
		
def rasterPlot(fileName):
	
	spikes = getDataFromFile(fileName).data
	
	namesListPrime, times = zip(*spikes)
	namesList = list(set(namesListPrime))
	nameTimeDict = {}
	for name in namesList:
		nameTimeDict[name]=[]
		
	for n, t in spikes:
		nameTimeDict[n].append(t)
		
	pl.figure(4)
	for i in range(len(namesList)):
		pl.plot(nameTimeDict[namesList[i]],[i]*len(nameTimeDict[namesList[i]]),'b.')
	
	return

#-------------------------------------------------------------------------------

def speedAccTradeoffSpikeCountUnCorr(C, N=240, r0 = 40, bP = .4, bN = .4, dots = True, thetaVals=[0,15], plotStyle='-', subject='both'):

    rP = N*(r0 + bP*C)/1000
    rN = N*(r0 - bN*C)/1000

    if len(thetaVals) == 2:
        thetaVals = np.linspace(thetaVals[0],thetaVals[1],100)
    else:
        thetaVals = np.array(thetaVals)
    
    def FC(theta):
        return 1/(1+(rN/rP)**theta)
    
    def RT(theta):
        return theta/(rP-rN)*(1-(rN/rP)**theta)/(1+(rN/rP)**theta)
    
    RTVals = RT(thetaVals)
    FCVals = FC(thetaVals)

    pl.plot(RTVals, FCVals, plotStyle)
    if dots == True:
        CohVals, RTValsRoitman, FCValsRoitman, CohDict = getRoitmanData(subject=subject)
        
        if dots == True:
            if subject == 'both':
                pl.plot(CohDict[0][C][0], CohDict[0][C][1],'bo')
                pl.plot(CohDict[1][C][0], CohDict[1][C][1],'ro')
            elif subject == 'b':
                pl.plot(CohDict[0][C][0], CohDict[0][C][1],'bo')
            elif subject == 'n':
                pl.plot(CohDict[0][C][0], CohDict[0][C][1],'ro')
        pl.ylim([.49,1.01])

    return


















