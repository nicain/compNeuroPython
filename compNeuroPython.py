#!/usr/bin/env python
#
# compNeuroPython.py
# Created by Nicholas Cain on 9/14/11.

import os
import numpy as np
import pbsTools as pt
import pylab as pl
import time
import random

def tic():
    return time.clock()

def toc(timeIn):
    print time.clock() - timeIn

#-------------------------------------------------------------------------------

class NTL(object):

    def __init__(self, nameTimeListIn=None):

        # Null creator:
        if nameTimeListIn == None:
            data = []
        
        self.data = nameTimeListIn

        self._whenSpiked = {}

    def numberOfNeurons(self):

        return len(set(zip(*self.data)[0]))

    def getNameList(self):
        
        return set(zip(*self.data)[0])

    def numberOfNeurons(self):
    
        return len(set(zip(*self.data)[0]))
        
    def getNeuronSpikeTimes(self, n):
        
        if not isinstance(n,str):
            n = str(n)
        
        if n in self._whenSpiked.keys():
            return self._whenSpiked[n]

        spikeTimeList = []

        for currN,currT in self.data:
            if n == currN:
                spikeTimeList.append(currT)

        self._whenSpiked[n] = spikeTimeList

        return spikeTimeList

        
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
               ylabel='Firing Rate',
               causal=True):
    
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
    if t==None:
        tMin = validTimesList[0]
        tMax = validTimesList[-1]
        t = arange(tMin,tMax,dt)
        events, edges = histogram(validTimesList, bins=t)
    else:
        tMin=min(validTimesList[0],t[0])
        tMax=max(validTimesList[-1],t[-1])
        dt=t[1]-t[0]
        histT =  arange(tMin-window,tMax,dt)
        events, edges = histogram(validTimesList, bins=histT)
    y = zeros(len(t))
    
    # Create firing rate function:
    maxInd = len(events)-1
    def getFR(t):
        if causal == False:
            if t < tMin + window/2:
                lInd = 0
            else:
                lInd = floor((t-window/2)/dt)
            if t > tMax - window/2:
                rInd = maxInd
            else:
                rInd = floor((t+window/2)/dt)
        else:
            if t < tMin + window:
                lInd = 0
            else:
                lInd = floor((t-window)/dt)
            if t > tMax - window:
                rInd = maxInd
            else:
                rInd = floor(t/dt)
            
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
               dt=1,
                        tOn=0):
    
    # Get Firing Rate Data:
    ytVals1 = firingRate(spikeList1, window=window, nameList=nameList, plotsOn=False,
                        dt=dt, xlabel='t', ylabel='Firing Rate')
    ytVals2 = firingRate(spikeList2, t=ytVals1[0], window=window, nameList=nameList, plotsOn=False,
                        dt=dt, xlabel='t', ylabel='Firing Rate')
    
    return firstCrossing(ytVals1, ytVals2, thetaList, tOn=tOn)
    
#-------------------------------------------------------------------------------

def firstCrossing(ytVals1, ytVals2, thetaList, tOn = 0):
               
    if not isinstance(thetaList,list):
        thetaList = [thetaList]
    thetaList.sort()
    
    # Get Firing Rate Data:
    t1, y1 = ytVals1
    t2, y2 = ytVals2



    tMin = min(min(t1),min(t2)) 
    tMax = max(max(t1),max(t2))
    tLen = max(len(t1), len(t2))
    tNew = np.linspace(tMin, tMax, tLen)
    
    y1New = np.zeros(len(tNew))
    y2New = np.zeros(len(tNew))
    for i in range(len(y1New)):
        currT = tNew[i]
        tmpInds = np.nonzero(t1<=currT)[0]
        if len(tmpInds) == 0:
            y1LInd = 0
        else:
            y1LInd = tmpInds[-1]
        tmpInds = np.nonzero(t1>=currT)[0]
        if len(tmpInds) == 0:
            y1RInd = len(t1)-1
        else:
            y1RInd = tmpInds[0]
        
        tmpInds = np.nonzero(t2<=currT)[0]
        if len(tmpInds) == 0:
            y2LInd = 0
        else:
            y2LInd = tmpInds[-1]
        tmpInds = np.nonzero(t2>=currT)[0]
        if len(tmpInds) == 0:
            y2RInd = len(t2)-1
        else:
            y2RInd = tmpInds[0]
        
        y1New[i] = (y1[y1RInd]+y1[y1LInd])/2            
        y2New[i] = (y2[y2RInd]+y2[y2LInd])/2



    ti = np.nonzero(np.array(tNew)<=tOn)[0][-1] + 1
    currTime = tNew[ti]
    FCList = []
    RTList = []
    for theta in thetaList:
        try:
            while y1New[ti] < theta and y2New[ti] < theta:
                ti += 1
                currTime = tNew[ti]
               
            print y1New[ti], y2New[ti]
            if y1New[ti] >= theta and y2New[ti] >= theta:
                if random.random() < .5:
                    FCList.append(1)
                else:
                    FCList.append(0)
            elif y1New[ti] >= theta:
                FCList.append(1)
            else:
                FCList.append(0)
            RTList.append(currTime)
        
        except:
            FCList.append(-1)
            RTList.append(float("inf"))

                

    
    return RTList, FCList
    
            
#-------------------------------------------------------------------------------

#def loadDir(dir='./'):
#    
#    # Walk the load directory:
#    run = []
#    UUIDHash = {}
#    for root, dirs, files in os.walk(dir):
#        for f in files:
#        
#            fullFileName = os.path.join(os.path.abspath(dir),f)
#        
#            # Extract data from file name:
#            fileNamePrefix,fileNameSuffix = f.rsplit('.',1)
#            whichPool, UUIDAndSettings = fileNamePrefix.split('_',1)
#            UUIDAndSettingsList = UUIDAndSettings.split('_')
#            UUID = UUIDAndSettingsList[0]
#            settingsList = UUIDAndSettingsList[1:]
#            
#            # Load data:
#            spikeData = getDataFromFile(fullFileName)
#            
#            # Add
#            if not UUID in UUIDHash.keys():
#                UUIDHash[UUID]=len(run)
#                run.append({'UUID':UUID, 
#                            'Coh':settingsList[0],
#                            'TOn':settingsList[1],
#                            'TOff':settingsList[2],
#                            'TMax':settingsList[3],
#                            'inputCorr':settingsList[4]})
#                run[-1]['poolData']={whichPool:{'spikes':spikeData.data,
#                                                'FR':firingRate(spikeData)}}
#            else:
#                thisUUIDInd = UUIDHash[UUID]
#                run[thisUUIDInd]['poolData'][whichPool] = {'spikes':spikeData.data,
#                                                           'FR':firingRate(spikeData)}
#
#    return run

#-------------------------------------------------------------------------------


    
#-------------------------------------------------------------------------------

#def psychoChrono(data, saveResults=1):
#
#    # Import packages:
#    import numpy as np
#    
#    # Get Coherence Values:
#    CVals = []
#    for trial in data:
#        if not float(trial['Coh']) in CVals:
#            CVals.append(float(trial['Coh']))
#            
#    # Get theta values:
#    thetaVals = data[0]['theta']
#    
#    # Sort and make a hash table:
#    CVals.sort()
#    CValsHash = {}
#    for i in range(len(CVals)):
#        CValsHash[CVals[i]]=i
#        
#    # Compute number of samples:
#    N = len(data)/len(CVals)
#            
#    # Organize each trial by Coh type:
#    FCRTDict = {}
#    for i in range(len(thetaVals)):
#        RTVals = np.zeros(len(CVals))
#        FCVals = np.zeros(len(CVals))
#        
#        for trial in data:
#            RTVals[CValsHash[float(trial['Coh'])]] += float(trial['RT'][i])/N
#            FCVals[CValsHash[float(trial['Coh'])]] += float(trial['FC'][i])/N
#        
#        FCRTDict[thetaVals[i]] = {'FC':FCVals,'RT':RTVals}
#        
#    if saveResults:
#        
#        import pbsTools as pt
#        pt.pickle(FCRTDict, "../psychoChronoAnalysis.dat")
#    
#    return FCRTDict

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

def ntfToFRDir(dir="./", verbose=True, tIn=None, causal=True, window=50):
    
    for f in dirList(dir):
        ntfToFRFile(f, window=window, causal=causal, tIn=tIn)

    return

#-------------------------------------------------------------------------------

def ntfToFRFile(fileName,window=50, causal=True, tIn=None):

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
        t,y = firingRate(data,window=window, causal=causal, t=tIn)

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
        
def doubleListFromFile(fileName, isFloat=False, whichDir='./', delim="\t"):
    
    f = open(os.path.join(whichDir,fileName), 'r')
    L1 = []
    L2 = []
    try:    
        for line in f:
            if line == "\n":
                break
            else:
                a,b = line.strip().split(delim)
                L1.append(a)
                L2.append(b)

        if isFloat==True:
            L1 = np.array(map(float,L1))
            L2 = np.array(map(float,L2))
            
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

def tripleListFromFile(fileName, isFloat=False):
    
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

        if isFloat==True:
            L1 = np.array(map(float,L1))
            L2 = np.array(map(float,L2))
            L3 = np.array(map(float,L3))
            
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
    
def findFileName(substringList, N=100000000,whichDir="./"):

    validFileNames = []
    for f in dirList(whichDir):
        if all(map(f.count, substringList)):
            validFileNames.append(os.path.join(whichDir,f))
            if len(validFileNames) >= N:
                break
            
    return validFileNames

#-------------------------------------------------------------------------------

#def thresholdTestCoh(Coh, thetaList):
#    UUIDList = getUUIDListCoh(Coh)
#    RTList = []
#    FCList = []
#    for UUID in UUIDList:
#        RT, FC = thresholdTestUUID(UUID, thetaList)
#        RTList.append(RT)
#        FCList.append(FC)
#        
#    return RTList, FCList

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
    
    UUIDcounter = 0
    for UUID in UUIDList:
        UUIDcounter = UUIDcounter + 1
        if verbose == 1:
            print "Processing UUID " + UUID + "  (" + str(UUIDcounter) + "/" + str(len(UUIDList)) + ")"
        f = findFileName(["thresholdTest",UUID], N=1)[0]
        fSetting = findFileName([".ntf",UUID], N=1,whichDir="../NTF")[0]
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

def thresholdTestDir(thetaList, verbose=1, tOn = 0):
    
    UUIDList = getUUIDList(dir="./")
    
    for UUID in UUIDList:
        thresholdTestUUID(UUID, thetaList, verbose=verbose, tOn=tOn)
    return
        
        
#-------------------------------------------------------------------------------
    
def thresholdTestUUID(UUID, thetaList, verbose=1, tOn = 0):
    
    fileName = "thresholdTestFR_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
    
        if verbose:
            print "Testing UUID (firing rate): " + UUID
        
        try:
            GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
            GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
        except IndexError:
            ntfToFRFile(findFileName([UUID, ".ntf", "GESel1"])[0])
            GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
                
            ntfToFRFile(findFileName([UUID, ".ntf", "GESel2"])[0])
            GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
    
                        
        t1,y1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2,y2 = doubleListFromFile(GESel2FileName, isFloat=True)       

        RTList, FCList = firstCrossing([t1,y1], [t2,y2], thetaList, tOn=tOn)
        
        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList

#-------------------------------------------------------------------------------

def thresholdTestUUIDDiff(UUID, thetaList, verbose=1, tOn = 0):
    
    fileName = "thresholdTestFRDiff_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (firing rate, diff): " + UUID
        
        try:
            GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
            GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
        except IndexError:
            ntfToFRFile(findFileName([UUID, ".ntf", "GESel1"])[0])
            GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
            
            ntfToFRFile(findFileName([UUID, ".ntf", "GESel2"])[0])
            GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
        
        
        t1,y1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2,y2 = doubleListFromFile(GESel2FileName, isFloat=True)
        
        tMin = min(min(t1),min(t2)) 
        tMax = max(max(t1),max(t2))
        tLen = max(len(t1), len(t2))
        t3 = np.linspace(tMin, tMax, tLen)
    
        y3 = np.zeros(len(t3))
        for i in range(len(y3)):
            currT = t3[i]
            tmpInds = np.nonzero(t1<=currT)[0]
            if len(tmpInds) == 0:
                y1LInd = 0
            else:
                y1LInd = tmpInds[-1]
            tmpInds = np.nonzero(t1>=currT)[0]
            if len(tmpInds) == 0:
                y1RInd = len(t1)-1
            else:
                y1RInd = tmpInds[0]

            tmpInds = np.nonzero(t2<=currT)[0]
            if len(tmpInds) == 0:
                y2LInd = 0
            else:
                y2LInd = tmpInds[-1]
            tmpInds = np.nonzero(t2>=currT)[0]
            if len(tmpInds) == 0:
                y2RInd = len(t2)-1
            else:
                y2RInd = tmpInds[0]
                
            y1Tmp = (y1[y1RInd]+y1[y1LInd])/2            
            y2Tmp = (y2[y2RInd]+y2[y2LInd])/2
                
            y3[i] = y1Tmp - y2Tmp


        y3i = np.nonzero(t3<=tOn)[0][-1]
        FCList = []
        RTList = []
        for theta in thetaList:
            
            try:

                while abs(y3[y3i]) < theta:
                    y3i += 1
                
                if y3[y3i] >= theta:
                    FCList.append(1)
                else:
                    FCList.append(0)
                RTList.append(t3[y3i])
            
            except:
                FCList.append(-1)
                RTList.append(float("inf"))

        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."

    return RTList, FCList

#-------------------------------------------------------------------------------

def thresholdTestICUUID(UUID, thetaList, verbose=1, tOn = 0, FRDelta=.2):
    
    fileName = "thresholdTestFRIC_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (firing rate, delta): " + UUID
        
        try:
            GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
            GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
        except IndexError:
            ntfToFRFile(findFileName([UUID, ".ntf", "GESel1"])[0])
            GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
            
            ntfToFRFile(findFileName([UUID, ".ntf", "GESel2"])[0])
            GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
        
        
        t1,y1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2,y2 = doubleListFromFile(GESel2FileName, isFloat=True)
        
        ti1 = np.nonzero(np.array(t1)<=tOn)[0][-1] + 1
        ti2 = np.nonzero(np.array(t1)<=tOn)[0][-1] + 1
    
        yDelta = abs(y1[ti1] - y2[ti2])
        if yDelta < FRDelta:
        
        
            RTList, FCList = firstCrossing([t1,y1], [t2,y2], thetaList, tOn=tOn)
            
            tripleListToFile(thetaList, RTList, FCList, fileName)
            if verbose:
                print "UUID " + UUID + " tested and saved."
            return RTList, FCList

        else:

            if verbose:
                print "UUID " + UUID + " Rejected (deltaFR)."
            return 0

#-------------------------------------------------------------------------------

def thresholdTestICUUIDDiff(UUID, thetaList, verbose=1, tOn = 0, FRDelta=.2):
    
    fileName = "thresholdTestFRICDiff_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (firing rate, diff): " + UUID
        
        try:
            GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
            GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
        except IndexError:
            ntfToFRFile(findFileName([UUID, ".ntf", "GESel1"])[0])
            GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
            
            ntfToFRFile(findFileName([UUID, ".ntf", "GESel2"])[0])
            GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
        
        
        t1,y1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2,y2 = doubleListFromFile(GESel2FileName, isFloat=True)
        
        tMin = min(min(t1),min(t2)) 
        tMax = max(max(t1),max(t2))
        tLen = max(len(t1), len(t2))
        t3 = np.linspace(tMin, tMax, tLen)
        
        y3 = np.zeros(len(t3))
        for i in range(len(y3)):
            currT = t3[i]
            tmpInds = np.nonzero(t1<=currT)[0]
            if len(tmpInds) == 0:
                y1LInd = 0
            else:
                y1LInd = tmpInds[-1]
            tmpInds = np.nonzero(t1>=currT)[0]
            if len(tmpInds) == 0:
                y1RInd = len(t1)-1
            else:
                y1RInd = tmpInds[0]
            
            tmpInds = np.nonzero(t2<=currT)[0]
            if len(tmpInds) == 0:
                y2LInd = 0
            else:
                y2LInd = tmpInds[-1]
            tmpInds = np.nonzero(t2>=currT)[0]
            if len(tmpInds) == 0:
                y2RInd = len(t2)-1
            else:
                y2RInd = tmpInds[0]
            
            y1Tmp = (y1[y1RInd]+y1[y1LInd])/2            
            y2Tmp = (y2[y2RInd]+y2[y2LInd])/2
            
            y3[i] = y1Tmp - y2Tmp
        
        
        y3i = np.nonzero(t3<=tOn)[0][-1]
                
        yDelta = abs(y3[y3i])
        print yDelta, t3[y3i]
        if yDelta < FRDelta:
            
            
            FCList = []
            RTList = []
            for theta in thetaList:
                
                try:
                    
                    while abs(y3[y3i]) < theta:
                        y3i += 1
                    
                    if y3[y3i] >= theta:
                        FCList.append(1)
                    else:
                        FCList.append(0)
                    RTList.append(t3[y3i])
                
                except:
                    FCList.append(-1)
                    RTList.append(float("inf"))

            tripleListToFile(thetaList, RTList, FCList, fileName)
            if verbose:
                print "UUID " + UUID + " tested and saved."
        
        else:
            
            if verbose:
                print "UUID " + UUID + " Rejected (deltaFR)."
            return 0
                
                
       

#-------------------------------------------------------------------------------

def thresholdTestSpikesUUID(UUID, thetaList, verbose=1, tOn = 0, namePrefix="InputSel"):
    
    fileName = "thresholdTestSpikes_" + UUID + ".dat"

    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (spikes): " + UUID
        
        GESel1FileName = findFileName([UUID, ".ntf", namePrefix+"1"])[0]
        GESel2FileName = findFileName([UUID, ".ntf", namePrefix+"2"])[0]
        who1,t1 = doubleListFromFile(GESel1FileName)
        who2,t2 = doubleListFromFile(GESel2FileName)


        t1 = [float(val) for val in t1]
        t2 = [float(val) for val in t2]

        spikeCounter = 0
        t1i = 0
        t2i = 0
        currTime = 0
        FCList = []
        RTList = []
        spikeCounter = 0
        for theta in thetaList:

            try:
                while abs(spikeCounter) < theta:
                    if t1[t1i] < t2[t2i]:
                        currTime = t1[t1i]
                        t1i += 1
                        if currTime > tOn:
                            spikeCounter += 1
                    
                    else:
                        currTime = t2[t2i]
                        t2i += 1
                        if currTime > tOn:
                            spikeCounter -= 1
                    
                if spikeCounter >= theta:
                    FCList.append(1)
                else:
                    FCList.append(0)
                RTList.append(currTime)

            except:
                FCList.append(-1)
                RTList.append(float("inf"))
        
        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList


#-------------------------------------------------------------------------------

def thresholdTestSpikesBGTooUUID(UUID, thetaList, verbose=1, tOn = 0, tau=100, beta=0):
    
    betaDict = {}
    betaDict[0]=''
    betaDict[.01]='0.01'
    betaDict[.05]='0.05'
    betaDict[.1]='0.10'
    betaDict[1]='1.00'
    
    fileName = "thresholdTestSpikesBGToo" +betaDict[beta]+ "_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (spikes, BG): " + UUID
        
        GESel1FileName = findFileName([UUID, ".ntf", "InputSel1"])[0]
        GESel2FileName = findFileName([UUID, ".ntf", "InputSel2"])[0]
        BGESel1FileName = findFileName([UUID, ".ntf", "BGESel1"])[0]
        BGESel2FileName = findFileName([UUID, ".ntf", "BGESel2"])[0]
        who1,t1 = doubleListFromFile(GESel1FileName)
        who2,t2 = doubleListFromFile(GESel2FileName)
        whoBG1,tBG1 = doubleListFromFile(BGESel1FileName)
        whoBG2,tBG2 = doubleListFromFile(BGESel2FileName)
        
        t1 = [float(val) for val in t1]
        t2 = [float(val) for val in t2]
        tBG1 = [float(val) for val in tBG1]
        tBG2 = [float(val) for val in tBG2]
        
        t1i = 0
        t2i = 0
        tBG1i = 0
        tBG2i = 0
        currTime = 0
        FCList = []
        RTList = []
            
        spikeCounter = 0        
        for theta in thetaList:
            try:
                while abs(spikeCounter) < theta:
                    oldTime = currTime
                    currTime = min(t1[t1i], t2[t2i], tBG1[tBG1i], tBG2[tBG2i])
                    if currTime == t1[t1i]:
                        t1i += 1
                        if currTime > tOn:
                            spikeCounter += 1
                    if currTime == t2[t2i]:
                        t2i += 1
                        if currTime > tOn:
                            spikeCounter -= 1
                    if currTime == tBG1[tBG1i]:
                        tBG1i += 1
                        if currTime > tOn:
                            spikeCounter += 1
                    if currTime == tBG2[tBG2i]:
                        tBG2i += 1
                        if currTime > tOn:
                            spikeCounter -= 1

                    spikeCounter *= np.exp((currTime-oldTime)/tau*beta)                   
            
                if spikeCounter >= theta:
                    FCList.append(1)
                else:
                    FCList.append(0)
                RTList.append(currTime)
            
            except:
                FCList.append(-1)
                RTList.append(float("inf"))
        
        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList

#-------------------------------------------------------------------------------
def thresholdTestCurrentUUID(UUID, thetaList, verbose=1, tOn = 0):
    
    fileName = "thresholdTestCurrent_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (current): " + UUID
        
        GESel1FileName = findFileName([UUID, ".dat", "GESel1PoolInput"])[0]
        GESel2FileName = findFileName([UUID, ".dat", "GESel2PoolInput"])[0]
        t1,x1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2,x2 = doubleListFromFile(GESel2FileName, isFloat=True)
        
        y1 = -x1.cumsum()
        y2 = -x2.cumsum()

        y3 = y1-y2
        
        if tOn == 0:
            ti = 0
        else:
            ti = np.nonzero(np.array(t1)<tOn)[0][-1] + 1
        
        currTime = 0
        FCList = []
        RTList = []
        
        for theta in thetaList:
            
            try:
                while abs(y3[ti]) < theta:
                    ti += 1
                
                if y3[ti] >= theta:
                    FCList.append(1)
                else:
                    FCList.append(0)
                RTList.append(t1[ti])
            
            except:
                FCList.append(-1)
                RTList.append(float("inf"))
        
        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList

#-------------------------------------------------------------------------------
def thresholdTestSynapseUUID(UUID, thetaList, verbose=1, tOn = 0):
    
    fileName = "thresholdTestSynapse_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (current): " + UUID
        
        GESel1FileName = findFileName([UUID, ".dat", "GESel1SInputSum"])[0]
        GESel2FileName = findFileName([UUID, ".dat", "GESel2SInputSum"])[0]
        t1,x1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2,x2 = doubleListFromFile(GESel2FileName, isFloat=True)

        if tOn == 0:
            ti = 0
        else:
            ti = np.nonzero(np.array(t1)<tOn)[0][-1] + 1
    
        y1 = x1[ti:].cumsum()
        y2 = x2[ti:].cumsum()
        
        y3 = y1-y2
    
                
        
        currTime = 0
        FCList = []
        RTList = []
        ti = 0
        
        for theta in thetaList:
            
            try:
                while abs(y3[ti]) < theta:
                    ti += 1
                
                if y3[ti] >= theta:
                    FCList.append(1)
                else:
                    FCList.append(0)
                RTList.append(t1[ti])
            
            except:
                FCList.append(-1)
                RTList.append(float("inf"))
        
        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList


#-------------------------------------------------------------------------------
def thresholdTestCurrentBGTooUUID(UUID, thetaList, verbose=1, tOn = 0):
    
    fileName = "thresholdTestCurrentBGToo_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (current BG): " + UUID
        
        GESel1FileName = findFileName([UUID, ".dat", "GESel1PoolInput"])[0]
        GESel2FileName = findFileName([UUID, ".dat", "GESel2PoolInput"])[0]
        GESel1BGFileName = findFileName([UUID, ".dat", "GESel1PoolBG"])[0]
        GESel2BGFileName = findFileName([UUID, ".dat", "GESel2PoolBG"])[0]
        t1,x1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2,x2 = doubleListFromFile(GESel2FileName, isFloat=True)
        tBG1,xBG1 = doubleListFromFile(GESel1BGFileName, isFloat=True)
        tBG2,xBG2 = doubleListFromFile(GESel2BGFileName, isFloat=True)
        
        x1 += xBG1
        x2 += xBG2
        
        y1 = -x1.cumsum()
        y2 = -x2.cumsum()

        y3 = y1-y2
        
        if tOn == 0:
            ti = 0
        else:
            ti = np.nonzero(np.array(t1)<tOn)[0][-1] + 1
        
        currTime = 0
        FCList = []
        RTList = []
        
        for theta in thetaList:
            
            try:
                while abs(y3[ti]) < theta:
                    ti += 1
                
                if y3[ti] >= theta:
                    FCList.append(1)
                else:
                    FCList.append(0)
                RTList.append(t1[ti])
            
            except:
                FCList.append(-1)
                RTList.append(float("inf"))
        
        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList

#-------------------------------------------------------------------------------
def thresholdTestSynapseBGTooUUID(UUID, thetaList, verbose=1, tOn = 0):
    
    fileName = "thresholdTestSynapseBGToo_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (current BG): " + UUID
        
        GESel1FileName = findFileName([UUID, ".dat", "GESel1SInputSum"])[0]
        GESel2FileName = findFileName([UUID, ".dat", "GESel2SInputSum"])[0]
        GESel1BGFileName = findFileName([UUID, ".dat", "GESel1SBGSum"])[0]
        GESel2BGFileName = findFileName([UUID, ".dat", "GESel2SBGSum"])[0]
        t1,x1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2,x2 = doubleListFromFile(GESel2FileName, isFloat=True)
        tBG1,xBG1 = doubleListFromFile(GESel1BGFileName, isFloat=True)
        tBG2,xBG2 = doubleListFromFile(GESel2BGFileName, isFloat=True)
        
        x1 += xBG1
        x2 += xBG2

        if tOn == 0:
            ti = 0
        else:
            ti = np.nonzero(np.array(t1)<tOn)[0][-1] + 1
    
    
        y1 = x1[ti:].cumsum()
        y2 = x2[ti:].cumsum()
        
        y3 = y1-y2
        
        currTime = 0
        FCList = []
        RTList = []
            
        ti=0
        for theta in thetaList:
            
            try:
                while abs(y3[ti]) < theta:
                    ti += 1
                
                if y3[ti] >= theta:
                    FCList.append(1)
                else:
                    FCList.append(0)
                RTList.append(t1[ti])
            
            except:
                FCList.append(-1)
                RTList.append(float("inf"))
        
        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList

#-------------------------------------------------------------------------------
def splitCurrentAnalysisBG(UUID, thetaList, verbose=1, tBGOn=2000, tInputOn = 6000):
    
    fileName = "splitCurrentAnalysisBG_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (current BG): " + UUID
        
        GESel1FileName = findFileName([UUID, ".dat", "GESel1PoolInput"])[0]
        GESel2FileName = findFileName([UUID, ".dat", "GESel2PoolInput"])[0]
        GESel1BGFileName = findFileName([UUID, ".dat", "GESel1PoolBG"])[0]
        GESel2BGFileName = findFileName([UUID, ".dat", "GESel2PoolBG"])[0]
        
        t1,x1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2,x2 = doubleListFromFile(GESel2FileName, isFloat=True)
        tBG1,xBG1 = doubleListFromFile(GESel1BGFileName, isFloat=True)
        tBG2,xBG2 = doubleListFromFile(GESel2BGFileName, isFloat=True)
        
        tiInputOn = np.nonzero(np.array(t1)<tInputOn)[0][-1] + 1
        tiBGOn = np.nonzero(np.array(t1)<tBGOn)[0][-1] + 1
        
        x1 = -x1[tiInputOn:]
        x2 = -x2[tiInputOn:]
        xBG1 = -xBG1[tiBGOn:(len(x1)+tiBGOn)]
        xBG2 = -xBG2[tiBGOn:(len(x2)+tiBGOn)]
        
        
        t3 = tBG1[0:len(xBG1)]
        y3 = ((x1+xBG1)-(x2+xBG2)).cumsum()
        
        FCList = []
        RTList = []
        ti = 0
        for theta in thetaList:
            
            try:
                while abs(y3[ti]) < theta:
                    ti += 1
                
                if y3[ti] >= theta:
                    FCList.append(1)
                else:
                    FCList.append(0)
                RTList.append(t3[ti])
            
            except:
                FCList.append(-1)
                RTList.append(float("inf"))
        
        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList

#-------------------------------------------------------------------------------
def splitCurrentAnalysisAll(UUID, thetaList, verbose=1, tBGOn=2000, tInputOn = 6000):
    
    fileName = "splitCurrentAnalysisAll_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (current BG): " + UUID
        
        GESel1FileName = findFileName([UUID, ".dat", "GESel1PoolInput"])[0]
        GESel2FileName = findFileName([UUID, ".dat", "GESel2PoolInput"])[0]
        GESel1BGFileName = findFileName([UUID, ".dat", "GESel1PoolBG"])[0]
        GESel2BGFileName = findFileName([UUID, ".dat", "GESel2PoolBG"])[0]
        GESel1RecAMPAFileName = findFileName([UUID, ".dat", "GESel1PoolRecAMPA"])[0]
        GESel2RecAMPAFileName = findFileName([UUID, ".dat", "GESel2PoolRecAMPA"])[0]
        GESel1RecNMDAFileName = findFileName([UUID, ".dat", "GESel1PoolRecNMDA"])[0]
        GESel2RecNMDAFileName = findFileName([UUID, ".dat", "GESel2PoolRecNMDA"])[0]
        GESel1RecGABAFileName = findFileName([UUID, ".dat", "GESel1PoolRecGABA"])[0]
        GESel2RecGABAFileName = findFileName([UUID, ".dat", "GESel2PoolRecGABA"])[0]

        t1,x1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2,x2 = doubleListFromFile(GESel2FileName, isFloat=True)
        tBG1,xBG1 = doubleListFromFile(GESel1BGFileName, isFloat=True)
        tBG2,xBG2 = doubleListFromFile(GESel2BGFileName, isFloat=True)
        tAMPA1, xAMPA1 = doubleListFromFile(GESel1RecAMPAFileName, isFloat=True)
        tAMPA2, xAMPA2 = doubleListFromFile(GESel2RecAMPAFileName, isFloat=True)
        tNMDA1, xNDMA1 = doubleListFromFile(GESel1RecNMDAFileName, isFloat=True)
        tNMDA2, xNDMA2 = doubleListFromFile(GESel2RecNMDAFileName, isFloat=True)
        tGABA1, xGABA1 = doubleListFromFile(GESel1RecGABAFileName, isFloat=True)
        tGABA2, xGABA2 = doubleListFromFile(GESel2RecGABAFileName, isFloat=True)

        tiInputOn = np.nonzero(np.array(t1)<tInputOn)[0][-1] + 1
        tiBGOn = np.nonzero(np.array(t1)<tBGOn)[0][-1] + 1
        
        x1 = -x1[tiInputOn:]
        x2 = -x2[tiInputOn:]
        xBG1 = -xBG1[tiBGOn:(len(x1)+tiBGOn)]
        xBG2 = -xBG2[tiBGOn:(len(x2)+tiBGOn)]
        xAMPA1 = -xAMPA1[tiBGOn:(len(x1)+tiBGOn)]
        xAMPA2 = -xAMPA2[tiBGOn:(len(x1)+tiBGOn)]
        xNDMA1 = -xNDMA1[tiBGOn:(len(x1)+tiBGOn)]
        xNDMA2 = -xNDMA2[tiBGOn:(len(x1)+tiBGOn)]    
        xGABA1 = -xGABA1[tiBGOn:(len(x1)+tiBGOn)]
        xGABA2 = -xGABA2[tiBGOn:(len(x1)+tiBGOn)]


        t3 = tBG1[0:len(xBG1)]
        y3 = ((x1+xBG1+xAMPA1+xNDMA1+xGABA1)-(x2+xBG2+xAMPA2+xNDMA2+xGABA2)).cumsum()

        FCList = []
        RTList = []
        ti = 0
        for theta in thetaList:
            
            try:
                while abs(y3[ti]) < theta:
                    ti += 1
                
                if y3[ti] >= theta:
                    FCList.append(1)
                else:
                    FCList.append(0)
                RTList.append(t3[ti])
            
            except:
                FCList.append(-1)
                RTList.append(float("inf"))
        
        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList

#-------------------------------------------------------------------------------
def totalSpikeDiffCompare(UUID, thetaList, verbose=1, tOn = 2000):
    
    fileName = "totalSpikeDiffCompare_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        thetaList, RTList, FCList = tripleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (current BG): " + UUID
        
        GESel1FileName = findFileName([UUID, ".ntf", "GESel1Dummy"])[0]
        GESel2FileName = findFileName([UUID, ".ntf", "GESel2Dummy"])[0]
        
        who1,t1 = doubleListFromFile(GESel1FileName, isFloat=True)
        who2,t2 = doubleListFromFile(GESel2FileName, isFloat=True)
        
    

    
        ti1 = np.nonzero(np.array(t1)<tOn)[0][-1] + 1
        ti2 = np.nonzero(np.array(t2)<tOn)[0][-1] + 1
        
        print ti1, ti2
        

        FCList = []
        RTList = []
        spikeCounter = 0
        for theta in thetaList:
            
            try:
                while abs(spikeCounter) < theta:
                    if t1[ti1] < t2[ti2]:
                        currTime = t1[ti1]
                        ti1 += 1
                        spikeCounter += 1
                    
                    else:
                        currTime = t2[ti2]
                        ti2 += 1
                        spikeCounter -= 1
                
                if spikeCounter >= theta:
                    FCList.append(1)
                else:
                    FCList.append(0)
                RTList.append(currTime)
        
            except:
                FCList.append(-1)
                RTList.append(float("inf"))    
        
        tripleListToFile(thetaList, RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList

#-------------------------------------------------------------------------------
def totalSpikeIntCompare(UUID, RTList, verbose=1):
    
    fileName = "totalSpikeIntCompare_" + UUID + ".dat"
    
    if os.path.isfile(fileName):
        RTList, FCList = doubleListFromFile(fileName)
        if verbose:
            print "UUID " + UUID + " loaded."
    else:
        
        if verbose:
            print "Testing UUID (current BG): " + UUID
        
        GESel1FileName = findFileName([UUID, ".ntf", "GESel1Dummy"])[0]
        GESel2FileName = findFileName([UUID, ".ntf", "GESel2Dummy"])[0]
        
        who1,t1 = doubleListFromFile(GESel1FileName, isFloat=True)
        who2,t2 = doubleListFromFile(GESel2FileName, isFloat=True)

        tOn = RTList[0]
        ti1 = np.nonzero(np.array(t1)<tOn)[0][-1] + 1
        ti2 = np.nonzero(np.array(t2)<tOn)[0][-1] + 1
        
        FCList = []
        spikeCounter1 = 0
        spikeCounter2 = 0
        for RT in RTList:
            
            while t1[ti1] < RT and ti1 < (len(t1)-1):
                ti1 += 1
                spikeCounter1 += 1
                
            while t2[ti2] < RT and ti2 < (len(t2)-1):
                ti2 += 1
                spikeCounter2 += 1
                
            if spikeCounter1 == spikeCounter2:
                if random.random() < .5:
                    FCList.append(1)
                else:
                    FCList.append(0)
            elif spikeCounter1 > spikeCounter2:
                FCList.append(1)
            else:
                FCList.append(0)
        
        doubleListToFile(RTList, FCList, fileName)
        if verbose:
            print "UUID " + UUID + " tested and saved."
    return RTList, FCList

##-------------------------------------------------------------------------------
#def thresholdTestCurrentUUID(UUID, thetaList, verbose=1, tOn = 0):
#    
#    fileName = "thresholdTestCurrent_" + UUID + ".dat"
#    
#    if os.path.isfile(fileName):
#        thetaList, RTList, FCList = tripleListFromFile(fileName)
#        if verbose:
#            print "UUID " + UUID + " loaded."
#    else:
#        
#        if verbose:
#            print "Testing UUID (current): " + UUID
#
#        GESel1FileName = findFileName([UUID, ".dat", "GESel1PoolInput"])[0]
#        GESel2FileName = findFileName([UUID, ".dat", "GESel2PoolInput"])[0]
#        t1,x1 = doubleListFromFile(GESel1FileName, isFloat=True)
#        t2,x2 = doubleListFromFile(GESel2FileName, isFloat=True)
#
#        y1 = -x1.cumsum()
#        y2 = -x2.cumsum()
#        
#        if tOn == 0:
#            ti = 0
#        else:
#            ti = np.nonzero(np.array(t1)<tOn)[0][-1] + 1
#
#        currTime = 0
#        FCList = []
#        RTList = []
#
#        for theta in thetaList:
#            
#            try:
#                while y1[ti] < theta and y2[ti] < theta:
#                    ti += 1
#                
#                if y1[ti] >= theta:
#                    FCList.append(1)
#                else:
#                    FCList.append(0)
#                RTList.append(t1[ti])
#            
#            except:
#                FCList.append(-1)
#                RTList.append(float("inf"))
#        
#        tripleListToFile(thetaList, RTList, FCList, fileName)
#        if verbose:
#            print "UUID " + UUID + " tested and saved."
#    return RTList, FCList
#
#
##-------------------------------------------------------------------------------
#def thresholdTestCurrentBGTooUUID(UUID, thetaList, verbose=1, tOn = 0):
#    
#    fileName = "thresholdTestCurrentBGToo_" + UUID + ".dat"
#    
#    if os.path.isfile(fileName):
#        thetaList, RTList, FCList = tripleListFromFile(fileName)
#        if verbose:
#            print "UUID " + UUID + " loaded."
#    else:
#        
#        if verbose:
#            print "Testing UUID (current BG): " + UUID
#        
#        GESel1FileName = findFileName([UUID, ".dat", "GESel1PoolInput"])[0]
#        GESel2FileName = findFileName([UUID, ".dat", "GESel2PoolInput"])[0]
#        GESel1BGFileName = findFileName([UUID, ".dat", "GESel1PoolBG"])[0]
#        GESel2BGFileName = findFileName([UUID, ".dat", "GESel2PoolBG"])[0]
#        t1,x1 = doubleListFromFile(GESel1FileName, isFloat=True)
#        t2,x2 = doubleListFromFile(GESel2FileName, isFloat=True)
#        tBG1,xBG1 = doubleListFromFile(GESel1BGFileName, isFloat=True)
#        tBG2,xBG2 = doubleListFromFile(GESel2BGFileName, isFloat=True)
#        
#        x1 += xBG1
#        x2 += xBG2
#
#        y1 = -x1.cumsum()
#        y2 = -x2.cumsum()
#        
#        if tOn == 0:
#            ti = 0
#        else:
#            ti = np.nonzero(np.array(t1)<tOn)[0][-1] + 1
#        
#        currTime = 0
#        FCList = []
#        RTList = []
#        
#        for theta in thetaList:
#            
#            try:
#                while y1[ti] < theta and y2[ti] < theta:
#                    ti += 1
#                
#                if y1[ti] >= theta:
#                    FCList.append(1)
#                else:
#                    FCList.append(0)
#                RTList.append(t1[ti])
#            
#            except:
#                FCList.append(-1)
#                RTList.append(float("inf"))
#        
#        tripleListToFile(thetaList, RTList, FCList, fileName)
#        if verbose:
#            print "UUID " + UUID + " tested and saved."
#    return RTList, FCList



















#-------------------------------------------------------------------------------

def getRoitmanData(subject='n'):

    CohVals = [0, 3.2, 6.4, 12.8, 25.6, 51.2]
    if subject == 'both':
        CohValsB, RTValsB, FCValsB, CohDictB = getRoitmanData(subject='b')
        CohValsN, RTValsN, FCValsN, CohDictN = getRoitmanData(subject='n')
        

        
        CohVals = []
        for i in range(len(CohValsB[0])):
            CohVals.append((CohValsB[0][i] + CohValsN[0][i])/2)
        
        RTVals = []
        for i in range(len(RTValsB[0])):
            RTVals.append((RTValsB[0][i] + RTValsN[0][i])/2)
                
        FCVals = []
        for i in range(len(FCValsB[0])):
            FCVals.append((FCValsB[0][i] + FCValsN[0][i])/2)
        
        CohDict = {}
        for key in CohDictB[0].keys():
            CohDict[key] = [0,0]
            for i in [0,1]:
                CohDict[key][i] = (CohDictN[0][key][i]+CohDictB[0][key][i])/2



        return [CohVals], [RTVals], [FCVals], [CohDict]
        
    elif subject == 'b':
        FCVals = [0.50462962962963, 0.615560640732265, 0.738532110091743, 0.93348623853211, 0.995412844036697, 1    ]
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
    
def speedAccTradeoff(data, C, plotStyle='-', dots=True, highlightTheta=None, subject = 'n', TND=350,thetaMin=1):

    thetaVals = data.keys()
    CVals = data[thetaVals[0]][0]
    CInd = CVals.index(C)
    RT = []
    FC = []
    
    if dots == True:
        CohVals, RTVals, FCVals, CohDict = getRoitmanData(subject=subject)
        
    if not highlightTheta == None:
        specialThetaRT = data[highlightTheta][1][CInd]
        specialThetaFC = data[highlightTheta][2][CInd]
    
    for theta in thetaVals:
        if theta >= thetaMin:
            RT.append(TND + data[theta][1][CInd])
            FC.append(data[theta][2][CInd])
    
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

    return RT,FC
    
#-------------------------------------------------------------------------------
        
def speedAccTradeoffFull(data, plotStyle='-', dots=True, highlightTheta=None, subject = 'n',TND=350,thetaMin=1):

    CohVals = [0, 3.2, 6.4, 12.8, 25.6, 51.2]
    
    for CohVal in CohVals:
        speedAccTradeoff(data, CohVal, plotStyle=plotStyle, dots=dots, highlightTheta=highlightTheta, subject=subject, TND=TND,thetaMin=thetaMin)
    
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
def plotFR(filename,plotstyle = '-', figure=1,tIn=None):

    if isType(filename,'fr'):
        t,y = doubleListFromFile(filename)
        pl.figure(figure)
        pl.plot(t,y,plotstyle)
    elif isType(filename,'ntf'):
        data = getDataFromFile(filename)
        t,y = firingRate(data, 
                         t=tIn,
                         window=50, 
                         nameList='All', 
                         plotsOn=False, 
                         dt=1,
                         xlabel='t',
                         ylabel='Firing Rate')
        pl.plot(t,y,plotstyle)


    return 0

#-------------------------------------------------------------------------------
def plotPhaseSpace(filename1,filename2,plotstyle = '-', figure=1):
    
    if isType(filename1,'fr'):
        t1,y1 = doubleListFromFile(filename1)
        
    if isType(filename2,'fr'):
        t2,y2 = doubleListFromFile(filename2)
        
        
        
        pl.figure(figure)
        pl.plot(y1[0:len(y1)],y2[0:len(y1)],plotstyle)
    
    return 0

#-------------------------------------------------------------------------------
def plotFRUUID(UUID,plotstyle = '-', figure=1):
    
    fileNameList = findFileName([".fr",UUID])
    
    for f in fileNameList:
        plotFR(f,plotstyle = plotstyle, figure=figure)
    
    return 0

#-------------------------------------------------------------------------------
def plotFRDir(whichDir="./", plotstyle = '-', figure=1):

    fileNameList = findFileName([".fr"], whichDir=whichDir)
    
    for f in fileNameList:
        plotFR(f,plotstyle = plotstyle, figure=figure)
    
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
        
    for i in range(len(namesList)):
        pl.plot(nameTimeDict[namesList[i]],[i]*len(nameTimeDict[namesList[i]]),'b.',markersize=1)
    
    return

#-------------------------------------------------------------------------------
def speedAccTradeoffSpikeIntUnCorr(C, N=240, r0 = 40, bP = .4, bN = .4, 
                                     dots = True, thetaMax=15, RTMax = 900,
                                     plotStyle='-', TND=350, highlightTheta=False, makePlot=True):

    from math import log
    
    rP = (r0 + bP*C)
    rN = (r0 - bN*C)
    
    h0 = log(rN/rP)
    Ez = N*(rP-rN)
    
    if makePlot == True:
        SpeedAccTradeoffExact(h0, Ez, thetaMax, RTMax=RTMax, TND=TND, plotStyle = plotStyle)

        if dots == True:
            plotDots(C, subject='both', RTMin = TND, RTMax = RTMax)
        
        if highlightTheta != False:
            FCSpecial = FCvsTheta([highlightTheta],h0,Ez)
            RTSpecial = RTvsTheta([highlightTheta],h0,Ez,TND=TND)
            pl.plot(RTSpecial,FCSpecial,'*')
        return
    else:
        theta=pl.linspace(0,thetaMax, 1000)
        RT=RTvsTheta(theta,h0,Ez,TND=TND)
        FC=FCvsTheta(theta,h0,Ez)

        if highlightTheta != False:
            FCSpecial = FCvsTheta([highlightTheta],h0,Ez)
            RTSpecial = RTvsTheta([highlightTheta],h0,Ez,TND=TND)
            
            return (RT, FC, theta), (RTSpecial, FCSpecial, highlightTheta)
        else:
        
            return RT, FC, theta

#-------------------------------------------------------------------------------
def speedAccTradeoffSPRTUnCorr(C, N=240, r0 = 40, bP = .4, bN = .4, 
                                     dots = True, thetaMax=15, RTMax = 900,
                                     plotStyle='-', TND=350, highlightTheta=False, makePlot=True, thetaN=10):
    
    from math import log
    
    rP = (r0 + bP*C)
    rN = (r0 - bN*C)
    
    h0 = -1
    Ez = N*(rP-rN)*log(rP/rN)
    
    if makePlot == True:
    
        SpeedAccTradeoffExact(h0, Ez, thetaMax, RTMax=RTMax, TND=TND, plotStyle = plotStyle)
        
        if dots == True:
            plotDots(C, subject='both', RTMin = TND, RTMax = RTMax)
        
        if highlightTheta != False:
            FCSpecial = FCvsTheta([highlightTheta],h0,Ez)
            RTSpecial = RTvsTheta([highlightTheta],h0,Ez,TND=TND)
            pl.plot(RTSpecial,FCSpecial,'*')
        
        return

    else:
        theta = pl.linspace(.01, thetaMax, thetaN)
        return theta, RTvsTheta(theta,h0,Ez,TND=TND), FCvsTheta(theta,h0,Ez)

#-------------------------------------------------------------------------------
def speedAccTradeoffSpikeIntSIP(C, corr, N=240, r0 = 40, bP = .4, bN = .4, 
                                dots = True, thetaMax=15, thetaN=10,
                                RTMax = 900,
                                plotStyle='-', TND=350, makePlot = True, 
                                quickName = False,
                                whichRun = 0, saveResultDir = './savedResults'):
    
    from math import log
    from math import exp
    from scipy.optimize import bisect
    
    rP = (r0 + bP*C)
    rN = (r0 - bN*C)
    
    def phi(t):
        return (1-corr)*(exp(t)-1)*N*rP + corr*rP*(exp(N*t)-1) + (1-corr)*(exp(-t)-1)*N*rN + corr*rN*(exp(-N*t)-1)
    
    
    
    h0 = bisect(phi, -2,-.000001)
    Ez = N*(rP-rN)
    
    if makePlot == True:
        SpeedAccTradeoffExact(h0, Ez, thetaMax, RTMax=RTMax, TND=TND, plotStyle = plotStyle)
        if dots == True:
            plotDots(C, subject='both', RTMin = TND, RTMax = RTMax)
    else:
        if quickName == False:
            theta = pl.linspace(.01, thetaMax, thetaN)
            return theta, RTvsTheta(theta,h0,Ez,TND=TND), FCvsTheta(theta,h0,Ez)
        else:
            import analysisTools as at
            theta, overshootP, overshootN = at.overshootCorrection(
                                                                   {'XVar':'theta','N':N,'dt':.1,'corr':corr,'rP':rP,'rN':rN},
                                                                   saveResultDir = saveResultDir, 
                                                                   quickName = quickName, 
                                                                   whichRun = whichRun)
            
            return theta, overshootP, overshootN, RTvsTheta(theta+overshootP,h0,Ez,TND=TND), FCvsTheta(theta+overshootP,h0,Ez)

#            return theta, overshootP, overshootN, RTvsThetaOvershoot(theta,overshootP,overshootN,h0,Ez,TND=TND), FCvsThetaOvershoot(theta,overshootP,overshootN,h0,Ez)

    return

#-------------------------------------------------------------------------------
def speedAccTradeoffSPRTSIP(C, corr, N=240, r0 = 40, bP = .4, bN = .4, 
                                dots = True, thetaMax=15, thetaN=10,
                                RTMax = 900,
                                plotStyle='-', TND=350, makePlot = True, 
                                quickName = False,
                                whichRun = 0, saveResultDir = './savedResults'):
    
    from math import log
    
    rP = (r0 + bP*C)
    rN = (r0 - bN*C)
    
    h0 = -1
    Ez = (corr+N*(1-corr))*(rN*log(rN*1.0/rP) + rP*log(rP*1.0/rN))
    
#    print corr, h0, Ez
    
    if makePlot == True:
        SpeedAccTradeoffExact(h0, Ez, thetaMax, RTMax=RTMax, TND=TND, plotStyle = plotStyle)
        if dots == True:
            plotDots(C, subject='both', RTMin = TND, RTMax = RTMax)
    else:
        if quickName == False:
            theta = list(np.arange(0,thetaMax,np.log(rP/rN)))[1:]
#            theta = pl.linspace(.01, thetaMax, thetaN)
            return theta, RTvsTheta(theta,h0,Ez,TND=TND), FCvsTheta(theta,h0,Ez)
        else:
            import analysisTools as at
            theta, overshootP, overshootN = at.overshootCorrection(
                    {'XVar':'theta','N':N,'dt':.1,'corr':corr,'rP':rP,'rN':rN},
                     saveResultDir = saveResultDir, 
                     quickName = quickName, 
                     whichRun = whichRun)
            
            
            
            return theta, overshootP, overshootN, RTvsTheta(theta+np.mean(overshootP),h0,Ez,TND=TND), FCvsTheta(theta+np.mean(overshootP),h0,Ez)
    
    return

#-------------------------------------------------------------------------------
def speedAccTradeoffSpikeIntMIP(C, corr, N=240, r0 = 40, bP = .4, bN = .4, 
                                dots = True, thetaMax=15, thetaN=10,
                                RTMax = 900,
                                plotStyle='-', TND=350, makePlot = True, 
                                quickName = False,
                                whichRun = 0, saveResultDir = './savedResults'):
    
    from math import exp
    from scipy.optimize import bisect
    
    rP = (r0 + bP*C)
    rN = (r0 - bN*C)
    
    def phi(t):
        return 1/corr*( rN*((1+corr*(exp(-t)-1))**N-1)  +  rP*((1+corr*(exp(t)-1))**N-1) )
    
    h0 = bisect(phi, -1,-.0001)
    Ez = N*(rP-rN)

    if makePlot == True:
        SpeedAccTradeoffExact(h0, Ez, thetaMax, RTMax=RTMax, TND=TND, plotStyle = plotStyle)
        if dots == True:
            plotDots(C, subject='both', RTMin = TND, RTMax = RTMax)
    else:
        
        if quickName == False:
            theta = pl.linspace(.01, thetaMax, thetaN)
            return theta, RTvsTheta(theta,h0,Ez,TND=TND), FCvsTheta(theta,h0,Ez)
        else:
            import analysisTools as at
            theta, overshootP, overshootN = at.overshootCorrection(
                                                                   {'XVar':'theta','N':N,'dt':.1,'corr':corr,'rP':rP,'rN':rN},
                                                                   saveResultDir = saveResultDir, 
                                                                   quickName = quickName, 
                                                                   whichRun = whichRun)
            
            print np.mean(overshootP)
            
            return theta, overshootP, overshootN, RTvsTheta(theta+overshootP,h0,Ez,TND=TND), FCvsTheta(theta+overshootP,h0,Ez)
    
    return

#-------------------------------------------------------------------------------
def speedAccTradeoffSPRTMIP(C, corr, N=240, r0 = 40, bP = .4, bN = .4, 
                            dots = True, thetaMax=15, thetaN=10,
                            RTMax = 900,
                            plotStyle='-', TND=350, makePlot = True, 
                            quickName = False,
                            whichRun = 0, saveResultDir = './savedResults'):
    
    from math import log
    
    rP = (r0 + bP*C)
    rN = (r0 - bN*C)
    
    h0 = -1
    Ez = (1-(1-corr)**N)/corr*(rN*log(rN/rP) + rP*log(rP/rN))
    
#    print corr, h0, Ez
    
    if makePlot == True:
        SpeedAccTradeoffExact(h0, Ez, thetaMax, RTMax=RTMax, TND=TND, plotStyle = plotStyle)
        if dots == True:
            plotDots(C, subject='both', RTMin = TND, RTMax = RTMax)
    else:
        
        if quickName == False:
            theta = pl.linspace(.01, thetaMax, thetaN)
#            theta = list(np.arange(0,thetaMax,np.log(rP/rN)))[1:]
            return theta, RTvsTheta(theta,h0,Ez,TND=TND), FCvsTheta(theta,h0,Ez)
        else:
            import analysisTools as at
            theta, overshootP, overshootN = at.overshootCorrection(
                                                                   {'XVar':'theta','N':N,'dt':.1,'corr':corr,'rP':rP,'rN':rN},
                                                                   saveResultDir = saveResultDir, 
                                                                   quickName = quickName, 
                                                                   whichRun = whichRun)
            
            return theta, overshootP, overshootN, RTvsTheta(theta+overshootP,h0,Ez,TND=TND), FCvsTheta(theta+overshootP,h0,Ez)
    
    return

#-------------------------------------------------------------------------------
def FCvsTheta(theta,h0,Ez):

    from math import exp
    
    FC = [0]*len(theta)
    for i in range(len(FC)):
        FC[i] = 1/(1+exp(h0*theta[i]))
    
    return FC

#-------------------------------------------------------------------------------
def RTvsTheta(theta,h0,Ez, TND=350):

    from math import tanh
    
    RT = [0]*len(theta)
    for i in range(len(RT)):
        RT[i] = TND+ theta[i]/Ez*tanh(-h0*theta[i]/2)*1000
    
    return RT

##-------------------------------------------------------------------------------
def FCvsThetaOvershoot(theta,overshootP,overshootN,h0,Ez):
    
    from math import exp
    
    print theta,overshootP,overshootN
    
    FC = [0]*len(theta)
    for i in range(len(FC)):
        FC[i] = 1 - (exp(h0*(theta[i]+overshootP[i])) - 1)/(exp(h0*(theta[i]+overshootP[i]))-exp(h0*(-theta[i]+overshootN[i])))
    
#    print FC
    
    return FC

#-------------------------------------------------------------------------------
def RTvsThetaOvershoot(theta,overshootP,overshootN,h0,Ez, TND=350):
    
    from math import tanh
    
#    print theta,overshootP,overshootN
    
    FC = FCvsThetaOvershoot(theta,overshootP,overshootN,h0,Ez)
    
    RT = [0]*len(theta)
    for i in range(len(RT)):
        RT[i] = TND + ((theta[i]+overshootP[i])*(FC[i])+(-theta[i]+overshootN[i])*(1-FC[i]))/Ez*1000
    
#    print RT
    
    return RT

#-------------------------------------------------------------------------------
def SpeedAccTradeoffExact(h0, Ez, thetaMax, RTMax=900, TND=350, plotStyle = 'b-'):

    FCMin = .5
    FCMax = 1
    RTMax = RTMax
    RTMin = TND
    
    theta = pl.linspace(0,thetaMax, 1000)
    
    pl.plot(RTvsTheta(theta,h0,Ez,TND=TND), FCvsTheta(theta,h0,Ez))
    
    print theta[-1], RTvsTheta(theta,h0,Ez,TND=TND)[-1], FCvsTheta(theta,h0,Ez)[-1]
    
    pl.xlim([RTMin, RTMax])
    pl.ylim([FCMin, FCMax])

    return

#-------------------------------------------------------------------------------
def plotDots(C, subject='both', RTMin = 350, RTMax = 900):

    CohVals, RTValsRoitman, FCValsRoitman, CohDict = getRoitmanData(subject=subject)
    
    if subject == 'both':
        pl.plot(CohDict[0][C][0], CohDict[0][C][1],'bo')
        pl.plot(CohDict[1][C][0], CohDict[1][C][1],'ro')
    elif subject == 'b':
        pl.plot(CohDict[0][C][0], CohDict[0][C][1],'bo')
    elif subject == 'n':
        pl.plot(CohDict[0][C][0], CohDict[0][C][1],'ro')

    pl.xlim([RTMin, RTMax])
    pl.ylim([.5, 1])

    return



#-------------------------------------------------------------------------------
def speedAccTradeoffSPRTUnCorrBGToo(C, BGRate=2400, N=240, r0 = 40, bP = .4, bN = .4, 
                               dots = True, thetaMax=15, RTMax = 900,
                               plotStyle='-', TND=350):
    
    from math import log
    
    rP = BGRate + (r0 + bP*C)
    rN = BGRate + (r0 - bN*C)
    
    h0 = -1
    Ez = N*(rP-rN)*log(rP/rN)
    
    SpeedAccTradeoffExact(h0, Ez, thetaMax, RTMax=RTMax, TND=TND, plotStyle = plotStyle)
    
    if dots == True:
        plotDots(C, subject='both', RTMin = TND, RTMax = RTMax)
    
    return


#-------------------------------------------------------------------------------
def loadAllCurrents(UUID):
    
#    GESel1AllFileName = findFileName([UUID, ".dat", "GESel1DummyAll"])[0]
#    GESel2AllFileName = findFileName([UUID, ".dat", "GESel2DummyAll"])[0]
    GESel1RecAMPAFileName = findFileName([UUID, ".dat", "GESel1PoolRecAMPA"])[0]
    GESel2RecAMPAFileName = findFileName([UUID, ".dat", "GESel2PoolRecAMPA"])[0]
    GESel1RecNMDAFileName = findFileName([UUID, ".dat", "GESel1PoolRecNMDA"])[0]
    GESel2RecNMDAFileName = findFileName([UUID, ".dat", "GESel2PoolRecNMDA"])[0]
    GESel1RecGABAFileName = findFileName([UUID, ".dat", "GESel1PoolRecGABA"])[0]
    GESel2RecGABAFileName = findFileName([UUID, ".dat", "GESel2PoolRecGABA"])[0]
    GESel1InputFileName = findFileName([UUID, ".dat", "GESel1PoolInput"])[0]
    GESel2InputFileName = findFileName([UUID, ".dat", "GESel2PoolInput"])[0]
    GESel1BGFileName = findFileName([UUID, ".dat", "GESel1PoolBG"])[0]
    GESel2BGFileName = findFileName([UUID, ".dat", "GESel2PoolBG"])[0]
    
#    t, All1 = doubleListFromFile(GESel1AllFileName, isFloat=True)
#    t, All2 = doubleListFromFile(GESel2AllFileName, isFloat=True)
    t, xAMPA1 = doubleListFromFile(GESel1RecAMPAFileName, isFloat=True)
    t, xAMPA2 = doubleListFromFile(GESel2RecAMPAFileName, isFloat=True)
    t, xNMDA1 = doubleListFromFile(GESel1RecNMDAFileName, isFloat=True)
    t, xNMDA2 = doubleListFromFile(GESel2RecNMDAFileName, isFloat=True)
    t, xGABA1 = doubleListFromFile(GESel1RecGABAFileName, isFloat=True)
    t, xGABA2 = doubleListFromFile(GESel2RecGABAFileName, isFloat=True)
    t, xInput1 = doubleListFromFile(GESel1InputFileName, isFloat=True)
    t, xInput2 = doubleListFromFile(GESel2InputFileName, isFloat=True)
    t, xBG1 = doubleListFromFile(GESel1BGFileName, isFloat=True)
    t, xBG2 = doubleListFromFile(GESel2BGFileName, isFloat=True)

    return (t,
            (-xAMPA1,-xAMPA2),
            (-xNMDA1,-xNMDA2),
            (-xGABA1,-xGABA2),
            (-xInput1, -xInput2),
            (-xBG1, -xBG2))

#-------------------------------------------------------------------------------
def plotFRSel1Sel2(UUID, figureInd=1):

    GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
    GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]

    t1, fr1 = doubleListFromFile(GESel1FileName, isFloat=True)
    t2, fr2 = doubleListFromFile(GESel2FileName, isFloat=True)
    
    pl.figure(figureInd)
    if fr1[-1] > fr2[-1]:
        pl.plot(t1,fr1,'g')
        pl.plot(t2,fr2,'r')
    else:
        pl.plot(t1,fr1,'r')
        pl.plot(t2,fr2,'g')

#-------------------------------------------------------------------------------
def plotFRSel1Sel2Dir(dir="./",  figure="Same",N=-1):
    
    UUIDList = getUUIDList(dir=dir)
    
    print len(UUIDList)
    
    for i in range(len(UUIDList[0:N])):
        if figure == "Same":
            plotFRSel1Sel2(UUIDList[i], figureInd=1)
        else:
            plotFRSel1Sel2(UUIDList[i], figureInd=figure)


#-------------------------------------------------------------------------------
def doubleOnTrials(dir="./", FRDelta=10):
    
    UUIDList = getUUIDList(dir=dir)
    
    nTrials = len(UUIDList)
    
    print nTrials
    
    doubleOn = 0
    counter = 0
    for UUID in UUIDList:
        
        counter  += 1
        print counter, " (" + str(nTrials) + ")"
        
        GESel1FileName = findFileName([UUID, ".fr", "GESel1"])[0]
        GESel2FileName = findFileName([UUID, ".fr", "GESel2"])[0]
        
        t1, fr1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2, fr2 = doubleListFromFile(GESel2FileName, isFloat=True)

        if abs(fr1[-1]-fr2[-1]) < FRDelta:
            doubleOn += 1

    return (doubleOn*1.0/(nTrials*1.0), doubleOn)

#-------------------------------------------------------------------------------
def getSel1Sel2(UUID, tOn = 2000, what='FR', changeCoord=False,whichDir="./"):
    
    if what == 'FR':
        s1="GESel1"
        s2="GESel2"
        type=".fr"
    elif what == 'NMDA':
        s1="GESel1PoolRecNMDA"
        s2="GESel2PoolRecNMDA"
        type=".dat"

    GESel1FileName = findFileName([UUID, type, s1], whichDir=whichDir)[0]
    GESel2FileName = findFileName([UUID, type, s2], whichDir=whichDir)[0]
    
    t1, y1 = doubleListFromFile(GESel1FileName, isFloat=True, whichDir=whichDir)
    t2, y2 = doubleListFromFile(GESel2FileName, isFloat=True, whichDir=whichDir)


    if what == "FR":
        tMin = min(min(t1),min(t2)) 
        tMax = max(max(t1),max(t2))
        tLen = max(len(t1), len(t2))
        t3 = np.linspace(tMin, tMax, tLen)
        
        y3 = np.zeros(len(t3))
        y3neg = np.zeros(len(t3))
        y1New = np.zeros(len(t3))
        y2New = np.zeros(len(t3))
        for i in range(len(y3)):
            currT = t3[i]
            tmpInds = np.nonzero(t1<=currT)[0]
            if len(tmpInds) == 0:
                y1LInd = 0
            else:
                y1LInd = tmpInds[-1]
            tmpInds = np.nonzero(t1>=currT)[0]
            if len(tmpInds) == 0:
                y1RInd = len(t1)-1
            else:
                y1RInd = tmpInds[0]
            
            tmpInds = np.nonzero(t2<=currT)[0]
            if len(tmpInds) == 0:
                y2LInd = 0
            else:
                y2LInd = tmpInds[-1]
            tmpInds = np.nonzero(t2>=currT)[0]
            if len(tmpInds) == 0:
                y2RInd = len(t2)-1
            else:
                y2RInd = tmpInds[0]
            
            y1Tmp = (y1[y1RInd]+y1[y1LInd])/2            
            y2Tmp = (y2[y2RInd]+y2[y2LInd])/2
            
            y3[i] = y1Tmp - y2Tmp
            y3neg[i] = y1Tmp + y2Tmp

            y1New[i] = y1Tmp
            y2New[i] = y2Tmp
    elif what == "NMDA":
        y3=-(y1-y2)
        y3neg=-(y1+y2)
        t3=t1
    




    tOnInd = np.nonzero(t3<tOn)[0][-1] + 1

    if changeCoord == True:
        return t3[tOnInd:], y3[tOnInd:], y3neg[tOnInd:]
    elif changeCoord == False:
        return t3[tOnInd:], y1New[tOnInd:], y2New[tOnInd:]


#-------------------------------------------------------------------------------
def plotSel1Sel2DiffDir(dir="./",  figure=1, what='FR', ICDelta="Inf",N=-1):


    UUIDList = getUUIDList(dir=dir)
    
    print len(UUIDList)
    
    for i in range(len(UUIDList[0:N])):
        print i
        t,y,yneg=getSel1Sel2(UUIDList[i], what=what)
        if abs(y[0]) < float(ICDelta):
            print abs(y[0])
            pl.figure(figure)
#            if y[-1]>0:
#                pl.plot(y,yneg,'g')
#            else:  
#                pl.plot(y,yneg,'r')
            pl.plot(t,y-yneg)


    return

def smooth(x,window_len=3,window='hanning'):
    """smooth the data using a window with requested size.
        
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        
        input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
        flat window will produce a moving average smoothing.
        
        output:
        the smoothed signal
        
        example:
        
        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)
        
        see also: 
        
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter
        
        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        """
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    
    
    if window_len<3:
        return x
    
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    
    
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

def ICStudy(dir="./", what='FR', tOn=2000, plot=0, save=0, loadFile=None):

    if loadFile != None:
        result=pt.unpickle(loadFile)
    else:
    
        UUIDList = getUUIDList(dir=dir)

        if what == 'FR':
            s1="GESel1"
            s2="GESel2"
            type=".fr"
        elif what == 'NMDA':
            s1="GESel1PoolRecNMDA"
            s2="GESel2PoolRecNMDA"
            type=".dat"
        
        
        result = []
        for UUID in UUIDList:
            

            GESel1FileName = findFileName([UUID, type, s1])[0]
            GESel2FileName = findFileName([UUID, type, s2])[0]

            t1, y1 = doubleListFromFile(GESel1FileName, isFloat=True)
            t2, y2 = doubleListFromFile(GESel2FileName, isFloat=True)

            ti1 = np.nonzero(np.array(t1)<=tOn)[0][-1] + 1
            ti2 = np.nonzero(np.array(t2)<=tOn)[0][-1] + 1

            if y1[-1]-y2[-1] > 0:
                correct = 1
            else:
                correct = 0

            result.append([[y1[ti1]-y2[ti2],y1[ti1]+y2[ti2]],correct])
#            result.append([[y1[-1]-y2[-1],y1[-1]+y2[-1]],correct])

        if save == 1:
            pt.pickle(result,"ICStudy_"+what+".dat")

    if plot==1:
        for tuple, correct in result:
            if correct == 1:
                pl.plot([tuple[0]],[tuple[1]],'.g')
            else:
                pl.plot([tuple[0]],[tuple[1]],'.r')

    return result


#-------------------------------------------------------------------------------
def plotFRSel1Sel2SSDiffDir(dir="./",N=-1,figure=1):
    
    UUIDList = getUUIDList(dir=dir)
    
    print len(UUIDList)
    
    endValPlus=[]
    endValMinus=[]
    for i in range(len(UUIDList[0:N])):
    
        GESel1FileName = findFileName([UUIDList[i], ".fr", "GESel1"])[0]
        GESel2FileName = findFileName([UUIDList[i], ".fr", "GESel2"])[0]
    
        t1, fr1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2, fr2 = doubleListFromFile(GESel2FileName, isFloat=True)
    
        if fr1[-1] > fr2[-1]:
            endValPlus.append(fr1[-1]-fr2[-1])
        else:
            endValMinus.append(fr1[-1]-fr2[-1])

    endValPlus = np.array(endValPlus)
    endValMinus = np.array(endValMinus)

    print endValPlus.mean(), endValPlus.std(), endValPlus.std()/np.sqrt(len(endValPlus))
    print endValMinus.mean(), endValMinus.std(), endValMinus.std()/np.sqrt(len(endValMinus))
    
#    pl.figure(figure)
#    pl.hist(endValPlus,20)
#    pl.hist(endValMinus,20)

#-------------------------------------------------------------------------------
def meanFRmeanSBGmeanIBG(dir='./', skipT=500, labelList=["0"],type='Inh'):
    
    
    if type=='Inh':
        typeStr = "I"
    elif type=='Ex':
        typeStr = 'NSel'
    
#    UUIDList = getUUIDList(dir=dir)
#
##    print UUIDList
#    for UUID in UUIDList:
        
    for label in labelList:
#            GESel1FRFileName = findFileName([UUID, ".fr", "GESel"+str(i)+"_"])[0]
#            GESel1SBGSumFileName = findFileName([UUID, ".dat", "GESel"+str(i)+"ISGaba_"])[0]
#            GESel1IBGSumFileName = findFileName([UUID, ".dat", "GESel"+str(i)+"IBGSum_"])[0]

#            GESel1FRFileName = findFileName([UUID, ".fr", "GESel"+str(i)+"_"])[0]
        GESel1SBGSumFileName = findFileName([".dat", "G" + typeStr + "ShadowS" + label + "_"])[0]
        GESel1IBGSumFileName = findFileName([".dat", "G" + typeStr + "ShadowI" + label + "_"])[0]
        
#            tFR, FR = doubleListFromFile(GESel1FRFileName, isFloat=True)    
        tS, S = doubleListFromFile(GESel1SBGSumFileName, isFloat=True)
        tI, I = doubleListFromFile(GESel1IBGSumFileName, isFloat=True)


#            FRLi = np.nonzero(np.array(tFR)<=skipT/2)[0][-1] + 1
#            FRRi = np.nonzero(np.array(tFR)<=tFR[-1]-skipT/2)[0][-1] + 1
        SLi = np.nonzero(np.array(tS)<=skipT/2)[0][-1] + 1
        SRi = np.nonzero(np.array(tS)<=tS[-1]-skipT/2)[0][-1] + 1
        ILi = np.nonzero(np.array(tI)<=skipT/2)[0][-1] + 1
        IRi = np.nonzero(np.array(tI)<=tI[-1]-skipT/2)[0][-1] + 1
        
#            FRNew = FR[FRLi:FRRi]
#            tFRNew = tFR[FRLi:FRRi]

        SNew = S[SLi:SRi]
        tSNew = tS[SLi:SRi]
        
        INew = I[ILi:IRi]
        tINew = tI[ILi:IRi]
        
#            print FRNew.mean(), FRNew.var(), SNew.mean(), SNew.var(), INew.mean(), INew.var()
        print SNew.mean(), INew.mean(), INew.std()

    return 

#-------------------------------------------------------------------------------
def solve1DModel(dir='./',N=1, alpha=34.63, beta=6.4, betaGammaRatio=1.05):

    UUIDList = getUUIDList(dir=dir)
    
    
    endValPlus=[]
    endValMinus=[]
    counter = 0
    for UUID in UUIDList[0:N]:
        counter += 1
        print counter
        GESel1SBGSumFileName = findFileName([UUID, ".dat", "GESel1SBGSum"])[0]
        GESel2SBGSumFileName = findFileName([UUID, ".dat", "GESel2SBGSum"])[0]
        GESel1SInputSumFileName = findFileName([UUID, ".dat", "GESel1SInputSum"])[0]
        GESel2SInputSumFileName = findFileName([UUID, ".dat", "GESel2SInputSum"])[0]

        t, SB1 = doubleListFromFile(GESel1SBGSumFileName, isFloat=True)
        t, SB2 = doubleListFromFile(GESel2SBGSumFileName, isFloat=True)
        t, SI1 = doubleListFromFile(GESel1SInputSumFileName, isFloat=True)    
        t, SI2 = doubleListFromFile(GESel2SInputSumFileName, isFloat=True)
        
        Input=SI1+SB1-SI2-SB2

        ansSim=np.zeros(len(t))
        for i in range(1,len(t)):
            ansSim[i] = ansSim[i-1]+(t[1]-t[0])/1000*(-2*beta*ansSim[i-1]*(ansSim[i-1]/alpha-1)*(ansSim[i-1]/alpha+1)+(beta/1.05)*Input[i-1])
        
        if ansSim[-1] > 0:
            endValPlus.append(ansSim[-1])
        else:
            endValMinus.append(ansSim[-1])

        pl.plot(t,ansSim)
        pl.xlim([0,3000])
        pl.ylim([-50,50])


    endValPlus = np.array(endValPlus)
    endValMinus = np.array(endValMinus)
    print endValPlus.mean(), endValPlus.std(), endValPlus.std()/np.sqrt(len(endValPlus)) 
    print endValMinus.mean(), endValMinus.std(), endValMinus.std()/np.sqrt(len(endValMinus))
    print len(endValPlus)*1.0/(len(endValPlus)+len(endValMinus))

            




#-------------------------------------------------------------------------------
def plotdEdtVsE(dir="./",N=-1,figure=1):
    
    UUIDList = getUUIDList(dir=dir)
    
    print len(UUIDList)
    
    endValPlus=[]
    endValMinus=[]
    
    deVals = np.arange(-60,60,.1)[0:-1]+.05
    EVals = np.arange(-60,60,.1)[0:-1]+.05
    myBins = np.arange(-60,60,.1)
    
    for i in range(len(UUIDList[0:N])):
        
        GESel1FileName = findFileName([UUIDList[i], ".fr", "GESel1"])[0]
        GESel2FileName = findFileName([UUIDList[i], ".fr", "GESel2"])[0]
        
        t1, fr1 = doubleListFromFile(GESel1FileName, isFloat=True)
        t2, fr2 = doubleListFromFile(GESel2FileName, isFloat=True)
        
        frDiff = fr1-fr2
        dfrDiff = np.diff(frDiff)
    


        deVals += pl.histogram(dfrDiff,bins=myBins)[0]
        

        pl.figure(figure)
        pl.plot(EVals,deVals/N,'k')

#    pl.plot([-60,60],[0,0],'k--')









