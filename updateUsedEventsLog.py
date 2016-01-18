#! /usr/bin/env python
"""
    This script generates a log file which contains the ip-glasma
    events that have been used before.
"""
from os import path, remove
from glob import glob
from sys import exit
from ParameterDict import *

def prepareEventsList(numberOfEvents):
    ipglasmaDirectory = ipglasmaParameters['dataPath']
    # get data directory
    collisionSystem = ipglasmaParameters['projectile']+ipglasmaParameters['target']
    ecm = ipglasmaParameters['ecm']
    centrality = centralityParameters['centrality'][:-1]
    ipglasmaDataDirectory = path.join(ipglasmaDirectory, '%s_%d'%(collisionSystem, ecm),
        '%s%s'%(collisionSystem, centrality))
    if not path.exists(ipglasmaDataDirectory):
        print "No initial files for specified collision system or centrality:"
        print "%s at %d AGeV for %s"%(collisionSystem, ecm, centrality)
        exit()
    # get a list of all available files
    fileList =  glob(path.join(ipglasmaDataDirectory, ipglasmaParameters['filePattern']))
    fileList.sort() # sort files by name
    # get a list of used initial files
    try:
        with open('usedEventsLog.dat', 'r') as usedlogFile:
            usedFileList = [line.rstrip('\n') for line in usedlogFile]
    except:
        print "No usedEventsLog.dat found, use empty used events log."
        usedFileList = []
    # generate a list of up-to-use files
    unusedFileList = []
    for aFile in fileList:
        if aFile not in usedFileList:
            unusedFileList.append(aFile)
    # check if we have enough number of initial conditions
    if len(unusedFileList) < numberOfEvents:
        print "prepareEventsList: Not enough IP-Glasma events!"
        print "unused files: %d, require files: %d"%(len(unusedFileList), numberOfEvents)
        exit()
    # write out the up-to-use files
    with open('initEventsLog.dat', 'w') as f:
        for aFile in unusedFileList[0:numberOfEvents]:
            f.write(aFile + '\n')
    # keep a record of used files
    with open('usedEventsLog.dat', 'a+') as f:
        for aFile in unusedFileList[0:numberOfEvents]:
            f.write(aFile + '\n')


def cleanUsedEventsLog():
    if path.exists('usedEventsLog.dat'):
        remove('usedEventsLog.dat')
    # create an empty file
    f = open('usedEventsLog.dat', 'w')
    f.close()

