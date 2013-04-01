#! /usr/bin/env python
"""
    Search and submit (torque) all pbs files inside direct subdirectories under
    the given path.
"""

from sys import argv, exit
from os import path, listdir

# get optional arguments
if len(argv)>=2:
    targetWorkingDirectory = path.abspath(argv[1])
else:
    targetWorkingDirectory = "PlayGround"

# check existence of target working directory
if path.exists(targetWorkingDirectory):
    print("Start submitting jobs from %s" % targetWorkingDirectory)
else:
    print("Usage: submitJobs [from_directory=PlayGround]")
    exit()

from subprocess import call

for aFolder in listdir(targetWorkingDirectory):
    subFolder = path.join(targetWorkingDirectory, aFolder)
    for aFile in listdir(subFolder):
        if path.splitext(aFile)[1].lower() == '.pbs':
            commandString = "qsub %s" % aFile
            print("Submitting %s in %s..." % (aFile, subFolder))
            call(commandString, shell=True, cwd=subFolder)

print("Job submittion done. See RunRecord.txt file in each job folder for progress.")

