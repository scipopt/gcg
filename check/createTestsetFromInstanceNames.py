#!/usr/bin/python

import os
import sys

if len(sys.argv) >= 2:
	allInstancesFile = open(sys.argv[1], 'r')

	testFile = open("testset/"+sys.argv[2]+".test", 'w')

	allInstances = []

	for line in allInstancesFile:
		allInstances.append(line.split()[0])

	#print allInstances

	counter = 0
	
	for instance in allInstances:
		found = 0
		pathToInstance = ""
		for file in os.listdir("testset"):
			if file.endswith(".test"):
				currentFile = open("testset/"+file, 'r')
				for line in currentFile:
					if line.split("/")[len(line.split("/"))-1].split(".")[0] == instance and line.split("/")[0] == "instances":
						pathToInstance = line
						found = 1
						break
				if found == 1:
					break
		if found == 1:
			testFile.write(pathToInstance)
			counter = counter + 1
		else:
			print "Instance " + instance + " was not found!"

	print "Testfile "+ sys.argv[2]+".test" +" created in testfile folder! "+ str(counter) + " of " + str(len(allInstances)) + " instances contained."
			

else:
	print "Usage: ./createTestsetFromInstanceNames.py myFileWithInstanceNames.txt myTestsetName"
	



