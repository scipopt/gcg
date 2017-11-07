#!/usr/bin/python

import os
import sys


def main(argv):
	if len(argv) != 2:
		print "Usage: ./plotclassifier.py myOutfile.out classifiername"
		return

	#declaring needed data structures
	classnames = {}
	classnmembers = {}
	instancenames = []
	classifiernames = []
	blockcandidates = {}
	blockcandidatesnvotes = {}

	with open(argv[0]) as f:
		while True:
			line = f.readline()
			if not line: break
			if not line.startswith("Start writing complete"):
				continue
			else:
				#start handling information
				line = f.readline() #line now contains instance information
				print line
				line = line.split()
				instancename = line[1]
				instancename = instancename.split('/')
				instancename = instancename[len(instancename)-1]
				instancenames.append(instancename)
				classnames[instancename] = {}
				classnmembers[instancename] = {}
				blockcandidates[instancename] = []
				blockcandidatesnvotes[instancename] = []
				line = f.readline() #line now contains n blockcandidates
				nblockcandidates = int(line)
				for blockcand in range(nblockcandidates):
					#handle blockcandidates
					line = f.readline() #line now contains information for one blockcandidate and its number of votes
					line = line.split()
					blockcandidates[instancename].append(int(line[0]))
					blockcandidatesnvotes[instancename].append(int(line[2]))
				line = f.readline() # line now contains n cons classifer
				nconsclassifier = int(line)
				for consclassifier in range(nconsclassifier):
					line = f.readline()
					classifiername = line
					line = f.readline()
					nclasses = int(line)
					classnames[instancename][classifiername] = []
					classnmembers[instancename][classifiername] = []
					if classifiername not in classifiernames:
						classifiernames.append(classifiername)
					for classid in range(nclasses):
						line = f.readline()
						line = line.split(':')
						classname = line[0]
						line = f.readline()
						nmembers = int(line)
						classnames[instancename][classifiername].append(classname)
						classnmembers[instancename][classifiername].append(nmembers)
				continue
	print instancenames
	print blockcandidates
	print blockcandidatesnvotes
	print classnames
	print classnmembers


	

if __name__ == "__main__": 
	main(sys.argv[1:])

