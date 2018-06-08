#!/usr/bin/python

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.rcsetup as rcsetup

class Plotter:
	def __init__(self, filename):
		self.classnames = {}
		self.classnmembers = {}
		self.instancenames = []
		self.classifiernames = []
		self.blockcandidates = {}
		self.blockcandidatesnvotes = {}
		self.decompscores = {}
		self.decompssetpartmaster = {}
		self.decompnblocks = {}
		self.maxndecomps = 0

		self.filename = filename

		with open(filename) as f:
			while True:
				line = f.readline()
				if not line: break
				if not line.startswith("Start writing complete"):
					continue
				else:
					#start handling information
					line = f.readline() #line now contains instance information
	#				print line
					line = line.split()
					instancename = line[1]
					instancename = instancename.split('/')
					instancename = instancename[len(instancename)-1]
					self.instancenames.append(instancename)
					self.classnames[instancename] = {}
					self.classnmembers[instancename] = {}
					self.blockcandidates[instancename] = []
					self.blockcandidatesnvotes[instancename] = []
					self.decompnblocks[instancename] = []
					self.decompssetpartmaster[instancename] = []
					self.decompscores[instancename] = []

					line = f.readline() #line now contains n blockcandidates
					nblockcandidates = int(line)
					for blockcand in range(nblockcandidates):
						#handle blockcandidates
						line = f.readline() #line now contains information for one blockcandidate and its number of votes
						line = line.split()
						self.blockcandidates[instancename].append(int(line[0]))
						self.blockcandidatesnvotes[instancename].append(int(line[2]))
					line = f.readline() # line now contains n cons classifer
					nconsclassifier = int(line)
					for consclassifier in range(nconsclassifier):
						line = f.readline()
						classifiername = line
						classifiername = classifiername.strip(' \t\n')
						line = f.readline()
						nclasses = int(line)
						self.classnames[instancename][classifiername] = []
						self.classnmembers[instancename][classifiername] = []
						if classifiername not in self.classifiernames:
							self.classifiernames.append(classifiername)
						for classid in range(nclasses):
							line = f.readline()
							line = line.split(':')
							classname = line[0]
							line = f.readline()
							nmembers = int(line)
							self.classnames[instancename][classifiername].append(classname)
							self.classnmembers[instancename][classifiername].append(nmembers)
					line = f.readline() # line now contains n var classifer
					nvarclassifier = int(line)
					for varclassifier in range(nvarclassifier):
						line = f.readline()
						classifiername = line
						classifiername = classifiername.strip(' \t\n')
						line = f.readline()
						nclasses = int(line)
						self.classnames[instancename][classifiername] = []
						self.classnmembers[instancename][classifiername] = []
						if classifiername not in self.classifiernames:
							self.classifiernames.append(classifiername)
						for classid in range(nclasses):
							line = f.readline()
							line = line.split(':')
							classname = line[0]
							line = f.readline()
							nmembers = int(line)
							self.classnames[instancename][classifiername].append(classname)
							self.classnmembers[instancename][classifiername].append(nmembers)
					line = f.readline()
					ndecomps = int(line)
					for decomp in range(ndecomps):
						line = f.readline()
						nblocks = int(line)
						self.decompnblocks[instancename].append(nblocks)
						for block in range(nblocks):
							line = f.readline()
							nconss = int(line)
							line = f.readline()
							nvars = int(line)
						line = f.readline()
						nmasterconss = int(line)
						line = f.readline()
						nlinkingvars = int(line)
						line = f.readline()
						nmastervars = int(line)
						line = f.readline()
						ntotalstairlinking = int(line)
						line = f.readline()
						maxwhitescore = float(line)
						self.decompscores[instancename].append(maxwhitescore)
						line = f.readline()
						classicalscore = float(line)
						line = f.readline()
						setpartmaster = int(line)
						self.decompssetpartmaster[instancename].append(setpartmaster)
						line = f.readline()
						ndetectors = int(line)
						for detector in range(ndetectors):
							line = f.readline()
							detectorname = line
						line = f.readline()
						nconsclassifier = int(line)
						for consclassifier in range(nconsclassifier):
							line = f.readline()
							classifiernamedecomp = line
							line = f.readline()
							nmasterclasses = int(line)
							for masterclass in range(nmasterclasses):
								line = f.readline()
								line = line.split(':')
								masterclassname = line[0]
						line = f.readline()
						nvarclassifier = int(line)
						for varclassifier in range(nvarclassifier):
							line = f.readline()
							varclassifiernamedecomp = line
							line = f.readline()
							nmastervarclasses = int(line)
							for mastervarclass in range(nmastervarclasses):
								line = f.readline()
								line = line.split(':')
								mastervarclassname = line[0]
							line = f.readline()
							nlinkingvarclasses = int(line)
							for linkingvarclass in range(nlinkingvarclasses):
								line = f.readline()
								line = line.split(':')
								linkingvarclassname = line[0]
					continue



	def fractionofinstanceswithscoreatleast( self, decompscores, minscore):
		counter = 0
		for decompscore in decompscores:
		#	print decompscores[decompscore]
			if len(decompscores[decompscore]) == 0:
				continue
			if decompscores[decompscore][0] >= minscore:
				counter = counter + 1
		return float(counter)/float(len(decompscores))

	def fractionofinstanceswithatleasttaudecomps( self, decompscores, minscore):
		counter = 0
		for decompscore in decompscores:
		#	print decompscores[decompscore]
			if len(decompscores[decompscore]) >= minscore:
				counter = counter + 1
		return float(counter)/float(len(decompscores))



	def fractionofinstanceswithnblocksleast( self, decompnblocks, minblocks):
		counter = 0
		for decomp in decompnblocks:
		#	print decompscores[decompscore]
			if len(decompnblocks[decomp]) == 0:
				continue
			if decompnblocks[decomp][0] >= minblocks:
				counter = counter + 1
		return float(counter)/float(len(decompnblocks))

	def fractionofinstanceswithatleasttauclasses( self, classnames, minclasses, classifier):
		counter = 0
		for instance in classnames:
		#	print decompscores[decompscore]
			if len(classnames[instance][classifier]) >= minclasses:
				counter = counter + 1
		return float(counter)/float(len(classnames))


	def plotdetectionquality(self):
		tauvals = np.arange(0., 1., 0.01)
		instancefractions = []
		for tau in tauvals:
			instancefractions.append(self.fractionofinstanceswithscoreatleast(self.decompscores, tau) )
		plt.ylabel('fraction of instances')
		plt.xlabel('Whitest found decomp has at least this max white score')
	#	print tauvals
		#print instancefractions

		plt.plot(tauvals, instancefractions)

		plt.show()

	def plotnblocksofbest(self):
		maxnblocks = 0
		for decomp in self.decompnblocks:
			if self.decompnblocks[decomp][0] > maxnblocks:
				maxnblocks = self.decompnblocks[decomp][0]
		tauvals = np.arange(0., maxnblocks)
		instancefractions = []
		for tau in tauvals:
			instancefractions.append(self.fractionofinstanceswithnblocksleast(self.decompnblocks, tau) )
		plt.ylabel('fraction of instances')
		plt.xlabel('whitest found decomposition has at least this number of blocks ')
	#	print tauvals
		#print instancefractions

		plt.semilogx(tauvals, instancefractions)

		plt.show()

	def plotndecomps(self):
		maxndecomps = 0
		for decomp in self.decompscores:
			if len(self.decompscores[decomp]) > maxndecomps:
				maxndecomps = len(self.decompscores[decomp])
		tauvals = np.arange(0., maxndecomps)
		instancefractions = []
		for tau in tauvals:
			instancefractions.append(self.fractionofinstanceswithatleasttaudecomps(self.decompscores, tau) )
		plt.ylabel('fraction of instances')
		plt.xlabel('at least this number of decompositions is found')

	#	print tauvals
		#print instancefractions

		plt.plot(tauvals, instancefractions)

		plt.show()

	def plotnclassesforclassifier(self, classifier):
		maxnclasses = 0
		for instance in self.instancenames:
			if len(self.classnames[instance][classifier]) > maxnclasses:
				maxnclasses = len(self.classnames[instance][classifier])
		tauvals = np.arange(1., maxnclasses)
		instancefractions = []
		for tau in tauvals:
			instancefractions.append(self.fractionofinstanceswithatleasttauclasses(self.classnames, tau, classifier) )
		plt.ylabel('fraction of instances')
		plt.xlabel('at least this number of classes is found for classifier '+classifier)

	#	print tauvals
		#print instancefractions

		plt.plot(tauvals, instancefractions)

		plt.show()




	def getclassifiernames(self):
		return self.classifiernames

	def getfilename(self):
		return self.filename

	def getninstances(self):
		return len(self.instancenames)



def main(argv):
	if len(argv) != 2:
#		print "Usage: ./plotclassifier.py myOutfile.out classifiername"
		return

	plotty = Plotter(argv[0])
	#end file reading, start

	print(plotty.classifiernames)
	#plot decomp quality
	plotty.plotnclassesforclassifier("constypes")
	plotty.plotdetectionquality()
	plotty.plotnblocksofbest()
	plotty.plotndecomps()

	#print instancenames
	#print blockcandidates
	#print blockcandidatesnvotes
	#print classnames
	#print classnmembers
	#print decompnblocks
	#print decompscores



if __name__ == "__main__":
	main(sys.argv[1:])
