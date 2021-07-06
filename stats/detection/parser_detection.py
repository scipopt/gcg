#!/usr/bin/env python3

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.rcsetup as rcsetup

class Dataset:

	def checksection(self, line, keyword, warn = True):
		word = line.split()[0]
		if word == keyword:
			return  True
		else:
			if warn: print("ERROR: line is " + line + " but should be " + keyword)
			return False


	def __init__(self, filename, fromApp = False):
		self.fromApp = fromApp
		self.classnames = {}
		self.classnmembers = {}
		self.instancenames = []
		self.classifiernames = []
		self.blockcandidates = {}
		self.blockcandidatesnvotes = {}
		self.decompscores = {}
		self.decompmaxforwhitescores = {}
		self.decompids = {}
		self.classicalscores = {}
		self.decompssetpartmaster = {}
		self.decompnblocks = {}
		self.maxndecomps = 0
		self.detectiontimes = {}
		self.filename = filename
		nfound = 0
		nfoundnodec = 0
		with open(filename) as f:
			while True:
				line = f.readline()
				if not line: break
				#print(line)
				if line.startswith("Detection did not take place so far"):
					nfoundnodec += 1
				if not line.startswith("Start writing complete detection information"):
					continue
				else:
					#print("found ", nfound)
					nfound += 1
					#start handling information
					line = f.readline() #line now contains instance information
	#				print line
					line = line.split()
					instancename = line[1]
					# workaround for "filename: unknown" occuring multiple times
					if instancename == "unknown": instancename = "unknown_" + str(nfound)
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
					self.classicalscores[instancename] = []
					self.decompmaxforwhitescores[instancename] = []
					self.decompids[instancename] = []
					line = f.readline()
					if not self.checksection(line, "NBLOCKCANDIDATES"): return
					line = f.readline() #line now contains n blockcandidates on third position
					nblockcandidates = int(line.split()[2])
					for blockcand in range(nblockcandidates):
						#handle blockcandidates
						line = f.readline() #line now contains information for one blockcandidate and its number of votes
						line = line.split()
						self.blockcandidates[instancename].append(int(line[0]))
						if line[2] != "user":
							self.blockcandidatesnvotes[instancename].append(int(line[2]))
						else:
							self.blockcandidatesnvotes[instancename].append("user")
					line = f.readline()
					if not self.checksection(line, "DETECTIONTIME"): return
					line = f.readline()
					detectiontime = float(line)
					self.detectiontimes[instancename] = detectiontime
					line = f.readline() # line now contains keyword
					if not (self.checksection(line, "CONSPARTITION", warn=False) or self.checksection(line, "CONSCLASSIFIER",warn=False)): return
					line = f.readline() # line now contains n cons
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
					line = f.readline() # line now contains keyword
					if not (self.checksection(line, "VARPARTITION",warn=False) or self.checksection(line, "VARCLASSIFIER",warn=False)): return
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
					if not self.checksection(line, "DECOMPINFO"): return
					line = f.readline()
					ndecomps = int(line)-1
					for decomp in range(ndecomps):
						line = f.readline()
						if not self.checksection(line, "NEWDECOMP"): return
						line = f.readline()
						nblocks = int(line)
						self.decompnblocks[instancename].append(nblocks)
						line = f.readline()
						decompid = int(line)
						self.decompids[instancename].append(decompid)
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
						self.classicalscores[instancename].append(classicalscore)
						line = f.readline()
						decompmaxforwhitescore = float(line)
						self.decompmaxforwhitescores[instancename].append(decompmaxforwhitescore)
						line = f.readline()
						setpartmaster = int(line)
						self.decompssetpartmaster[instancename].append(setpartmaster)
						line = f.readline()
						ndetectors = int(line)
						for detector in range(ndetectors):
							line = f.readline()
							detectorname = line
						if line.startswith("@04"):
							continue
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

		if self.instancenames == []:
			print("Warning: Data could not be parsed.\n         Have you conducted the test with MODE=detectionstatistics?")
			if not self.fromApp:
				print("Terminating.")
				exit()
			else:
				print("Please choose a different file.")
				return
		print("File:      ", filename.split("/")[-1], "\nInstances: ", nfound+nfoundnodec, "(of which without detection:", nfoundnodec, ")")


	def getmaxdetectiontime(self):
		maxdetectiontime = 0.
		for instance in self.detectiontimes:
			if self.detectiontimes[instance] > maxdetectiontime:
				maxdetectiontime = self.detectiontimes[instance]
		return maxdetectiontime

	def getmaxnblocks(self):
		maxnblocks = 0.
		for instance in self.decompnblocks:
			if len(self.decompnblocks[instance]) == 0:
				continue
			if self.decompnblocks[instance][0] > maxnblocks:
				maxnblocks = self.decompnblocks[instance][0]
		return maxnblocks

	def getmaxnnontrivialdecomps(self):
		maxntrivialdecomps = 0
		for instance in self.decompscores :
			if len(self.decompscores[instance]) == 0:
				continue
			counter = self.getnnontrivialdecompsforinstance(instance)
			if counter > maxntrivialdecomps:
				maxntrivialdecomps = counter
		return maxntrivialdecomps

	def getnnontrivialdecompsforinstance(self, instance):
		counter = 0
		for score in self.decompscores[instance]:
			if score > 0.:
				counter = counter + 1
		return counter



	def getNNonTrivialDecomp( self):
		counter = 0
		for decompscore in self.decompscores:
		#	print decompscores[decompscore]
			if len(self.decompscores[decompscore]) == 0:
				continue
			if self.decompscores[decompscore][0] > 0.:
				counter = counter + 1
		return counter

	def getNNonTrivialDecompSetpartmaster( self):
		counter = 0
		for decompscore in self.decompscores:
		#	print decompscores[decompscore]
			if len(self.decompscores[decompscore]) == 0:
				continue
			decompid = 0
			while self.decompssetpartmaster[decompscore][decompid] != 1:
				decompid = decompid+1
				if decompid == len(self.decompscores[decompscore]):
					break
			if decompid == len(self.decompscores[decompscore]):
				continue
			if self.decompscores[decompscore][decompid] > 0.:
				counter = counter + 1
		return counter



	def fractionofinstanceswithscoreatleastsetpartmaster( self, minscore):
		counter = 0
		for decompscore in self.decompscores:
		#	print decompscores[decompscore]
			if len(self.decompscores[decompscore]) == 0:
				continue
			decompid = 0
			while self.decompssetpartmaster[decompscore][decompid] != 1:
				decompid = decompid+1
				if decompid == len(self.decompscores[decompscore]):
					break
			if decompid == len(self.decompscores[decompscore]):
				continue
			if self.decompscores[decompscore][decompid] >= minscore:
				counter = counter + 1
		return float(counter)/float(len(self.decompscores))


	def fractionofinstanceswithscoreatleast( self, decompscores, minscore):
		counter = 0
		for decompscore in decompscores:
		#	print decompscores[decompscore]
			if len(decompscores[decompscore]) == 0:
				continue
			if decompscores[decompscore][0] >= minscore:
				counter = counter + 1
		return float(counter)/float(len(decompscores))

	def fractionofinstanceswithatleasttaunontrivialdecomps( self, decompscores, tau):
		counter = 0
		for instance in decompscores:
		#	print decompscores[decompscore]
			if self.getnnontrivialdecompsforinstance(instance) >= tau:
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

	def fractionofmemberswithvalatmostwithscoreatleast( self, members, hashmapvalue, maxvalue):
		counter = 0
		for instance in members:
			if instance not in hashmapvalue:
				continue
			if hashmapvalue[instance] <= maxvalue:
				counter = counter + 1
		return float(counter)/float(len(members))


	def fractionofinstanceswithatleasttauclasses( self, classnames, minclasses, classifier):
		counter = 0
		for instance in classnames:
		#	print decompscores[decompscore]
			if len(classnames[instance][classifier]) >= minclasses:
				counter = counter + 1
		return float(counter)/float(len(classnames))


	def getsettingsname(self):
		settingname = self.filename.split(".")
		settingname = settingname[len(settingname)-3]
		return settingname

	def getclassifiernames(self):
		return self.classifiernames

	def getfilename(self):
		return self.filename

	def getninstances(self):
		return len(self.instancenames)



def main(argv):
	print("The parser is not intended to be called.")

if __name__ == "__main__":
	main(sys.argv[1:])
