#!/usr/bin/python

import os

print "Updating the all.solu file!"

try:
	os.remove('testset/all.solu')
	print "Removing all.solu and creating a new one."
except:
	print "Creating new all.solu file since none exists."

allSolu = open('testset/all.solu', 'w')

for file in os.listdir("testset"):
    if file.endswith(".solu"):
		currentFile = open("testset/"+file, 'r')
		for line in currentFile:
			allSolu.write(line)

allSolu.close()
