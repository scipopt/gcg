#!/usr/bin/python

import os

print "Updating the all.solu file!"

try:
	os.remove('testset/all.solu')
	print "Removing all.solu and creating a new one."
except:
	print "Creating new all.solu file since none exists."

allSolu = open('testset/all.solu', 'w')
info = dict() # for each instance, remember information as a tuple (solufile, status, objvalue)

for file in os.listdir("testset"):
    if file.endswith(".solu"):
		currentFile = open("testset/"+file, 'r')

		for line in currentFile:
                        # check for duplicate instances
                        tokens = line.split()
                        if len(tokens) > 2:
                                value = tokens[2]
                        else:
                                value = 0.0

                        if tokens[1] in info:
                                if tokens[0] != info[tokens[1]][1] or value != info[tokens[1]][2]:
                                        print "WARNING: double information for instance " + tokens[1]
                                        print "    in " + info[tokens[1]][0] + ": " + info[tokens[1]][1] + " " + str(info[tokens[1]][2])
                                        print "    in " + file + ": " + tokens[0] + " " + str(value)
                                continue

                        info[tokens[1]] = (file, tokens[0], value)
			allSolu.write(line)

                currentFile.close()

allSolu.close()
