#!/usr/bin/python

import os
import sys
from subprocess import call

if len(sys.argv) == 4 or len(sys.argv) == 5:

	abspathtototestset = sys.argv[1]
#	instancename = sys.argv[2]
	pathtogcg = sys.argv[2]
	ndecomps = -1
	testsettoread = sys.argv[3]
	testsetname = testsettoread
	if len(sys.argv ) >= 5:
		ndecomps = sys.argv[4]

	gcgexec = pathtogcg+"bin/gcg"
#	settings = pathtogcg+"settings/onlyConsclassConnected.set"
	settings = "default.set"
	testsetfile = pathtogcg+"check/testset/allDec_"+testsetname+".test"
	testsettoreadfile = pathtogcg + "check/testset/" + testsettoread+".test"
	tsfile = open(testsettoreadfile,"r")

	print pathtogcg


	# (1) create decs folder
	call(['mkdir', pathtogcg+"check/decs/"+testsetname+"All" ])

	for line in tsfile:

		help = line.split("/")
		help = help[len(help)-1]
		help = help.split(".")
		help2 =help[0]

		for i in range(1, len(help)-2):
			if( help[i] == "gz" or help[i] == "lp" or help[i] == "mps"):
				continue

			help2 = help2 + help[i]

		print help
		instancename = help2
#		abspathtoinstance = pathtogcg+"check/"+line
		abspathtoinstance = line

		abspathtoinstance = abspathtoinstance.rstrip()

		# (3) write all decomps
		print ""+gcgexec + ' -c '+ "read " + abspathtoinstance, '-c ' + "set load " + settings + ' -c ' + "detect " + '-c' + " write alldec " + pathtogcg+"check/decs/"+testsetname + " dec " + '-c' + " quit"
		call([gcgexec, '-c', "read " + abspathtoinstance, '-c', "set load " + settings, '-c', "set visual nmaxdecompstowrite " + ndecomps, '-c', "detect", '-c', "write alldec " + pathtogcg+"check/decs/"+testsetname + " dec", '-c', "quit" ] )

		# (4) create testset-file
		call(['touch', testsetfile])

		# (5) fill testset-file
		os.system('for i in `ls ' + pathtogcg+'check/decs/'+testsetname+'/*'+instancename+'*.dec | sort -V`;do echo "'+abspathtoinstance+';$i" >> '+testsetfile+'; done ')
		print 'for i in `ls ' + pathtogcg+'check/decs/'+testsetname+'/*'+instancename+'*.dec | sort -V`;do echo "'+abspathtoinstance+';$i" >> '+testsetfile+'; done '


    # (6) call make test
	os.system('cd ' + pathtogcg)

	print 'make test TEST=allDec_'+ testsetname+ ' LPS=spx OPT=opt TIME=1800 '
	#os.system('make test TEST=allDec_'+ instancename+ ' LPS=spx OPT=opt TIME=1800 ')
else:
	print "Usage: ./writeandtestalldecomps.py abspathtotestset testsetname abspathtogcg [ndecompsperinstance]"
