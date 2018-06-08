#!/usr/bin/python

import os
import sys
from subprocess import call

if len(sys.argv) == 4 or len(sys.argv) == 5:

	abspathtoinstance = sys.argv[1]
	instancename = sys.argv[2]
	pathtogcg = sys.argv[3]
	testsetname = instancename
	if len(sys.argv ) == 5:
		testsetname = sys.argv[4]

	gcgexec = pathtogcg+"bin/gcg"
#	settings = pathtogcg+"settings/onlyConsclassConnected.set"
        settings = "default.set"
	testsetfile = pathtogcg+"check/testset/allDec_"+testsetname+".test"

	print pathtogcg

	# (1) create decs folder
	call(['mkdir', pathtogcg+"check/decs/"+testsetname ])

	# (2) write all decomps
	print ""+gcgexec + ' -c '+ "read " + abspathtoinstance, '-c ' + "set load " + settings + ' -c ' + "detect " + '-c' + " write alldec " + pathtogcg+"check/decs/"+testsetname + " dec " + '-c' + " quit"
	call([gcgexec, '-c', "read " + abspathtoinstance, '-c', "set load " + settings, '-c', "detect", '-c', "write alldec " + pathtogcg+"check/decs/"+testsetname + " dec", '-c', "quit" ] )

	# (3) create testset-file
	call(['touch', testsetfile])

	# (4) fill testset-file
	os.system('for i in `ls ' + pathtogcg+'check/decs/'+testsetname+'/*'+instancename+'*.dec | sort -V`;do echo "'+abspathtoinstance+';$i" >> '+testsetfile+'; done ')
	print 'for i in `ls ' + pathtogcg+'check/decs/'+testsetname+'/*'+instancename+'*.dec | sort -V`;do echo "'+abspathtoinstance+';$i" >> '+testsetfile+'; done '
    # (5) call make test
	os.system('cd ' + pathtogcg)
#	print 'make test TEST=allDec_'+ testsetname+ ' LPS=spx OPT=opt TIME=1800 '
	#os.system('make test TEST=allDec_'+ instancename+ ' LPS=spx OPT=opt TIME=1800 ')
else:
	print "Usage: ./createTestsetAllDecomps.py abspathtoinstance instancename abspathtogcg [testsetname]"
