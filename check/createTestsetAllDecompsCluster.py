#!/usr/bin/python

import os
import sys
from subprocess import call

if len(sys.argv) == 4:

	abspathtoinstance = sys.argv[1]
	instancename = "cl_"+sys.argv[2]
	pathtogcg = sys.argv[3]

	gcgexec = pathtogcg+"bin/gcg"
	settings = pathtogcg+"settings/onlyConsclassConnected.set"
	testsetfile = pathtogcg+"check/testset/allDec_"+instancename+".test"

	print pathtogcg

	# (1) create decs folder
	call(['mkdir', pathtogcg+"check/decs/"+instancename ])

	# (2) write all decomps
	call([gcgexec, '-c', "read " + abspathtoinstance, '-c', "set load " + settings, '-c', "detect", '-c', "write alldec " + pathtogcg+"check/decs/"+instancename + " dec", '-c', "quit" ] )

	# (3) create testset-file
	call(['touch', testsetfile])

	# (4) fill testset-file
	os.system('for i in `ls ' + pathtogcg+'check/decs/'+instancename+'/*.dec | sort -V`;do echo "'+abspathtoinstance+';$i" >> '+testsetfile+'; done ')
    # (5) call make test
	os.system('cd ' + pathtogcg)
	os.system('make testcluster TEST=allDec_'+ instancename+ ' LPS=spx OPT=opt TIME=1800 ')
else:
	print "Usage: ./createTestsetAllDecomps.py abspathtoinstance instancename abspathtogcg"