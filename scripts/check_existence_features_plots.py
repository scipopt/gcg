#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import commands

def preambel(instancename):
    help = "#!/usr/bin/env zsh\n### Job name\n#BSUB -J  "
    help = help + instancename + "\n\n"
    return help


def midpart(instancename, counter):
    help = "### File / path where STDOUT & STDERR will be written\n###    %J is the job ID, %I is the array ID\n#BSUB -o "
    help = help + instancename + ".%" + str(counter+1) +".1\n"
    return  help +"### Request the time you need for execution in minutes\n### The format for the parameter is: [hour:]minute,\n### that means for 80 minutes you could also use this: 1:20\n#BSUB -W 0:20\n\n### Request memory you need for your job in TOTAL in MB\n#BSUB -M 2048\n\n### Change to the work directory\ncd /home/mb908627/gcg"

def finalpart(instancepath, instancename, gcgpath, settingspath, presolvedinfo, featurespath, matrixpath, decomppath, originstancepath):
    help = "\n\n./scripts/featuresPlots.py " + instancepath + " "  + instancename + " " + gcgpath
    help = help + " " + settingspath + " " + str(presolvedinfo) + " " + featurespath + " " + matrixpath
    help = help + " " +  decomppath + " " + originstancepath
    return  help


def getpath(instancepath):
    splitted = instancepath.split("/")
    newpath = ""
    for tokenid in range(len(splitted)-2):
        newpath += ( splitted[tokenid] + "/" )
    return newpath

def getinstancename(instancepath):
    instance = instancepath.strip()
    instance = instance.split("/")
    instancename = instance[len(instance)-1]
    if instancename.endswith(".gz"):
        instancename2 = instancename[:-3]
#        print instancename2
    else:
        instancename2 = instancename
    if instancename2.endswith(".mps"):
        instancename = instancename2[:-4]
    else:
        instancename = instancename2[:-3]
    return instancename



def main(argv):

    infile = open(argv[1], "r")
    instances = infile.readlines()
    counterfeatures = 0
    countermatrixplot = 0
    counterdecompplot = 0
    counterdecompdec = 0
    presolved = (argv[4] == 'presolved') or (argv[4] == 'True')
    failedfeatures = ""
    failedmatrixplot = ""
    faileddecompplot = ""
    faileddecompdec = ""
    overallfeatures = open("overallfeatures.csv", "w")
    failedinstances = open("check/testset/failedplots.test","w")

    for instance in instances:
        instance = instance.strip()
#        print instance
        pathtofile = getpath(instance)
#        print pathtofile
        instancename = getinstancename(instance)
        if presolved:
            instancename = instancename[2:]
        print instancename
        presolvedpath = pathtofile + "presolved_instances"
        decomppath = pathtofile + "decomps"
        matrixpath = pathtofile + "matrix"
        featurepath = pathtofile + "features_trivpresolved_component_decomps"
        originstancepath = pathtofile + "instances/" + instancename + ".mps.gz"


        jobscriptname = instancename + "job.sh"

#        pathtonewfile = pathtofile + "presolved_instances/p_" + instancename + ".mps"
#        if os.path.isfile(pathtonewfile+".gz"):
 #           continue

        if not os.path.isdir(presolvedpath):
            os.mkdir(presolvedpath)

        if not os.path.isdir(matrixpath):
            os.mkdir(matrixpath)

        if not os.path.isdir(decomppath):
            os.mkdir(decomppath)

        if not os.path.isdir(featurepath):
            os.mkdir(featurepath)

        decomppath = decomppath + "/" + instancename
        featurepath = featurepath + "/" + instancename
        featurepath = featurepath+".csv"
        matrixpath = matrixpath+"/" + instancename + ".gp"
        decompplotpath = decomppath + ".gp"
        decompdecpath = decomppath + ".dec"
        failed = False

        if not os.path.isfile(featurepath):
            counterfeatures = counterfeatures + 1
            failedfeatures = failedfeatures + instance + "\n"
            failed = True
        else:
            featurefileinst = open(featurepath, "r")
            line = featurefileinst.readlines()[1]
            overallfeatures.write(line+"\n")
        if not os.path.isfile(matrixpath):
            countermatrixplot = countermatrixplot + 1
            failedmatrixplot = failedmatrixplot + instance + "\n"
            failed = True

        if not os.path.isfile(decompplotpath):
            counterdecompplot = counterdecompplot + 1
            faileddecompplot = faileddecompplot + instance + "\n"
            failed = True

        if not os.path.isfile(decompdecpath):
            counterdecompdec = counterdecompdec + 1
            faileddecompdec = faileddecompdec + instance + "\n"
            failed = True
        if failed:
            failedinstances.write(instance+"\n")

    overallfeatures.close()
    failedinstances.close()

    print "nfails features: " + str(counterfeatures)
    print "nfails matrix: " + str(countermatrixplot)
    print "nfails decompplot: " + str(counterdecompplot)
    print "nfails decompdec: " + str(counterdecompdec)



    pass

if __name__ == "__main__":
    if len(sys.argv) == 5:
        main(sys.argv)
    else:
        print "usage: featuresPlotsOverall.py <testset> <gcg> <settings> <presolved>"
