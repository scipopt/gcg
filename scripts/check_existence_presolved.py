#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import commands

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
    else:
        instancename2 = instancename       #        print instancename2
    if instancename2.endswith(".mps"):
        instancename = instancename2[:-4]
    else:
        instancename = instancename2
    return instancename



def main(argv):

    infile = open(argv[1], "r")
    instances = infile.readlines()
    counter = 0
    failedones = ""
    for instance in instances:
        instance = instance.strip()
#        print instance
        pathtofile = getpath(instance)
 #       print pathtofile
        instancename = getinstancename(instance)
  #      print instancename
        presolvedpath = pathtofile + "presolved_instances"
        jobscriptname = instancename + "job.sh"
        jobscript = open(jobscriptname, "w")

        pathtonewfile = pathtofile + "presolved_instances/p_" + instancename + ".mps.gz"
        if not os.path.isfile(pathtonewfile):
            counter = counter + 1
            failedones = failedones + instance + "\n"
#        print str(counter) + " fails, "
    print "the following instances are missing:\n" +  failedones
    pass

if __name__ == "__main__":
    if len(sys.argv) == 2:
        main(sys.argv)
    else:
        print "usage: createBlockFiles.py <testset> <scip> <settings>"
