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
        print instancename2
    if instancename2.endswith(".mps"):
        instancename = instancename2[:-4]
    return instancename



def main(argv):

    infile = open(argv[1], "r")
    instances = infile.readlines()
    for instance in instances:
        instance = instance.strip()
        print instance
        pathtofile = getpath(instance)
        print pathtofile
        instancename = getinstancename(instance)
        print instancename
        presolvedpath = pathtofile + "presolved_instances"
        
        if not os.path.isdir(presolvedpath):
            os.mkdir(presolvedpath) 

        pathtonewfile = pathtofile + "presolved_instances/p_" + instancename + ".mps"

        scipcommand = argv[2] + " -c \"read " + instance + " set load " + argv[3] + " presolve write trans " + pathtonewfile + " quit \" "

        print scipcommand

        os.system(scipcommand)

#        os.system('rm -f check/%s.ref' %(mpsinstance))
        os.system('gzip %s' %(pathtonewfile))

    pass

if __name__ == "__main__":
    if len(sys.argv) == 4:
        main(sys.argv)
    else:
        print "usage: createBlockFiles.py <testset> <scip> <settings>"
