#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import commands


if len(sys.argv) == 2:
    infile = open(sys.argv[1], "r")
    instances = infile.readlines()

    for instance in instances:
        instance = instance.strip()
        print instance
        if instance.replace(".gz", "").endswith(".lp"):
            mpsinstance = instance.replace(".gz", "").replace(".lp",".mps")
            os.system('bin/gcg -c "read check/%s write problem check/%s quit"' %(instance, mpsinstance))
        else:
            mpsinstance = instance
        os.system('bin/gcg -c "read check/%s detect write problem check/%s.ref quit"' %(mpsinstance, mpsinstance))
        os.system('scripts/ref2block.sh check/%s.ref > check/%s.block' %(mpsinstance, mpsinstance))
        os.system('rm -f check/%s.ref' %(mpsinstance))
else:
    print "usage: createBlockFiles.py <testset>"
