#!/usr/bin/python

import os
import sys

f = open("instancesMittelmann.txt", "r")

f2 = open("instancesMittelmannClean.txt", "w")

for line in f:
	f2.write(line.split()[0] + "\n")
