#!/usr/bin/python

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.rcsetup as rcsetup
from detectionplotter import *

class Plotter:
    def __init__(self):
        print("Init of plotter /n")


    def plotdetectiontimes(self, datasets):
        maxdetectiontime = 0.
        labels = []
        for dataset in datasets:
            currtime = dataset.getmaxdetectiontime()
            if currtime > maxdetectiontime:
                maxdetectiontime = currtime
        tauvals = np.arange(0, maxdetectiontime*1.1, 1.1*maxdetectiontime/1000.)
        instfractsfordataset = []
        for dataset in datasets:
            instfract = []
            for tau in tauvals:
                instfract.append(dataset.fractionofmemberswithvalatmostwithscoreatleast(dataset.instancenames, dataset.detectiontimes, tau) )
            instfractsfordataset.append(instfract)
        plt.ylabel('fraction of instances')
        plt.xlabel('Detection time is at most')

        plt.axis([0., maxdetectiontime*1.1, 0., 1.])
    #   print tauvals
        #print instancefractions
        for datasetid in range(len(datasets)):
            plt.plot(tauvals, instancefractions[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='upper center', 
           bbox_to_anchor=[0.5, 1.1], 
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)
        
        plt.xscale('symlog')
        plt.show()




    def plotdetectionquality(self):
        tauvals = np.arange(0., 1., 0.01)
        instancefractions = []
        for tau in tauvals:
            instancefractions.append(self.fractionofinstanceswithscoreatleast(self.decompscores, tau) )
        plt.ylabel('fraction of instances')
        plt.xlabel('Whitest found decomp has at least this max white score')
    #   print tauvals
        #print instancefractions

        plt.plot(tauvals, instancefractions)

        plt.show()

    def plotdetectionqualitysetpartmaster(self):
        tauvals = np.arange(0., 1., 0.01)
        instancefractions = []
        for tau in tauvals:
            instancefractions.append(self.fractionofinstanceswithscoreatleastsetpartmaster(tau) )
        plt.ylabel('fraction of instances')
        plt.xlabel('Whitest found decomp with setpartitioning master has at least this max white score')
    #   print tauvals
        #print instancefractions

        plt.plot(tauvals, instancefractions)

        plt.show()



    def plotnblocksofbest(self):
        maxnblocks = 0
        for decomp in self.decompnblocks:
            if self.decompnblocks[decomp][0] > maxnblocks:
                maxnblocks = self.decompnblocks[decomp][0]
        tauvals = np.arange(0., maxnblocks)
        instancefractions = []
        for tau in tauvals:
            instancefractions.append(self.fractionofinstanceswithnblocksleast(self.decompnblocks, tau) )
        plt.ylabel('fraction of instances')
        plt.xlabel('whitest found decomposition has at least this number of blocks ')
        
    #   print tauvals
        #print instancefractions

        plt.semilogx(tauvals, instancefractions)

        plt.show()

    def plotndecomps(self):
        maxndecomps = 0
        for decomp in self.decompscores:
            if len(self.decompscores[decomp]) > maxndecomps:
                maxndecomps = len(self.decompscores[decomp])
        tauvals = np.arange(0., maxndecomps)
        instancefractions = []
        for tau in tauvals:
            instancefractions.append(self.fractionofinstanceswithatleasttaudecomps(self.decompscores, tau) )
        plt.ylabel('fraction of instances')
        plt.xlabel('at least this number of decompositions is found')

    #   print tauvals
        #print instancefractions

        plt.plot(tauvals, instancefractions)

        plt.show()

    def plotnclassesforclassifier(self, classifier):
        maxnclasses = 0
        for instance in self.instancenames:
            if len(self.classnames[instance][classifier]) > maxnclasses:
                maxnclasses = len(self.classnames[instance][classifier])
        tauvals = np.arange(1., maxnclasses)
        instancefractions = []
        for tau in tauvals:
            instancefractions.append(self.fractionofinstanceswithatleasttauclasses(self.classnames, tau, classifier) )
        plt.ylabel('fraction of instances')
        plt.xlabel('at least this number of classes is found for classifier '+classifier)

    #   print tauvals
        #print instancefractions

        plt.plot(tauvals, instancefractions)

        plt.show()




    

def main(argv):
    if len(argv) != 2:
#       print "Usage: ./plotclassifier.py myOutfile.out classifiername"
        return
    
    plotty = Plotter(argv[0])
    #end file reading, start 

    print(plotty.classifiernames)
    #plot decomp quality 
    plotty.plotnclassesforclassifier("constypes")
    plotty.plotdetectionquality()
    plotty.plotnblocksofbest()
    plotty.plotndecomps()

    #print instancenames
    #print blockcandidates
    #print blockcandidatesnvotes
    #print classnames
    #print classnmembers
    #print decompnblocks
    #print decompscores


    

if __name__ == "__main__": 
    main(sys.argv[1:])

