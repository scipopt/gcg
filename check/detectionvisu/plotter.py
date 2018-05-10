#!/usr/bin/python

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.rcsetup as rcsetup
from dataset import *

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
        print maxdetectiontime
        tauvals = np.arange(0, maxdetectiontime*1.1, 1.1*float(maxdetectiontime)/1000.)
        instfractsfordataset = []
        for dataset in datasets:
            instfract = []
            for tau in tauvals:
                instfract.append(dataset.fractionofmemberswithvalatmostwithscoreatleast(dataset.instancenames, dataset.detectiontimes, tau) )
            instfractsfordataset.append(instfract)
        plt.ylabel('fraction of instances')
        plt.xlabel('Detection time is at most')

        plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        plt.axis([0., maxdetectiontime*1.1, 0., 1.])
        for datasetid in range(len(datasets)):
            plt.plot(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='upper center',
           bbox_to_anchor=[0.5, 1.1],
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)

        plt.xscale('symlog')
        plt.show()




    def plotdetectionquality(self, datasets):
        tauvals = np.arange(0., 1., 0.01)
        instfractsfordataset = []
        labels = []
        for dataset in datasets:
            instfracts = []
            for tau in tauvals:
                instfracts.append(dataset.fractionofinstanceswithscoreatleast(dataset.decompscores, tau) )
            instfractsfordataset.append(instfracts)
        plt.ylabel('fraction of instances')
        plt.xlabel('Whitest found decomp has at least this max white score')
    #   print tauvals
        #print instancefractions
        plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        for datasetid in range(len(datasets)):
            plt.plot(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='upper center',
           bbox_to_anchor=[0.5, 1.1],
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)
        plt.show()


    def plotdetectionqualitysetpartmaster(self, datasets):
        tauvals = np.arange(0., 1., 0.01)
        instfractsfordataset = []
        labels = []
        for dataset in datasets:
            instfracts = []
            for tau in tauvals:
                instfracts.append(dataset.fractionofinstanceswithscoreatleastsetpartmaster(tau) )
            instfractsfordataset.append(instfracts)
        plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        plt.ylabel('fraction of instances')
        plt.xlabel('Whitest found decomp with setpartitioning master has at least this max white score')
    #   print tauvals
        #print instancefractions

        for datasetid in range(len(datasets)):
            plt.plot(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='upper center',
           bbox_to_anchor=[0.5, 1.1],
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)


        plt.show()



    def plotnblocksofbest(self, datasets):
        maxnblocks = 0
        for dataset in datasets:
            currblock = dataset.getmaxnblocks()
            if currblock > maxnblocks:
                maxnblocks = currblock
        tauvals = np.arange(0., maxnblocks)
        instfractsfordataset = []
        labels = []
        for dataset in datasets:
            instfracts = []
            for tau in tauvals:
                instfracts.append(dataset.fractionofinstanceswithnblocksleast(dataset.decompnblocks, tau) )
            instfractsfordataset.append(instfracts)
        plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        plt.ylabel('fraction of instances')
        plt.xlabel('whitest found decomposition has at least this number of blocks ')
    #   print tauvals
        #print instancefractions
        for datasetid in range(len(datasets)):
            plt.semilogx(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='upper center',
           bbox_to_anchor=[0.5, 1.1],
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)

        plt.show()

    def plotndecomps(self, datasets):
        maxndecomps = 0
        for dataset in datasets:
            currndecomps = dataset.getmaxnnontrivialdecomps()
            if currndecomps > maxndecomps:
                maxndecomps = currndecomps
        tauvals = np.arange(0., maxndecomps)
        instfractsfordataset = []
        labels = []
        for dataset in datasets:
            instfracts = []
            for tau in tauvals:
                instfracts.append(dataset.fractionofinstanceswithatleasttaunontrivialdecomps(dataset.decompscores, tau) )
            instfractsfordataset.append(instfracts)
        plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        plt.ylabel('fraction of instances')
        plt.xlabel('at least this number of nontrivial decompositions is found')

    #   print tauvals
        #print instancefractions

        for datasetid in range(len(datasets)):
            plt.plot(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='upper center',
           bbox_to_anchor=[0.5, 1.1],
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)

        plt.xscale('symlog')

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
        plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

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
