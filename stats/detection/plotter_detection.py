#!/usr/bin/env python3

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.rcsetup as rcsetup
import parser_detection as parser
import argparse


def parse_arguments(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', type=str,
                        default="plots",
                        help='output directory (default: "plots")')

    parser.add_argument('-c', '--classifier',
                        default="nonzeros",
                        help='classifier')

    parser.add_argument('filename', nargs='+',
                        help='.out-files to create plots with')
    args = parser.parse_args(args)
    return args

class Plotter:
    def __init__(self,fromApp=False):
        self.fromApp = fromApp

    def plotdetectiontimes(self, datasets, outdir="plots", filename="unknowntestset"):
        maxdetectiontime = 0.
        labels = []
        for dataset in datasets:
            currtime = dataset.getmaxdetectiontime()
            if currtime > maxdetectiontime:
                maxdetectiontime = currtime
        tauvals = np.arange(0, maxdetectiontime*1.1, 1.1*float(maxdetectiontime)/1000.)
        tauvals = np.insert(tauvals,len(tauvals),maxdetectiontime)
        instfractsfordataset = []
        for dataset in datasets:
            instfract = []
            for tau in tauvals:
                instfract.append(dataset.fractionofmemberswithvalatmostwithscoreatleast(dataset.instancenames, dataset.detectiontimes, tau) )
            instfractsfordataset.append(instfract)
        plt.ylabel('fraction of instances', size="small")
        plt.xlabel('Detection time is at most (seconds)', size="small")
        #plt.gca().set_prop_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        plt.axis([0., maxdetectiontime*1.1, 0., 1.])
        for datasetid in range(len(datasets)):
            plt.plot(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='lower left', bbox_to_anchor = (.0, 1.02, 1., 1.04), mode = 'expand', fontsize="small")

        if self.fromApp:
            plt.show()
        else:
            plt.savefig(os.path.join(outdir,'{}.detection.times.pdf'.format(filename)))
            plt.close()

    def plotdetectionquality(self, datasets, outdir="plots", filename="unknowntestset"):
        tauvals = np.arange(0., 1., 0.01)
        instfractsfordataset = []
        labels = []
        for dataset in datasets:
            instfracts = []
            for tau in tauvals:
                instfracts.append(dataset.fractionofinstanceswithscoreatleast(dataset.decompscores, tau) )
            instfractsfordataset.append(instfracts)
        plt.ylabel('fraction of instances', size="small")
        plt.xlabel('Whitest found decomp has at least this max white score', size="small")
        #plt.gca().set_prop_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        for datasetid in range(len(datasets)):
            plt.plot(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='lower left', bbox_to_anchor = (.0, 1.02, 1., 1.04), mode = 'expand', fontsize="small")

        if self.fromApp:
            plt.show()
        else:
            plt.savefig(os.path.join(outdir,'{}.detection.quality.pdf'.format(filename)))
            plt.close()

    def plotdetectionqualitysetpartmaster(self, datasets, outdir="plots", filename="unknowntestset"):
        tauvals = np.arange(0., 1., 0.01)
        instfractsfordataset = []
        labels = []
        for dataset in datasets:
            instfracts = []
            for tau in tauvals:
                instfracts.append(dataset.fractionofinstanceswithscoreatleastsetpartmaster(tau) )
            instfractsfordataset.append(instfracts)
        #plt.gca().set_prop_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        plt.ylabel('fraction of instances', size="small")
        plt.xlabel('Whitest found decomp by mastersetpart detector has at least this max white score', size="small")

        for datasetid in range(len(datasets)):
            plt.plot(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='lower left', bbox_to_anchor = (.0, 1.02, 1., 1.04), mode = 'expand', fontsize="small")

        if self.fromApp:
            plt.show()
        else:
            plt.savefig(os.path.join(outdir,'{}.detection.quality_SetPartMaster.pdf'.format(filename)))
            plt.close()

    def plotnblocksofbest(self, datasets, outdir="plots", filename="unknowntestset"):
        maxnblocks = 0
        for dataset in datasets:
            currblock = dataset.getmaxnblocks()
            if currblock > maxnblocks:
                maxnblocks = currblock
        tauvals = np.arange(0., maxnblocks)
        tauvals = np.insert(tauvals,len(tauvals),maxnblocks)
        instfractsfordataset = []
        labels = []
        for dataset in datasets:
            instfracts = []
            for tau in tauvals:
                instfracts.append(dataset.fractionofinstanceswithnblocksleast(dataset.decompnblocks, tau) )
            instfractsfordataset.append(instfracts)
        #plt.gca().set_prop_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        plt.ylabel('fraction of instances', size="small")
        plt.xlabel('whitest found decomposition has at least this number of blocks ', size="small")

        for datasetid in range(len(datasets)):
            plt.semilogx(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='lower left', bbox_to_anchor = (.0, 1.02, 1., 1.04), mode = 'expand', fontsize="small")

        if self.fromApp:
            plt.show()
        else:
            plt.savefig(os.path.join(outdir,'{}.detection.nBlocksOfBest.pdf'.format(filename)))
            plt.close()

    def plotndecomps(self, datasets, outdir="plots", filename="unknowntestset"):
        maxndecomps = 0
        for dataset in datasets:
            currndecomps = dataset.getmaxnnontrivialdecomps()
            if currndecomps > maxndecomps:
                maxndecomps = currndecomps
        tauvals = np.arange(0., maxndecomps)
        tauvals = np.insert(tauvals,len(tauvals),maxndecomps)
        instfractsfordataset = []
        labels = []
        for dataset in datasets:
            instfracts = []
            for tau in tauvals:
                instfracts.append(dataset.fractionofinstanceswithatleasttaunontrivialdecomps(dataset.decompscores, tau) )
            instfractsfordataset.append(instfracts)
        #plt.gca().set_prop_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        plt.ylabel('fraction of instances', size="small")
        plt.xlabel('at least this number of decompositions with score > 0 is found', size="small")

    #   print tauvals
        #print instancefractions

        for datasetid in range(len(datasets)):
            plt.semilogx(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='lower left', bbox_to_anchor = (.0, 1.02, 1., 1.04), mode = 'expand', fontsize="small")

        if self.fromApp:
            plt.show()
        else:
            plt.savefig(os.path.join(outdir,'{}.detection.decomps.pdf'.format(filename)))
            plt.close()

    def plotnclassesforclassifier(self, datasets, classifier, outdir="plots", filename="unknowntestset"):
        if self.fromApp:
            classifier = datasets[0].getclassifiernames()[classifier]
        maxnclasses = 0
        for dataset in datasets:
            for instance in dataset.instancenames:
                try:
                    if len(dataset.classnames[instance][classifier]) > maxnclasses:
                        maxnclasses = len(dataset.classnames[instance][classifier])
                except KeyError:
                    #print("Classifier {} did not work on instance {}. Skipping.".format(classifier, instance.split('.')[:-1]))
                    continue
        if maxnclasses == 0:
            print("Warning: No classifier worked, or data could not be read.")
        tauvals = np.arange(1., maxnclasses)
        tauvals = np.insert(tauvals,len(tauvals),maxnclasses)
        instfractsfordataset = []
        labels = []
        for dataset in datasets:
            instfracts = []
            for tau in tauvals:
                instfracts.append(dataset.fractionofinstanceswithatleasttaunontrivialdecomps(dataset.classnames, tau) )
            instfractsfordataset.append(instfracts)

        plt.ylabel('fraction of instances', size="small")
        plt.xlabel('at least this number of classes is found for classifier "'+str(classifier)+ '"', size="small")
        #plt.gca().set_prop_cycle(['red', 'green', 'blue', 'yellow', 'orange', 'pink', 'black', 'brown', 'magenta', 'purple', 'cyan', 'darkgreen'])

        for datasetid in range(len(datasets)):
            plt.semilogx(tauvals, instfractsfordataset[datasetid])
            labels.append(datasets[datasetid].getsettingsname())

        plt.legend(labels, ncol=4, loc='lower left', bbox_to_anchor = (.0, 1.02, 1., 1.04), mode = 'expand', fontsize="small")

        if self.fromApp:
            plt.show()
        else:
            plt.savefig(os.path.join(outdir,'{}.detection.classification_classes_{}.pdf'.format(filename,classifier)))
            plt.close()

def main():
    plotty = Plotter()
    datasets = []
    print("Parsing files...")
    args = sys.argv[1:]
    args = parse_arguments(args)

    for outfile in args.filename:
        datasets.append(parser.Dataset(outfile) )
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # take first testset's name
    args.filename = args.filename[0].split('/')[-1].split('.')[1]
    print("Generating visualizations...")
    plotty.plotdetectiontimes(datasets, outdir=args.outdir, filename=args.filename)
    plotty.plotdetectionquality(datasets, outdir=args.outdir, filename=args.filename)
    plotty.plotdetectionqualitysetpartmaster(datasets, outdir=args.outdir, filename=args.filename)
    plotty.plotnblocksofbest(datasets, outdir=args.outdir, filename=args.filename)
    plotty.plotndecomps(datasets, outdir=args.outdir, filename=args.filename)
    plotty.plotnclassesforclassifier(datasets,str(args.classifier), outdir=args.outdir, filename=args.filename)


if __name__ == "__main__":
    main()
