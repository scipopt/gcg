#!/usr/bin/env python3
# This script reads *.out files from whole testsets created by GCG's make test
# and parses its statistics into a pandas dataframe.
# Very hacky, beware!
import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# help functions
def ct(string):
    sim = string.strip()[1:-1]
    return sim.upper()

def real_path(string):
    splitted = string.split('/')
    check = False
    ret = ""
    for w in splitted:
        if w == "check":
            check = True
        elif check:
            ret += "/" + w
        else:
            continue
    return ret[1:]

def parseOutfiles(outfiles):
    # main data dictionary. Will contain data for ALL outfiles
    d = {
        'TOTAL TIME': [],
        'READING TIME': [],
        'COPYING TIME': [],
        'DETECTION TIME': [],
        'PRESOLVING TIME': [],
        'STATUS': [],
        'ROOT NODE TIME': [],
        'DUAL BOUNDS': [],
        'HEUR TIME ORIG': [],
        'HEUR CALLS ORIG': [],
        'HEUR FOUND ORIG': [],
        'HEUR TIME MASTER': [],
        'HEUR CALLS MASTER': [],
        'HEUR FOUND MASTER': [],
        'CUTS TIME MASTER': [],
        'CUTS CALLS MASTER': [],
        'CUTS FOUND MASTER': [],
        'CUTS APPLIED MASTER': [],
        'CUTS TIME ORIG': [],
        'CUTS CALLS ORIG': [],
        'CUTS FOUND ORIG': [],
        'CUTS APPLIED ORIG': [],
        'FARKAS TIME': [],
        'MASTER TIME': [],
        'PRICING TIME': [],
        'PRICING SOLVER TIME': [],
        'PRICING SOLVER TYPE': [],
        'DEGENERACY': [],
        'CONS LINEAR': [],
        'CONS KNAPSACK': [],
        'CONS LOGICOR': [],
        'CONS SETPPC': [],
        'CONS VARBOUND': [],
        'CONS AND': [],
        'LINKING VARS': [],
        'NBLOCKS': [],
        'NBLOCKSAGGR': [],
        'SOLUTIONS FOUND': [],
        'FIRST SOLUTION TIME': [],
        'BEST SOLUTION TIME': [],
        'PD INTEGRAL': [],
        'MASTER NCONSS': [],
        'MASTER NVARS': [],
        'BNB TREE NODES': [],
        'BNB TREE LEFT': [],
        'BNB TREE DEPTH': [],
        'BR RULE TIME ORIG': [],
        'BR RULE TIME GENERIC': [],
        'BR RULE TIME RELPSPROB': [],
        'BR RULE TIME RYANFOSTER': [],
        'BR RULE CALLS ORIG': [],
        'BR RULE CALLS GENERIC': [],
        'BR RULE CALLS RELPSPROB': [],
        'BR RULE CALLS RYANFOSTER': [],
        'RMP LP CALLS': [],
        'RMP LP TIME': [],
        'RMP LP ITERATIONS': [],
        'ORIGINAL LP CALLS': [],
        'ORIGINAL LP TIME': [],
        'ORIGINAL LP ITERATIONS': [],
        'LP FILE': [],
        'DEC FILE': [],
    }
    # "temporary" data dictionary. Will contain data for a single outfile and is reset afterwards
    data = {
        'TOTAL TIME': [],
        'READING TIME': [],
        'COPYING TIME': [],
        'DETECTION TIME': [],
        'PRESOLVING TIME': [],
        'STATUS': [],
        'ROOT NODE TIME': [],
        'DUAL BOUNDS': [],
        'HEUR TIME ORIG': [],
        'HEUR CALLS ORIG': [],
        'HEUR FOUND ORIG': [],
        'HEUR TIME MASTER': [],
        'HEUR CALLS MASTER': [],
        'HEUR FOUND MASTER': [],
        'CUTS TIME MASTER': [],
        'CUTS CALLS MASTER': [],
        'CUTS FOUND MASTER': [],
        'CUTS APPLIED MASTER': [],
        'CUTS TIME ORIG': [],
        'CUTS CALLS ORIG': [],
        'CUTS FOUND ORIG': [],
        'CUTS APPLIED ORIG': [],
        'FARKAS TIME': [],
        'MASTER TIME': [],
        'PRICING TIME': [],
        'PRICING SOLVER TIME': [],
        'PRICING SOLVER TYPE': [],
        'DEGENERACY': [],
        'CONS LINEAR': [],
        'CONS KNAPSACK': [],
        'CONS LOGICOR': [],
        'CONS SETPPC': [],
        'CONS VARBOUND': [],
        'CONS AND': [],
        'LINKING VARS': [],
        'NBLOCKS': [],
        'NBLOCKSAGGR': [],
        'SOLUTIONS FOUND': [],
        'FIRST SOLUTION TIME': [],
        'BEST SOLUTION TIME': [],
        'PD INTEGRAL': [],
        'MASTER NCONSS': [],
        'MASTER NVARS': [],
        'BNB TREE NODES': [],
        'BNB TREE LEFT': [],
        'BNB TREE DEPTH': [],
        'BR RULE TIME ORIG': [],
        'BR RULE TIME GENERIC': [],
        'BR RULE TIME RELPSPROB': [],
        'BR RULE TIME RYANFOSTER': [],
        'BR RULE CALLS ORIG': [],
        'BR RULE CALLS GENERIC': [],
        'BR RULE CALLS RELPSPROB': [],
        'BR RULE CALLS RYANFOSTER': [],
        'RMP LP CALLS': [],
        'RMP LP TIME': [],
        'RMP LP ITERATIONS': [],
        'ORIGINAL LP CALLS': [],
        'ORIGINAL LP TIME': [],
        'ORIGINAL LP ITERATIONS': [],
        'LP FILE': [],
        'DEC FILE': [],
    }

    # instance names
    index = []
    idx = []

    # initialize variables
    search = ""
    it = 0
    opstat = False
    ot = False
    status = False
    presolved = False
    read = False
    presolve = False
    copying = False
    pricersdone = False
    pricingsolversdone = False

    # write in dictionary
    for outfile in outfiles:
        SCIPlog = False
        #print(outfile)
        fh = open(outfile, 'r')
        for line in fh:
            # get instance name by read problem line
            #if line.startswith("read problem") and '.dec' not in line and '.blk' not in line:
            #    print(line.split()[2].split('/')[-1].split('.')[0])
            #    index.append(line.split()[2].split('/')[-1].split('.')[0])

            # get instance name by @01 tag (made by make test script)
            if line.startswith("@01"):
                #print(line.split()[1])
                index.append(line.split()[1])

            if not SCIPlog and line.startswith("SCIP>"):
                #print("Interpreting as SCIP (non-GCG) log!")
                SCIPlog = True
                opstat = True

            # get instance lp
            if line.startswith("read problem") and '.dec' not in line and '.blk' not in line:
                if line.split()[2][1:-1].startswith("/") and "/check/" in line.split()[2][1:-1]:
                    data['LP FILE'].append(line.split()[2][1:-1].split("/check/")[1])
                else:
                    data['LP FILE'].append(line.split()[2][1:-1])

            # get instance dec
            if line.startswith("read problem") and ('.dec' in line or '.blk' in line) and not SCIPlog:
                if line.split()[2][1:-1].startswith("/") and "/check/" in line.split()[2][1:-1]:
                    data['DEC FILE'].append(line.split()[2][1:-1].split("/check/")[1])
                else:
                    data['DEC FILE'].append(line.split()[2][1:-1])

            if line.startswith("Detection Time: "):
                data['DETECTION TIME'].append(float(line.split(':')[1].strip()))

            # reading of master stats finished
            if line.startswith("Original Program statistics:"):
                opstat = True
                continue
            elif line.startswith("Original Program Solution statistics:"):
                data['TOTAL TIME'].append(0.)
                ot = True
                data['READING TIME'].append(0.)
                read = True
                data['PRESOLVING TIME'].append(0.)
                data['COPYING TIME'].append(0.)
                copying = True
                opstat = True

            # get TOTAL TIME
            if line.startswith("Total Time         :") and opstat:
                data['TOTAL TIME'].append(float(line.split(':')[1]))
                ot = True
                continue

            # get status
            if line.startswith("SCIP Status") and not status:
                if line.split(':')[1].strip() == "problem is solved [optimal solution found]":
                    data['STATUS'].append(1)
                elif line.split(':')[1].strip() == "problem is solved [infeasible]":
                    data['STATUS'].append(2)
                elif line.split(':')[1].strip() == "solving was interrupted [time limit reached]":
                    data['STATUS'].append(3)
                elif line.split(':')[1].strip() == "solving was interrupted [memory limit reached]":
                    data['STATUS'].append(4)
                elif line.split(':')[1].strip() == "solving was interrupted [node limit reached]":
                    data['STATUS'].append(5)
                else:
                    data['STATUS'].append(0)
                status = True


            # get root node time
            if line.startswith("Time in root node:") or line.startswith("  time in root node:"):
                data['ROOT NODE TIME'].append(float(line.split(':')[1]))
                continue

            # get degeneracy
            if line.startswith("Degeneracy:"):
                search = "DEGENERACY"
                data['DEGENERACY'].append([])
                continue

            if search == "DEGENERACY":
                if line.startswith("Dual Bounds:"):
                    search = "DUALS"    # no empty string here!
                    data['DUAL BOUNDS'].append([])
                    continue
                data['DEGENERACY'][-1].append((int(line.split(':')[0]), float(line.split(':')[1])))
                continue

            # get dual bound development
            if line.startswith("Dual Bounds:"):
                search = "DUALS"
                data['DUAL BOUNDS'].append([])
                continue

            if search == "DUALS":
                if line.startswith("GCG"):
                    search = ""
                    continue
                data['DUAL BOUNDS'][-1].append((int(line.split(':')[0]), float(line.split(':')[1])))
                continue

            # get successful heuristics (master)
            if line.startswith("Primal Heuristics") and opstat:
                search = "HEURISTICS MASTER"
                data['HEUR TIME MASTER'].append(0.)
                data['HEUR CALLS MASTER'].append(0)
                data['HEUR FOUND MASTER'].append(0)
                continue

            if search == "HEURISTICS MASTER":
                if not line.startswith("Diving Statistics"):
                    if line.split(':')[1].split()[0].replace('.', '', 1).isdigit():
                        data['HEUR TIME MASTER'][-1] += float(line.split(':')[1].split()[0])
                        if line.split(':')[1].split()[1].replace('.', '', 1).isdigit():
                            data['HEUR TIME MASTER'][-1] += float(line.split(':')[1].split()[1])
                        if line.split(':')[1].split()[2].isdigit():
                            data['HEUR CALLS MASTER'][-1] += int(line.split(':')[1].split()[2])
                        if line.split(':')[1].split()[3].isdigit():
                            data['HEUR FOUND MASTER'][-1] += int(line.split(':')[1].split()[3])
                    else:
                        search = ""

            # get successful heuristics (master)
            if line.startswith("Primal Heuristics"):
                search = "HEURISTICS ORIG"
                data['HEUR TIME ORIG'].append(0.)
                data['HEUR CALLS ORIG'].append(0)
                data['HEUR FOUND ORIG'].append(0)
                continue

            if search == "HEURISTICS ORIG":
                if not line.startswith("Diving Statistics"):
                    if line.split(':')[1].split()[0].replace('.', '', 1).isdigit():
                        data['HEUR TIME ORIG'][-1] += float(line.split(':')[1].split()[0])
                        if line.split(':')[1].split()[1].replace('.', '', 1).isdigit():
                            data['HEUR TIME ORIG'][-1] += float(line.split(':')[1].split()[1])
                        if line.split(':')[1].split()[2].isdigit():
                            data['HEUR CALLS ORIG'][-1] += int(line.split(':')[1].split()[2])
                        if line.split(':')[1].split()[3].isdigit():
                            data['HEUR FOUND ORIG'][-1] += int(line.split(':')[1].split()[3])
                    else:
                        search = ""

            # get branching rule statistics
            if line.startswith("Branching Rules") and opstat == False:
                search = "BRANCHINGRULES"
                #data['BR RULE TIME EMPTY'].append(0.)
                data['BR RULE TIME GENERIC'].append(0.)
                data['BR RULE TIME ORIG'].append(0.)
                data['BR RULE TIME RELPSPROB'].append(0.)
                data['BR RULE TIME RYANFOSTER'].append(0.)
                # calls (lp, ext and ps)
                #data['BR RULE CALLS EMPTY'].append(0)
                data['BR RULE CALLS GENERIC'].append(0)
                data['BR RULE CALLS ORIG'].append(0)
                data['BR RULE CALLS RELPSPROB'].append(0)
                data['BR RULE CALLS RYANFOSTER'].append(0)
                continue

            if search == "BRANCHINGRULES":
                if not line.startswith("Primal Heuristics"):
                    if line.split(':')[1].split()[0].replace('.', '', 1).isdigit():
                        if line.lstrip().startswith("generic"):
                            data['BR RULE TIME GENERIC'][-1] += float(line.split(':')[1].split()[0])
                            data['BR RULE TIME GENERIC'][-1] += float(line.split(':')[1].split()[1])
                            data['BR RULE CALLS GENERIC'][-1] += int(line.split(':')[1].split()[2])
                            data['BR RULE CALLS GENERIC'][-1] += int(line.split(':')[1].split()[3])
                            data['BR RULE CALLS GENERIC'][-1] += int(line.split(':')[1].split()[4])
                        elif line.lstrip().startswith("orig"):
                            data['BR RULE TIME ORIG'][-1] += float(line.split(':')[1].split()[0])
                            data['BR RULE TIME ORIG'][-1] += float(line.split(':')[1].split()[1])
                            data['BR RULE CALLS ORIG'][-1] += int(line.split(':')[1].split()[2])
                            data['BR RULE CALLS ORIG'][-1] += int(line.split(':')[1].split()[3])
                            data['BR RULE CALLS ORIG'][-1] += int(line.split(':')[1].split()[4])
                        elif line.lstrip().startswith("relpsprob"):
                            data['BR RULE TIME RELPSPROB'][-1] += float(line.split(':')[1].split()[0])
                            data['BR RULE TIME RELPSPROB'][-1] += float(line.split(':')[1].split()[1])
                            data['BR RULE CALLS RELPSPROB'][-1] += int(line.split(':')[1].split()[2])
                            data['BR RULE CALLS RELPSPROB'][-1] += int(line.split(':')[1].split()[3])
                            data['BR RULE CALLS RELPSPROB'][-1] += int(line.split(':')[1].split()[4])
                        elif line.lstrip().startswith("ryanfoster"):
                            data['BR RULE TIME RYANFOSTER'][-1] += float(line.split(':')[1].split()[0])
                            data['BR RULE TIME RYANFOSTER'][-1] += float(line.split(':')[1].split()[1])
                            data['BR RULE CALLS RYANFOSTER'][-1] += int(line.split(':')[1].split()[2])
                            data['BR RULE CALLS RYANFOSTER'][-1] += int(line.split(':')[1].split()[3])
                            data['BR RULE CALLS RYANFOSTER'][-1] += int(line.split(':')[1].split()[4])
                    else:
                        search = ""

            # get cutting plane statistics
            if line.startswith("Separators") and opstat:
                search = "CUTS MASTER"
                data['CUTS TIME MASTER'].append(0.)
                data['CUTS CALLS MASTER'].append(0)
                data['CUTS FOUND MASTER'].append(0)
                data['CUTS APPLIED MASTER'].append(0)
                continue

            if search == "CUTS MASTER":
                if not line.startswith("Pricers") and not line.startswith("Cutselectors"):
                    offset = 0 if line.split(':')[0].strip() != "cut pool" else -1
                    if line.split(':')[1].split()[0].isdigit():
                        data['CUTS TIME MASTER'][-1] += float(line.split(':')[1].split()[0])
                    if line.split(':')[1].split()[2+offset].isdigit():
                        data['CUTS CALLS MASTER'][-1] += int(line.split(':')[1].split()[2+offset])
                    if line.split(':')[1].split()[5+offset].isdigit():
                        data['CUTS FOUND MASTER'][-1] += int(line.split(':')[1].split()[5+offset])
                    if line.split(':')[1].split()[6+offset].isdigit():
                        data['CUTS APPLIED MASTER'][-1] += int(line.split(':')[1].split()[6+offset])
                else:
                    search = ""

            # get cutting plane statistics
            if line.startswith("Separators") and not opstat:
                search = "CUTS ORIG"
                data['CUTS TIME ORIG'].append(0.)
                data['CUTS CALLS ORIG'].append(0)
                data['CUTS FOUND ORIG'].append(0)
                data['CUTS APPLIED ORIG'].append(0)
                continue

            if search == "CUTS ORIG":
                if not line.startswith("Pricers") and not line.startswith("Cutselectors"):
                    offset = 0 if line.split(':')[0].strip() != "cut pool" else -1
                    if line.split(':')[1].split()[0].isdigit():
                        data['CUTS TIME ORIG'][-1] += float(line.split(':')[1].split()[0])
                    if line.split(':')[1].split()[2+offset].isdigit():
                        data['CUTS CALLS ORIG'][-1] += int(line.split(':')[1].split()[2+offset])
                    if line.split(':')[1].split()[5+offset].isdigit():
                        data['CUTS FOUND ORIG'][-1] += int(line.split(':')[1].split()[5+offset])
                    if line.split(':')[1].split()[6+offset].isdigit():
                        data['CUTS APPLIED ORIG'][-1] += int(line.split(':')[1].split()[6+offset])
                else:
                    search = ""

            if search == "PRICING SOLVER":
                if line.lstrip().startswith("Solving Details"):
                    search = ""
                elif sum([int(x) for x in line.split(':')[1].split()[:4]]) > 0.02:
                    if line.lstrip().startswith("knapsack"):
                        # type: Knapsack (=1)
                        data['PRICING SOLVER TYPE'][-1].append("Knapsack")
                        data['FARKAS TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:6]])
                        data['PRICING SOLVER TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:]])
                    elif line.lstrip().startswith("cliquer"):
                        # type: Cliquer (=2)
                        data['PRICING SOLVER TYPE'][-1].append("Cliquer")
                        data['FARKAS TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:6]])
                        data['PRICING SOLVER TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:]])
                    elif line.lstrip().startswith("mip"):
                        # type: CLIQUER (=4)
                        data['PRICING SOLVER TYPE'][-1].append("MIP")
                        data['FARKAS TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:6]])
                        data['PRICING SOLVER TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:]])
                    else:
                        print(line)
                        try:
                            # type: own solver (=8)
                            data['PRICING SOLVER TYPE'][-1].append("Custom")
                            data['FARKAS TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:6]])
                            data['PRICING SOLVER TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:]])
                        except:
                            if line.startswith("SCIP Status"):
                                search = ""
                            else:
                                if (len(str(line.strip())) != 0): print(f"Unable to handle pricing solver '{line.lstrip().split(':')[0]}'")

            # get Farkas Time and type of Pricing (Cliquer / Knapsack)
            if line.startswith("Pricing Solver"):
                search = "PRICING SOLVER"
                # initialize values and say that pricing solvers are done because we now collect and are done afterwards
                if not pricingsolversdone:
                    # append 0 for pricing solver type because pricing took place (else there were no pricing solver section),
                    # but not neccessarily pricing solvers were needed
                    data['PRICING SOLVER TYPE'].append([])
                    data['FARKAS TIME'].append(0.)
                    data['PRICING SOLVER TIME'].append(0.)
                pricingsolversdone = True

            if search == "PRICING":
                if line.lstrip().startswith("problem variables") or line.lstrip().startswith("gcg"):
                    data['PRICING TIME'][-1] += float(line.split(':')[1].split()[0])
                else:
                    search = ""

            if line.startswith("Pricers"):
                search = "PRICING"
                if not pricersdone:
                    data['PRICING TIME'].append(0.)
                pricersdone = True

            # get Master time
            if line.startswith("Master Program statistics:"):
                search = "MASTER"

            if search == "MASTER":
                try:
                    if line.split(':')[1].strip() == "problem creation / modification":
                        #then there will be no master solving time, thus we set it to 0
                        #print("No master time found for instance %s" % index[-1])
                        data['MASTER TIME'].append(0.)
                        search = ""
                except:
                    pass

            if search == "MASTER" and line.split(':')[0].strip() == "solving":
                data['MASTER TIME'].append(float(line.split(':')[1]))
                search = ""

            # get reading time
            if line.startswith("  reading") and not read and opstat:
                data['READING TIME'].append(float(line.split(':')[1]))
                read = True
                continue

            # get presolving time
            if line.startswith("  presolving") and not presolve and opstat:
                line = line.split('(')[0]
                data['PRESOLVING TIME'].append(float(line.split(':')[1]))
                presolve = True
                continue

            # get copying time
            if line.startswith("  copying") and not copying and opstat:
                line = line.split('(')[0]
                data['COPYING TIME'].append(float(line.split(':')[1]))
                copying = True
                continue

            # get constraints
            if line.startswith("presolved problem has") and not presolved and not SCIPlog:
                search = "CONSS"
                data['CONS LINEAR'].append(0)
                data['CONS KNAPSACK'].append(0)
                data['CONS LOGICOR'].append(0)
                data['CONS SETPPC'].append(0)
                data['CONS VARBOUND'].append(0)
                data['CONS AND'].append(0)
                continue

            if search == "CONSS":
                res = re.search("constraints of type", line)
                if res:
                    constype = ct(line[res.end():-1])
                    data['CONS ' + constype][-1] = int(line[:res.start()])
                else:
                    search = ""
                    presolved = True

            # get number of blocks
            if line.startswith("Decomp statistics"):
                search = "BLOCKS"

            if search == "BLOCKS":
                if line.lstrip().startswith("blocks"):
                    data['NBLOCKS'].append(int(line.split(':')[1]))
                if line.lstrip().startswith("aggr. blocks"):
                    data['NBLOCKSAGGR'].append(int(line.split(':')[1]))
                    search = ""

            # get solution statistics
            if line.startswith("Solution") and not opstat:
                search = "SOLUTION"
                continue

            if search == "SOLUTION":
                if line.lstrip().startswith("Solutions found"):
                    data['SOLUTIONS FOUND'].append(int(line.split(':')[1].split()[0]))
                elif line.lstrip().startswith("First Solution"):
                    data['FIRST SOLUTION TIME'].append(float(line.split(':')[1].split()[7]))
                elif line.lstrip().startswith("Primal Bound") and data['SOLUTIONS FOUND'][-1] > 0:
                    data['BEST SOLUTION TIME'].append(float(line.split(':')[1].split()[7]))
                elif line.lstrip().startswith("Avg. Gap") and data['SOLUTIONS FOUND'][-1] > 0:
                    data['PD INTEGRAL'].append(float(line.split(':')[1].split('%')[1].split()[0][1:]))
                    search = ""
                else:
                    continue

            # get master statistics
            if line.startswith("Master statistics"):
                search = "MASTER STATS"
                continue

            if search == "MASTER STATS":
                if line.lstrip().startswith("master"):
                    data['MASTER NCONSS'].append(int(line.split(':')[1].split()[6]))
                    data['MASTER NVARS'].append(int(line.split(':')[1].split()[0]))
                    search = ""

            if line.startswith("Number of LinkingVars:"):
                data['LINKING VARS'].append(int(line.split(':')[1]))

            # get Branch-and-Bound Tree stats
            if line.startswith("B&B Tree") and opstat:
                search = "BNB"
                continue

            if search == "BNB":
                if line.lstrip().startswith("nodes (total)"):
                    data['BNB TREE NODES'].append(int(line.split(':')[1].split()[0]))
                elif line.lstrip().startswith("nodes left"):
                    data['BNB TREE LEFT'].append(int(line.split(':')[1]))
                elif line.lstrip().startswith("max depth (total)"):
                    data['BNB TREE DEPTH'].append(int(line.split(':')[1]))

            # get LP stats
            if line.startswith("LP") and not opstat:
                search = "RMP LP"
                continue

            if search == "RMP LP":
                if line.lstrip().startswith("primal LP"):
                    data['RMP LP CALLS'].append(int(line.split(':')[1].split()[1]))
                    data['RMP LP TIME'].append(float(line.split(':')[1].split()[0]))
                    data['RMP LP ITERATIONS'].append(int(line.split(':')[1].split()[2]))
                elif line.lstrip().startswith("dual LP"):
                    data['RMP LP CALLS'][-1] += int(line.split(':')[1].split()[1])
                    data['RMP LP TIME'][-1] += float(line.split(':')[1].split()[0])
                    data['RMP LP ITERATIONS'][-1] += int(line.split(':')[1].split()[2])
                elif line.lstrip().startswith("lex dual LP"):
                    data['RMP LP CALLS'][-1] += int(line.split(':')[1].split()[1])
                    data['RMP LP TIME'][-1] += float(line.split(':')[1].split()[0])
                    data['RMP LP ITERATIONS'][-1] += int(line.split(':')[1].split()[2])
                elif line.lstrip().startswith("barrier LP"):
                    data['RMP LP CALLS'][-1] += int(line.split(':')[1].split()[1])
                    data['RMP LP TIME'][-1] += float(line.split(':')[1].split()[0])
                    data['RMP LP ITERATIONS'][-1] += int(line.split(':')[1].split()[2])
                elif line.lstrip().startswith("resolve instable"):
                    data['RMP LP CALLS'][-1] += int(line.split(':')[1].split()[1])
                    data['RMP LP TIME'][-1] += float(line.split(':')[1].split()[0])
                    data['RMP LP ITERATIONS'][-1] += int(line.split(':')[1].split()[2])

            # get LP stats
            if line.startswith("LP") and opstat:
                search = "ORIGINAL LP"
                continue

            if search == "ORIGINAL LP":
                if line.lstrip().startswith("primal LP"):
                    data['ORIGINAL LP CALLS'].append(int(line.split(':')[1].split()[1]))
                    data['ORIGINAL LP TIME'].append(float(line.split(':')[1].split()[0]))
                    data['ORIGINAL LP ITERATIONS'].append(int(line.split(':')[1].split()[2]))
                elif line.lstrip().startswith("dual LP"):
                    data['ORIGINAL LP CALLS'][-1] += int(line.split(':')[1].split()[1])
                    data['ORIGINAL LP TIME'][-1] += float(line.split(':')[1].split()[0])
                    data['ORIGINAL LP ITERATIONS'][-1] += int(line.split(':')[1].split()[2])
                elif line.lstrip().startswith("lex dual LP"):
                    data['ORIGINAL LP CALLS'][-1] += int(line.split(':')[1].split()[1])
                    data['ORIGINAL LP TIME'][-1] += float(line.split(':')[1].split()[0])
                    data['ORIGINAL LP ITERATIONS'][-1] += int(line.split(':')[1].split()[2])
                elif line.lstrip().startswith("barrier LP"):
                    data['ORIGINAL LP CALLS'][-1] += int(line.split(':')[1].split()[1])
                    data['ORIGINAL LP TIME'][-1] += float(line.split(':')[1].split()[0])
                    data['ORIGINAL LP ITERATIONS'][-1] += int(line.split(':')[1].split()[2])
                elif line.lstrip().startswith("resolve instable"):
                    data['ORIGINAL LP CALLS'][-1] += int(line.split(':')[1].split()[1])
                    data['ORIGINAL LP TIME'][-1] += float(line.split(':')[1].split()[0])
                    data['ORIGINAL LP ITERATIONS'][-1] += int(line.split(':')[1].split()[2])


            # sync point
            if line.startswith("=ready="):
                it += 1
                search = ""
                if SCIPlog: opstat = True
                else: opstat = False
                ot = False
                read = False
                status = False
                presolved = False
                presolve = False
                copying = False
                pricersdone = False
                pricingsolversdone = False


                if len(data['TOTAL TIME']) < it:
                    data['TOTAL TIME'].append(float('NaN'))
                if len(data['READING TIME']) < it:
                    data['READING TIME'].append(float('NaN'))
                if len(data['PRESOLVING TIME']) < it:
                    data['PRESOLVING TIME'].append(float('NaN'))
                if len(data['COPYING TIME']) < it:
                    data['COPYING TIME'].append(float('NaN'))
                if len(data['DETECTION TIME']) < it:
                    data['DETECTION TIME'].append(float('NaN'))
                if len(data['STATUS']) < it:
                    data['STATUS'].append(0)
                if len(data['ROOT NODE TIME']) < it:
                    data['ROOT NODE TIME'].append(float('NaN'))

                if len(data['DUAL BOUNDS']) < it:
                    data['DUAL BOUNDS'].append([float('NaN')])

                if len(data['HEUR TIME MASTER']) < it:
                    data['HEUR TIME MASTER'].append(float('NaN'))
                elif len(data['HEUR TIME MASTER']) == it + 1:
                    data['HEUR TIME MASTER'][-2] += data['HEUR TIME MASTER'][-1]
                    data['HEUR TIME MASTER'] = data['HEUR TIME MASTER'][:-1]
                if len(data['HEUR CALLS MASTER']) < it:
                    data['HEUR CALLS MASTER'].append(-1)
                elif len(data['HEUR CALLS MASTER']) == it + 1:
                    data['HEUR CALLS MASTER'][-2] += data['HEUR CALLS MASTER'][-1]
                    data['HEUR CALLS MASTER'] = data['HEUR CALLS MASTER'][:-1]
                if len(data['HEUR FOUND MASTER']) < it:
                    data['HEUR FOUND MASTER'].append(-1)
                elif len(data['HEUR FOUND MASTER']) == it + 1:
                    data['HEUR FOUND MASTER'][-2] += data['HEUR FOUND MASTER'][-1]
                    data['HEUR FOUND MASTER'] = data['HEUR FOUND MASTER'][:-1]

                if len(data['HEUR TIME ORIG']) < it:
                    data['HEUR TIME ORIG'].append(float('NaN'))
                elif len(data['HEUR TIME ORIG']) == it + 1:
                    data['HEUR TIME ORIG'][-2] += data['HEUR TIME ORIG'][-1]
                    data['HEUR TIME ORIG'] = data['HEUR TIME ORIG'][:-1]
                if len(data['HEUR CALLS ORIG']) < it:
                    data['HEUR CALLS ORIG'].append(-1)
                elif len(data['HEUR CALLS ORIG']) == it + 1:
                    data['HEUR CALLS ORIG'][-2] += data['HEUR CALLS ORIG'][-1]
                    data['HEUR CALLS ORIG'] = data['HEUR CALLS ORIG'][:-1]
                if len(data['HEUR FOUND ORIG']) < it:
                    data['HEUR FOUND ORIG'].append(-1)
                elif len(data['HEUR FOUND ORIG']) == it + 1:
                    data['HEUR FOUND ORIG'][-2] += data['HEUR FOUND ORIG'][-1]
                    data['HEUR FOUND ORIG'] = data['HEUR FOUND ORIG'][:-1]

                if len(data['CUTS TIME MASTER']) < it:
                    data['CUTS TIME MASTER'].append(float('NaN'))
                elif len(data['CUTS TIME MASTER']) == it + 1:
                    data['CUTS TIME MASTER'][-2] += data['CUTS TIME MASTER'][-1]
                    data['CUTS TIME MASTER'] = data['CUTS TIME MASTER'][:-1]
                if len(data['CUTS CALLS MASTER']) < it:
                    data['CUTS CALLS MASTER'].append(-1)
                elif len(data['CUTS CALLS MASTER']) == it + 1:
                    data['CUTS CALLS MASTER'][-2] += data['CUTS CALLS MASTER'][-1]
                    data['CUTS CALLS MASTER'] = data['CUTS CALLS MASTER'][:-1]
                if len(data['CUTS FOUND MASTER']) < it:
                    data['CUTS FOUND MASTER'].append(-1)
                elif len(data['CUTS FOUND MASTER']) == it + 1:
                    data['CUTS FOUND MASTER'][-2] += data['CUTS FOUND MASTER'][-1]
                    data['CUTS FOUND MASTER'] = data['CUTS FOUND MASTER'][:-1]
                if len(data['CUTS APPLIED MASTER']) < it:
                    data['CUTS APPLIED MASTER'].append(-1)
                elif len(data['CUTS APPLIED MASTER']) == it + 1:
                    data['CUTS APPLIED MASTER'][-2] += data['CUTS APPLIED MASTER'][-1]
                    data['CUTS APPLIED MASTER'] = data['CUTS APPLIED MASTER'][:-1]

                if len(data['CUTS TIME ORIG']) < it:
                    data['CUTS TIME ORIG'].append(float('NaN'))
                elif len(data['CUTS TIME ORIG']) == it + 1:
                    data['CUTS TIME ORIG'][-2] += data['CUTS TIME'][-1]
                    data['CUTS TIME ORIG'] = data['CUTS TIME ORIG'][:-1]
                if len(data['CUTS CALLS ORIG']) < it:
                    data['CUTS CALLS ORIG'].append(-1)
                elif len(data['CUTS CALLS ORIG']) == it + 1:
                    data['CUTS CALLS ORIG'][-2] += data['CUTS CALLS ORIG'][-1]
                    data['CUTS CALLS ORIG'] = data['CUTS CALLS ORIG'][:-1]
                if len(data['CUTS FOUND ORIG']) < it:
                    data['CUTS FOUND ORIG'].append(-1)
                elif len(data['CUTS FOUND ORIG']) == it + 1:
                    data['CUTS FOUND ORIG'][-2] += data['CUTS FOUND ORIG'][-1]
                    data['CUTS FOUND ORIG'] = data['CUTS FOUND ORIG'][:-1]
                if len(data['CUTS APPLIED ORIG']) < it:
                    data['CUTS APPLIED ORIG'].append(-1)
                elif len(data['CUTS APPLIED ORIG']) == it + 1:
                    data['CUTS APPLIED ORIG'][-2] += data['CUTS APPLIED ORIG'][-1]
                    data['CUTS APPLIED ORIG'] = data['CUTS APPLIED ORIG'][:-1]

                if len(data['FARKAS TIME']) < it:
                    data['FARKAS TIME'].append(float('NaN'))
                if len(data['MASTER TIME']) < it:
                    data['MASTER TIME'].append(float('NaN'))
                if len(data['PRICING TIME']) < it:
                    data['PRICING TIME'].append(float('NaN'))
                if len(data['PRICING SOLVER TIME']) < it:
                    data['PRICING SOLVER TIME'].append(float('NaN'))

                if len(data['PRICING SOLVER TYPE']) < it:
                    data['PRICING SOLVER TYPE'].append(-1)

                if len(data['DEGENERACY']) < it:
                    data['DEGENERACY'].append(float('NaN'))

                if len(data['CONS LINEAR']) < it:
                    data['CONS LINEAR'].append(-1)
                if len(data['CONS KNAPSACK']) < it:
                    data['CONS KNAPSACK'].append(-1)
                if len(data['CONS LOGICOR']) < it:
                    data['CONS LOGICOR'].append(-1)
                if len(data['CONS SETPPC']) < it:
                    data['CONS SETPPC'].append(-1)
                if len(data['CONS VARBOUND']) < it:
                    data['CONS VARBOUND'].append(-1)
                if len(data['CONS AND']) < it:
                    data['CONS AND'].append(-1)

                if len(data['NBLOCKS']) < it:
                    data['NBLOCKS'].append(-1)

                if len(data['NBLOCKSAGGR']) < it:
                    data['NBLOCKSAGGR'].append(-1)

                if len(data['SOLUTIONS FOUND']) < it:
                    data['SOLUTIONS FOUND'].append(-1)
                if len(data['FIRST SOLUTION TIME']) < it:
                    data['FIRST SOLUTION TIME'].append(float('NaN'))
                if len(data['BEST SOLUTION TIME']) < it:
                    data['BEST SOLUTION TIME'].append(float('NaN'))
                if len(data['PD INTEGRAL']) < it:
                    data['PD INTEGRAL'].append(float('NaN'))

                if len(data['MASTER NCONSS']) < it:
                    data['MASTER NCONSS'].append(-1)
                if len(data['MASTER NVARS']) < it:
                    data['MASTER NVARS'].append(-1)
                if len(data['LINKING VARS']) < it:
                    data['LINKING VARS'].append(float('NaN'))

                if len(data['BNB TREE NODES']) < it:
                    data['BNB TREE NODES'].append(-1)
                if len(data['BNB TREE LEFT']) < it:
                    data['BNB TREE LEFT'].append(-1)
                if len(data['BNB TREE DEPTH']) < it:
                    data['BNB TREE DEPTH'].append(-1)

                # times (exectime and setuptime)
                if len(data['BR RULE TIME GENERIC']) < it:
                    data['BR RULE TIME GENERIC'].append(-1)
                if len(data['BR RULE TIME ORIG']) < it:
                    data['BR RULE TIME ORIG'].append(-1)
                if len(data['BR RULE TIME RELPSPROB']) < it:
                    data['BR RULE TIME RELPSPROB'].append(-1)
                if len(data['BR RULE TIME RYANFOSTER']) < it:
                    data['BR RULE TIME RYANFOSTER'].append(-1)

                # calls (lp, ext and ps)
                if len(data['BR RULE CALLS GENERIC']) < it:
                    data['BR RULE CALLS GENERIC'].append(-1)
                if len(data['BR RULE CALLS ORIG']) < it:
                    data['BR RULE CALLS ORIG'].append(-1)
                if len(data['BR RULE CALLS RELPSPROB']) < it:
                    data['BR RULE CALLS RELPSPROB'].append(-1)
                if len(data['BR RULE CALLS RYANFOSTER']) < it:
                    data['BR RULE CALLS RYANFOSTER'].append(-1)

                if len(data['RMP LP CALLS']) < it:
                    data['RMP LP CALLS'].append(-1)
                if len(data['RMP LP TIME']) < it:
                    data['RMP LP TIME'].append(0.)
                if len(data['RMP LP ITERATIONS']) < it:
                    data['RMP LP ITERATIONS'].append(-1)

                if len(data['ORIGINAL LP CALLS']) < it:
                    data['ORIGINAL LP CALLS'].append(-1)
                if len(data['ORIGINAL LP TIME']) < it:
                    data['ORIGINAL LP TIME'].append(float('NaN'))
                if len(data['ORIGINAL LP ITERATIONS']) < it:
                    data['ORIGINAL LP ITERATIONS'].append(-1)

                if len(data['LP FILE']) < it:
                    data['LP FILE'].append(-1)
                if len(data['DEC FILE']) < it:
                    data['DEC FILE'].append(-1)

        datalengths = []
        for key in data:
            datalengths.append(len(data[key]))
            d[key] += data[key]
            data[key] = []

        if len(set(datalengths)) == 2:
            print("One error in input. Possibly unrecoverable.")
            differentlen = [count(x) for x in set(datalengths)].min()
            differentlenidx = datalenghts.index(differentlen)
            for l, key in enumerate(data):
                print(outfile, data[data.keys[differentlenidx]], datalengths(differentlenidx))
        elif len(set(datalengths)) > 2:
            print("Multiple errors in input. Unrecoverable.")
            for l, key in enumerate(data):
                print(outfile, key, datalengths[l])


        idx += index
        index = []
        it = 0


    # build pandas data frame
    #pd.set_option("max_columns", 999)
    try: df = pd.DataFrame(index=idx, data=d)
    except:
        for key in d.keys():
            if len(d[key]) != len(idx):
                print(f"Fatal: Not as many entries as expected for key '{key}' (expected: {len(idx)}, got: {len(d[key])}).\n"+\
                       "       This could be due to a keyword appearing multiple times (e.g. once in original and once in master problem).\n"+\
                       f"       Instances: {idx}\n"+\
                       f"       {key}: {d[key]}\n"+\
                       "       Terminating.")
                exit()

    df.drop_duplicates(subset=["LP FILE", "DEC FILE"], keep="last", inplace=True)
    return df

def main(outfiles, save=True, path=""):
    # check command line arguments
    if len(outfiles) < 1:
        sys.exit("Usage: ./parser_general.py OUTFILE1 [OUTFILE2 ...]")
    df = parseOutfiles(outfiles)

    if save:
        print("Saving to", os.path.join(path,'{}.general.pkl'.format(outfiles[0].split('/')[-1])))
        df.to_pickle(os.path.join(path,'{}.general.pkl'.format(outfiles[0].split('/')[-1])))
    return df

if __name__ == '__main__':
    args = sys.argv[1:]
    main(args)
