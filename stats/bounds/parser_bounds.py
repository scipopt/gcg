#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import pickle as pickler
#import shutil

def generate_files(files,params):
    """
    Parse the files and generate temporary files containing only problem name and execution time
    :param files: List of files to be parsed
    :return: A list of all the generated files to be deleted after performance profiling
    """

    # Create a dictionary, where all the dataframes, that are generated in the following, are 'globally' stored for comparison
    df_dict = {}
    set_dict = {}

    for file in files:
        # file = os.path.join(DIR, filename)
        with open(file) as _file:
            df = None
            dfvar = None
            orig = False
            name = None
            problemFileName = None
            rootbounds = False
            vardetails = False
            settings = 'default'
            varlines = {}
            varheader = None
            boundlines = {}
            boundheader = None
            scip_status = ""
            for line in _file:
                if line.startswith("@04"):
                    # store the current dataframe globally
                    if params['compare'] and df is not None:
                        if not (name in df_dict):
                            df_dict[name] = []
                            set_dict[name] = []
                        df_dict[name].append(df.copy())
                        set_dict[name].append(settings)
                    # Save dataframe for finished instance - maybe a feature for the future, if one wants to have pickles of single instances
                    #if params['save']:
                    #    if df is not None:
                    #        # print("Exporting dataframe for instance {}".format(name))
                    #        #df.to_pickle("./{}/{}.bounds.pkl".format(params['dataframedir'],name))
                    #        df.to_pickle("./{}/{}.bounds.pkl".format(params['dataframedir'],name))
                    #    else:
                    #        print("Dataframe for instance {} is empty. Have you tested with STATISTICS=true?".format(name))

                if line.startswith("@01"):
                    # reset python variables for next instance
                    problemFileName = None
                    df = None
                    dfvar = None
                    boundheader = None
                    varheader = None
                    varlines = {}
                    boundlines = {}
                if line.startswith("GCG> set load"):
                    # store current settings
                    settings=line.split()[-1]
                    settings=settings.split("/")[-1]
                    settings = os.path.splitext(settings)[0]
                elif not problemFileName and line.startswith("read problem "):
                    # get the problem name from the file name as in "check.awk", in case it is "BLANK" in the actual "Problem name"-line
                    tmparray = line.split("<")[-1].replace(">","").replace("\n","").split("/")[-1].split(".")
                    problemFileName = tmparray[0]
                    if tmparray[-1] == "gz" or tmparray[-1] == "z" or tmparray[-1] == "GZ" or tmparray[-1] == "Z":
                        tmparray.pop()
                    for i in range(1,len(tmparray)-1):
                        problemFileName += "." + tmparray[i]
                elif not orig and line.startswith("Original Program statistics:"):
                    orig = True
                elif orig and line.startswith("Master Program statistics:"):
                    orig = False
                elif orig and line.startswith("Presolved Problem  :"):
                    orig = False
                elif orig and line.startswith("SCIP Status        :"):
                    scip_status = line.split(":")[-1].strip()
                elif orig and line.startswith("  Problem name     :"):
                    # store problem name
                    name = line.split()[3]
                    name = name.split("/")[-1]
                    tmp_name = name.split(".")[-1]
                    if tmp_name[-1] == "gz" or tmp_name[-1] == "z" or tmp_name[-1] == "GZ" or tmp_name[-1] == "Z":
                        name = os.path.splitext(name)[0]
                    name = os.path.splitext(name)[0]
                    if name == 'BLANK':
                        name = problemFileName
                elif not rootbounds and line.startswith("Root bounds"):
                    # prepare storage of root bounds
                    rootbounds = True
                elif rootbounds and line.startswith(" iter            pb            db"):
                    # store root bounds header
                    line_array = line.split()
                    # add gap column
                    line_array.append('gap')
                    boundheader = line_array
                elif rootbounds and line.startswith("Pricing Summary:"):
                    # finished with storing root bounds
                    rootbounds = False
                    # correct "local" gap to "global" gap (current best bound)
                    if True: #not params['allgaps']:
                        for currKey in boundlines:
                            for prevKeys in boundlines:
                                if int(prevKeys) < int(currKey) and float(boundlines[prevKeys][-1]) < float(boundlines[currKey][-1]):
                                    boundlines[currKey][-1] = boundlines[prevKeys][-1]
                    #exit()
                elif rootbounds:
                    # store root bound line
                    line_array = line.split()
                    # Calculate gap
                    z_lower = float(sorted([line_array[1],line_array[2]],reverse=True)[0])
                    z_upper = float(sorted([line_array[1],line_array[2]],reverse=True)[1])
                    if z_upper*z_lower > 0:
                        gapvar = np.absolute( (z_upper - z_lower) / np.absolute(z_lower))
                        #if name == "I06": print("({}-{})/{} = {}".format(z_upper , z_lower,z_lower,gapvar))
                    else:
                        gapvar = np.inf
                    # print("Appending gapvar {} for instance {} at time {}".format(gapvar,name,line_array[3]))
                    line_array.append(gapvar)
                    boundlines[line_array[0]] = line_array
                elif not vardetails and line.startswith("AddedVarDetails:"):
                    # prepare storage of var details
                    vardetails = True
                elif vardetails and line.startswith("VAR: name	node	time") and vardetails:
                    # store var details header
                    line_array = line.split()
                    varheader = line_array[1:]
                elif vardetails and line.startswith("VAR:") and not int(line.split()[2]) == 1:
                    # ignore variables that were not create in the root node
                    continue
                elif vardetails and line.startswith("Root node:"):
                    # finished reading var details (and root bounds)
                    vardetails = False

                    # create dict with root bounds header
                    boundmap = {}
                    try:
                        for i in range(len(boundheader)):
                            boundmap[i] = boundheader[i]
                    except TypeError:
                        print("   Warning: Root bounds could not be found.\n            Perhaps you did not compile GCG *and* SCIP with STATISTICS=true.\n   Terminating.")
                        exit()

                    # use boundlines dict to create data frame
                    df = pd.DataFrame.from_dict(data = boundlines, orient = 'index', dtype = float)

                    # if no root bounds are present, ignore instance
                    if len(df) == 0:
                        print(name)
                        print("   -> ignored")
                        print("   -> SCIP Status : {}".format(scip_status))
                        continue

                    # use root bounds header to rename columns of data frame
                    df.rename(columns = boundmap, inplace=True)

                    # sort lines according to iteration
                    df.sort_values(by='iter', inplace=True)

                    # create var data frame from varlines dict
                    dfvar = pd.DataFrame.from_dict(data = varlines, orient = 'index', dtype = float)

                    # create dict with var header
                    varmap = {}
                    for i in range(len(varheader)):
                        varmap[i] = varheader[i]

                    # use var header to rename columns of var data frame
                    dfvar.rename(columns = varmap, inplace=True)

                    # set index of var data frame to name of var
                    dfvar = dfvar.set_index(keys='name')

                    # set type of var data frame
                    dfvar=dfvar.astype(float)

                    # create new column in data frame containing the number of lp vars generated in each iteration
                    df['nlpvars'] = 0
                    for i in range(len(df)):
                        df.at[str(i), 'nlpvars'] = len(dfvar[(dfvar['rootlpsolval'] != 0) & (dfvar['rootredcostcall'] == i)])

                    # add the number of all lp-variables, not created by reduced cost pricing (e.g. by Farkas-Pricing)
                    if params['farkas']:
                        df.at[str(0),'nlpvars'] = df['nlpvars'][0] + len(dfvar[(dfvar['rootlpsolval'] != 0) & (dfvar['rootredcostcall'] == -1.)])

                    # create new column in data frame containing the number of lp vars generated until each iteration
                    df['nlpvars_cum'] = df[(df['iter'] < len(df))].cumsum(axis=0)['nlpvars']

                    # compute total number of vars in root lp solution, that are included in the plot
                    nlpvars_total = df['nlpvars_cum'].iloc[-1]

                    # create new column in data frame containing the percentage of lp vars generated until each iteration
                    df['lpvars'] = df['nlpvars_cum']/nlpvars_total

                    # repeat this for the vars (generated at the root) in ip solution
                    df['nipvars'] = 0

                    for i in range(len(df)):
                        df.at[str(i), 'nipvars'] = len(dfvar[(dfvar['solval'] > 0) & (dfvar['rootredcostcall'] == i)])

                    if params['farkas']:
                        df.at[str(0),'nipvars'] = df['nipvars'][0] + len(dfvar[(dfvar['solval'] > 0) & (dfvar['rootredcostcall'] == -1.)])

                    df['nipvars_cum'] = df[(df['iter'] < len(df))].cumsum(axis=0)['nipvars']

                    nipvars_total = df['nipvars_cum'].iloc[-1]

                    df['ipvars'] = df['nipvars_cum']/nipvars_total

                    # set type of data frame
                    df=df.astype(float)

                    # set infty
                    infty = 10.0 ** 20

                    # set dual bounds of -infinity to NAN
                    df['db'][df['db'] <= -infty] = np.nan
                    #df=df.dropna()

                    # workaround for iterations that were done at the same time (SCIP only counts in 1/100 of a second)
                    df['time_count'] = df.groupby('time')['time'].transform('count')
                    df['time_first'] = df.groupby('time')['iter'].transform('first')
                    df['time'] = df['time'] + 0.01*(df['iter'] - df['time_first'])/df['time_count']
                    df['time_diff'] = df["time"].diff(1)
                    # Chained assignment incoming
                    pd.options.mode.chained_assignment = None
                    df['time_diff'][0] = df['time'][0]

                    df['db_ma'] = df['db'].rolling(window=5,center=False).mean()
                elif vardetails:
                    # store details of variable
                    line_array = line.split()
                    varlines[line_array[1]] = line_array[1:]
                    ## End of pickling

    if params['save']:
        with open(os.path.join(params['outdir'],"{}.boundsdict.pkl".format(file.split('/')[-1])), 'wb') as handle:
            #print("Dumping Boundsdict to {}".format(str(handle)))
            pickler.dump(df_dict, handle, protocol=pickler.HIGHEST_PROTOCOL)
        with open(os.path.join(params['outdir'],"{}.boundsset.pkl".format(file.split('/')[-1])), 'wb') as handle:
            #print("Dumping Boundsset to {}".format(str(handle).split('=')[-1]))
            pickler.dump(set_dict, handle, protocol=pickler.HIGHEST_PROTOCOL)
        exit()
    else: return df_dict, set_dict
