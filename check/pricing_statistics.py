#!/usr/bin/env python2

import sys
import os
import argparse
import time

import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import lines
from matplotlib import transforms
from matplotlib import ticker

# Define global variables
params = {}

def parse_arguments(args):
    """
    Parse the command-line arguments
    :param args: Command-line arguments provided during execution
    :return: Parsed argumetns
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--outdir', type=str,
                        default="plots",
                        help='output directory (default: "plots")')

    parser.add_argument('-a', '--allnodes', action='store_true',
                        help='set this flag, to collect data for all nodes, not just the root')

    parser.add_argument('-p', '--png', action='store_true',
                        help='set this flag for a non-zoomable png-plot as output')

    parser.add_argument('-s', '--splitrounds', type=int,
                        default=0,
                        help='the plot can be split into pieces containing a maximum of SPLITROUNDS rounds each (default is no splitting)')

    parser.add_argument('-m', '--minround', type=int,
                        default=1,
                        help='start the data-collection with pricing-round MINROUND (default is the first round)')

    parser.add_argument('-M', '--maxround', type=int,
                        default=0,
                        help='end the data-collection with pricing-round MAXROUND (default is MAXROUND=0, which will be the last round)')

    parser.add_argument('-c', '--colors', type=str,
                        default="nipy_spectral",
                        help='name of the color-map, that is used for the bars (see matplotlib documentation for maps, default is nipy_spectral)')

    parser.add_argument('-l', '--lines', action='store_true',
                        help='enforce lines between pricing-rounds on the plots (default is to not draw lines for rounds, that are too short)')

    parser.add_argument('-i', '--instances', type=str, nargs = '*',
                        default="",
                        help='names of the instances to be included in the plot/data-collection (default is all instances in FILENAMES)')

    parser.add_argument('-z', '--summary-only', action='store_true',
                        help='create only summary plots')

    parser.add_argument('-Z', '--no-summary', action='store_true',
                        help='create no summary plots')

    parser.add_argument('-t', '--no-text', action='store_true',
                        help='do not write any text on the plots (such as node or round numbers)')

    parser.add_argument('-S', '--save', action='store_true',
                        help='saves the collected data in a pickle-file, in the OUTDIR, instead of plotting it (see --load)')

    parser.add_argument('-L', '--load', action='store_true',
                        help='loads earlier collected data from a pickle-file instead of parsing a GCG-outfile (see --save)')

    parser.add_argument('filenames', nargs='+',
                        help='Names of the files to be used for the bound plots')

    parsed_args = parser.parse_args(args)

    # check if the provided arguments are consistent
    if (0 < parsed_args.maxround and parsed_args.maxround < parsed_args.minround) or parsed_args.minround <= 0 or parsed_args.maxround < 0:
        print 'please make sure that 1 <= MINROUND <= MAXROUND (or MAXROUND == 0 for the last possible round)'
        exit()
    if not parsed_args.colors in plt.cm.datad:
        print 'please use a colormap that is supported by pyplot (' + parsed_args.colors + ' is not supported)'
        exit()

    return parsed_args

def set_params(args):
    """
    Set the global parameters from the parsed command-line arguments
    :param args: parsed command-line arguments
    :return:
    """
    params['outdir'] = args.outdir
    params['root_only'] = not args.allnodes
    params['details'] = not args.png
    params['nSplit'] = args.splitrounds
    params['minRound'] = args.minround
    params['maxRound'] = args.maxround
    params['colors'] = args.colors
    params['lines'] = args.lines
    params['instances'] = args.instances
    params['summary_only'] = args.summary_only
    params['no_summary'] = args.no_summary
    params['no_text'] = args.no_text
    params['save'] = args.save
    params['load'] = args.load

def get_colmap(pricers):
    """
    Returns a list of colors, with same length as pricers, that can be used for the bar plot
    Each pricing problem has its own color
    :param pricers: a list with pricing_problem ids
    :return: a list of colors as demanded by pyplot.bar()
    """
    pricer_to_color = {}
    col_ind = 0

    # build the mapping pricer id -> color id
    for p in pricers:
        if not p in pricer_to_color:
            pricer_to_color[p] = col_ind
            col_ind += 1

    # get a color map of the right length, so that each color-id gets its own color
    cmap = plt.get_cmap(params['colors'],len(pricer_to_color))

    # build a list of colors and return it
    colors = [cmap(pricer_to_color[p]) for p in pricers]
    return colors

def remove_overlapping_texts(figure, texts):
    """
    Removes all texts in figure from the list, that overlap, by setting their visibility to False
    :param figure: the figure to which the texts belong
    :param texts: list of texts, that are to be checked
    :return:
    """
    rend = figure.canvas.get_renderer()
    for i,txt in enumerate(texts):
        if not txt.get_visible():
            continue
        bb = txt.get_window_extent(renderer=rend)
        for nxt_txt in texts[(i+1):]:
            nxt_bb = nxt_txt.get_window_extent(renderer=rend)
            if bb.overlaps(nxt_bb):
                nxt_txt.set_visible(False)
            else:
                break

def make_plot(data, name):
    """
    Make the plots from the structured data
    :param data: dataframe with the collected data
    :param name: the problemname
    :return:
    """
    start_time = time.time()

    # set the heights of zero
    ymin = -0.15

    # flat out the data again (maybe in the final version, there should be no dataframe to begin with)
    flat_data = data.reset_index()
    flat_data.time = flat_data.time + 0.01 # workaround for most times being zero (0,01s is SCIPs smallest time-interval)

    # define position, width and height of the peaks; the former are defined by time, the latter by nVars
    x = (flat_data.time.cumsum()-flat_data.time).values
    y = (flat_data.nVars - ymin).values
    widths = flat_data.time.values
    colors = get_colmap(flat_data['pricing_prob'].values)

    print '    data restructured:', time.time() - start_time
    start_time = time.time()

    # make the bar plot
    if params['details']:
        lw = 0.01
    else:
        lw = 1.0
    plt.bar(x, y, widths, bottom = ymin, align = 'edge', linewidth = lw, edgecolor = 'k', color = colors, label='pricing problems')
    fig = plt.gcf()
    ax = plt.gca()

    print '    data plotted:', time.time() - start_time
    start_time = time.time()

    # set parameters
    if params['details']:
        fig.set_size_inches(8*11.7,8*8.3)
        textsize = ax.get_window_extent().height * 0.013
    else:
        fig.set_size_inches(11.7,8.3)
        textsize = 12
    totalTime = max(x)
    ymax = max(y)
    ax.set_ylim([ymin,ymax])

    # formatting
    ax.set_xlim([0,x[-1]+widths[-1]])
    old_yticks = ax.get_yticks()
    new_yticks = []
    for i,n in enumerate(old_yticks):
        if abs(n - int(n)) < 0.01 and n >= 0:
            new_yticks.append(n)
    ax.set_yticks(new_yticks)
    ax.tick_params(axis = 'both', length = textsize/2, width = textsize/40, labelsize = textsize*0.9)
    ax.set_xlabel('Time / s', size = 1.15*textsize)
    ax.set_ylabel('\# of variables', size = 1.15*textsize)
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

    print '    data formatted:', time.time() - start_time
    start_time = time.time()

    if params['no_text']:
        # save the figure
        if params['details']:
            plt.savefig(params['outdir'] + '/' + name + '.pdf')
        else:
            plt.savefig(params['outdir'] + '/' + name + '.png')
        plt.close()

        print '    save:', time.time() - start_time
        print '    saved figure'
    else:
        # add information about the stabilization & pricing rounds
        prev_rnd = flat_data['pricing_round'][0]
        prev_stab = flat_data['stab_round'][0]
        prev_x = 0
        texts = []
        enfLine = False
        for i in range(len(x)):
            rnd = flat_data['pricing_round'][i]
            stab = flat_data['stab_round'][i]
            if stab > prev_stab or rnd > prev_rnd:
                if rnd > prev_rnd:
                    # bold line for a new pricing round
                    if params['lines'] or (x[i] - prev_x)/totalTime > 0.005 or enfLine:
                        line = lines.Line2D([x[i],x[i]],[0,1],color='r',linewidth=1.0, transform = trans)
                        ax.add_line(line)
                        if (x[i] - prev_x)/totalTime > 0.025:
                            enfLine = True
                        else:
                            enfLine = False
                    texts.append(ax.text(prev_x, 1.01, 'Round '+str(prev_rnd), rotation='vertical',va='bottom', ha='left', size = textsize, transform = trans))
                    prev_rnd = rnd
                    prev_stab = 1
                    prev_x = x[i]
                else:
                    # dashed line for a new stabilization round
                    line = lines.Line2D([x[i],x[i]],[0,1],color='orange',linestyle='--',linewidth=0.8, transform = trans)
                    ax.add_line(line)
                    prev_stab = stab
        texts.append(ax.text(prev_x, 1.01, 'Round '+str(prev_rnd), rotation='vertical',va='bottom', ha='left', size = textsize, transform = trans))
        text_height = [t for t in texts if t.get_visible()][-1].get_window_extent(renderer = fig.canvas.get_renderer()).transformed(ax.transAxes.inverted()).y1

        # check for overlapping texts
        remove_overlapping_texts(fig,texts)

        print '    stab- and pricing-round information:', time.time() - start_time
        start_time = time.time()

        # add information about the nodes
        prev_node = flat_data['node'][0]
        prev_x = 0
        texts = []
        for i in range(len(x)):
            node = flat_data['node'][i]
            if node > prev_node:
                line = lines.Line2D([x[i],x[i]],[1,text_height+0.01],color='r',linewidth=1.0, transform = trans)
                line.set_clip_on(False)
                ax.add_line(line)
                texts.append(ax.text(prev_x, text_height+0.02, 'Node '+str(prev_node), ha='left', size = textsize, style='italic', transform = trans))
                prev_node = node
                prev_x = x[i]
        texts.append(ax.text(prev_x, text_height+0.02, 'Node '+str(prev_node), ha='left', size = textsize, style='italic', transform = trans))

        # check for overlapping texts
        remove_overlapping_texts(fig,texts)

        text_height = [t for t in texts if t.get_visible()][-1].get_window_extent(renderer = fig.canvas.get_renderer()).transformed(ax.transAxes.inverted()).y1

        print '    node information:', time.time() - start_time
        start_time = time.time()

        # save the figure
        fig.subplots_adjust(top=0.98/text_height)
        if params['details']:
            plt.savefig(params['outdir'] + '/' + name + '.pdf')
        else:
            plt.savefig(params['outdir'] + '/' + name + '.png')
        plt.close()

        print '    save:', time.time() - start_time
        print '    saved figure'

def make_summary_plot(data, name):
    """
    For each problem create one summary plot, which shows for each pricing round the
    cumulated running time of the pricers in this round as well as the fraction of
    pricers, that did not find a variable.
    :param data: dataframe with the collected data
    :param name: the problemname
    :return:
    """
    start_time = time.time()

    summary = pd.DataFrame()
    summary['time'] = data.groupby(level=['node','pricing_round','stab_round']).sum().time
    summary['found_frac'] = data.astype(bool).groupby(level=['node','pricing_round','stab_round']).sum().nVars/data.groupby(level=['node','pricing_round','stab_round']).count().nVars*100
    summary = summary.reset_index()

    print '    extracted summary data:', time.time() - start_time
    start_time = time.time()

    fig,ax1 = plt.subplots()
    ax2 = ax1.twinx()

    # get the data for the plot
    x = [i+1 for i in summary.index.values]
    y_time = summary.time.values
    y_found_frac = summary.found_frac.values

    # format the plot
    ax1.set_xlabel('Pricing Round', size='large')
    ax1.set_ylabel('Time / s', color='k', size='large')
    ax2.set_ylabel('Fraction of successfull pricers / \%', color='r')

    ax1.get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    ax1.set_xlim([0, max(x) + 0.9])
    if max(y_time) > 0:
        ax1.set_ylim([-max(y_time)*0.1,max(y_time)*1.1])
    else:
        ax1.set_ylim([-0.001,0.01])
    ax2.set_ylim([-max(y_found_frac)*0.15,max(y_found_frac)*1.15])

    p1,p2 = ax1.transData.transform([(1,0),(2,0)])
    perimeter = max([3.5,(p2[0] - p1[0])/5.])

    # plot the data
    ax1.scatter(x,y_time, color='k', s=perimeter**2)
    ax2.scatter(x,y_found_frac, color='r', s=perimeter**2)

    # add a line after the root-node
    if summary.node.max() > 1:
        x_line = (summary[summary.node > 1].index[0] + 1 + summary[summary.node == 1].index[-1] + 1)/2.
        line = lines.Line2D([x_line,x_line],[0,1],color='orange',linewidth=.5,linestyle='--', transform = transforms.blended_transform_factory(ax1.transData, ax1.transAxes))
        ax1.add_line(line)
        ax2.text(x_line, 50, "\it{End of Root}", va = 'center', ha = 'left', rotation = 90, size = 'smaller', color = 'orange', zorder = 1)
    print '    plotted summary:', time.time() - start_time
    start_time = time.time()

    # save the plot
    fig.set_size_inches(11.7,8.3)
    plt.tight_layout()
    if params['details']:
        plt.savefig(params['outdir'] + '/' + name + '_summary.pdf')
    else:
        plt.savefig(params['outdir'] + '/' + name + '_summary.png')
    plt.close()
    print '    saved summary:', time.time() - start_time

def make_bubble_plot(data, name):
    data = data.reset_index()
    # add a column 'round' to the dataframe, that counts each stab_round and pricing_round
    rnd = []
    prev_rnd = 1
    for ind in data.index:
        if ind == 0:
            rnd.append(1)
            continue
        if data.loc[ind].pricing_round > data.loc[ind-1].pricing_round or data.loc[ind].stab_round > data.loc[ind-1].stab_round:
            rnd.append(prev_rnd + 1)
            prev_rnd = prev_rnd + 1
        else:
            rnd.append(prev_rnd)
    data['round'] = rnd

    pricers = data.pricing_prob.unique().tolist()

    # add x and y data to plot for every pricer
    x = {}
    y = {}
    for p in pricers:
        x[p] = data[(data.pricing_prob == p) & (data.nVars >= 1)]['round'].values
        y[p] = [p for i in x[p]]

    colors = get_colmap(pricers)

    fig = plt.gcf()
    ax = plt.gca()

    ax.set_xlabel('Round')
    ax.set_ylabel('Pricer ID')
    ax.set_xlim([0,data['round'].max()+0.9])
    ax.set_ylim([min(pricers) - 0.5, max(pricers) + 0.5])
    ax.get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    ax.get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    p1,p2 = ax.transData.transform([(1,min(pricers)),(2,min(pricers))])
    perimeter = max([3.5,(p2[0] - p1[0])/2.2])

    for p in pricers:
        ax.scatter(x[p],y[p], color = colors[pricers.index(p)], s=perimeter**2)

    # add a line after the root-node
    if data.node.max() > 1:
        x_line = (data[data.node > 1]['round'].iloc[0] + data[data.node == 1]['round'].iloc[-1])/2.
        line = lines.Line2D([x_line,x_line],[0,1],color='orange',linewidth=.5,linestyle='--', transform = transforms.blended_transform_factory(ax.transData, ax.transAxes))
        ax.add_line(line)
        ax.text(x_line, (max(pricers) + min(pricers))/2., "\it{End of Root}", va = 'center', ha = 'left', rotation = 90, size = 'smaller', color = 'orange', zorder = 1)

    # save the plot
    fig.set_size_inches(11.7,8.3)
    plt.tight_layout()
    plt.savefig(params['outdir'] + '/' + name + '_bubble.pdf')
    plt.close()

def plots(data, name):
    """
    Master-function for plotting. Splits the data, if necessary and calls make_plot()
    :param data: collected data as dataframe
    :param name: name of the problem
    :return:
    """
    # use tex to render the text output
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    make_bubble_plot(data, name)

    if not params['no_summary']:
        # build the summary plot
        make_summary_plot(data, name)
        if params['summary_only']:
            return

    if params['maxRound'] <= 0:
        maxRnd = max(data.index.get_level_values('pricing_round').values)
    else:
        maxRnd = min(params['maxRound'],max(data.index.get_level_values('pricing_round').values))
    if params['minRound'] <= 1:
        minRnd  = min(data.index.get_level_values('pricing_round').values)
    else:
        minRnd  = max(params['minRound'],min(data.index.get_level_values('pricing_round').values))
    if params['root_only']:
        maxNode = 1
    else:
        maxNode = max(data.index.get_level_values('node').values)

    if params['nSplit'] <= 0:
        # do not split the plot, but still check if rounds were neglected
        if params['minRound'] > 1 or params['maxRound'] > 0:
            name += '_rounds' + str(params['minRound']) + 'to' + str(maxRnd)
        make_plot(data.query('node <= @maxNode'), name)
    else:
        # split the plot by rounds
        fromRnd = minRnd - 1
        for i in range(1,(maxRnd-minRnd)+1):
            if i % params['nSplit'] <> 0 and i <> (maxRnd-minRnd):
                continue
            toRnd = i+minRnd
            make_plot(data.query('@fromRnd < pricing_round <= @toRnd & node <= @maxNode'), name + '_rounds' + str(fromRnd + 1) + 'to' + str(toRnd))
            fromRnd = toRnd

    # this line does not do anything, but removing code-analysis warnings, as the variable is not used where the parser could see it..
    maxNode += 1

def load_data(files):
    """
    Plots data, that was parsed and collected earlier from the generate_files() method and saved in a pickle-file
    :param files: the pickle files, from which the dataframes are to load (one file per instance)
    :return:
    """
    for file in files:
        if not os.path.exists(file):
            print 'there is no file ' + file
            continue
        name = os.path.splitext(os.path.basename(file))[0]
        if params['instances'] <> '' and not (name in params['instances']):
            print 'skipping', name
            continue
        print 'entering', name
        start_time = time.time()
        df = pd.read_pickle(file)
        print '    loading data:', time.time() - start_time
        if df.empty:
            print '    no data found'
        else:
            start_time = time.time()
            plots(df, name)
            print '    total plotting:', time.time() - start_time
        print '    leaving', name

def collect_data(name, ind_node, ind_pricing_round, ind_stab_round, ind_pricing_prob, val_time, val_nVars):
    index = pd.MultiIndex.from_arrays([ind_node, ind_pricing_round, ind_stab_round, ind_pricing_prob],
                                      names=["node", "pricing_round", "stab_round", "pricing_prob"])
    data = {'time': val_time, 'nVars': val_nVars}
    df = pd.DataFrame(data=data, index = index)

    if params['save']:
        start_time = time.time()
        df.to_pickle(params['outdir'] + '/' + name + '.pkl')
        print '    total saving:', time.time() - start_time

    else:
        start_time = time.time()
        plots(df, name)
        print '    total plotting:', time.time() - start_time

def generate_files(files):
    """
    Parse the files and structure the pricing-data in a dataframe then make the plots
    :param files: List of files to be parsed
    :return:
    """
    for file in files:
        with open(file) as _file:
            # initialize all index-lists
            ind_node = []
            ind_pricing_round = []
            ind_stab_round = []
            ind_pricing_prob = []

            # initialize the value-lists
            val_time = []
            val_nVars = []

            # initialize all counters
            node = 0
            pricing_round = 0
            stab_round = 0
            pricing_prob = 0

            # initialize all other variables
            problemFileName = None
            done = False

            for line in _file:
                if line.find("@0")<>-1:
                    # if the file is a out-file, generated by the check-script, reset the variables whenever a new instance starts
                    if line.startswith("@01"):
                        # print message, if the previous problem is not done yet
                        if not done and problemFileName:
                            print '    ended abruptly'
                            collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_pricing_prob, val_time, val_nVars)
                            print '    leaving', problemFileName
                            done = True

                        # initialize all index-lists
                        ind_node = []
                        ind_pricing_round = []
                        ind_stab_round = []
                        ind_pricing_prob = []

                        # initialize the value-lists
                        val_time = []
                        val_nVars = []

                        # initialize all counters
                        node = 0
                        pricing_round = 0
                        stab_round = 0
                        pricing_prob = 0

                        # initialize all other variables
                        problemFileName = None
                        done = False
                    else:
                        # ignore lines, where the output ends abrubtly (e.g. when the hard limit of the check-script is reached)
                        continue

                # if the problem is already processed, continue
                elif done:
                    continue

                elif not problemFileName and line.startswith("read problem "):
                    # get the problem name from the file name as in "check.awk"
                    tmparray = line.split("<")[-1].replace(">","").replace("\n","").split("/")[-1].split(".")
                    problemFileName = tmparray[0]
                    if tmparray[-1] == "gz" or tmparray[-1] == "z" or tmparray[-1] == "GZ" or tmparray[-1] == "Z":
                        tmparray.pop()
                    for i in range(1,len(tmparray)-1):
                        problemFileName += "." + tmparray[i]
                    if params['instances'] <> '' and not (problemFileName in params['instances']):
                        print 'skipping', problemFileName
                        done = True
                        continue
                    print 'entering', problemFileName

                # pricer statistics end
                elif line.startswith("SCIP Status        :"):
                    # continue if no data is found
                    if not ind_node or not ind_pricing_round or not ind_stab_round or not ind_pricing_prob or not val_time or not val_nVars:
                        print '    no pricing data found'
                        done = True
                        continue

                    collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_pricing_prob, val_time, val_nVars)

                    done = True

                    print '    leaving', problemFileName
                    continue

                # ignore all other lines, that do not contain pricer statistics messages
                elif not line.startswith("[src/pricer_gcg.cpp:"):
                    continue

                # extract the pricing-statistics message
                message = line.split("] statistic: ")[-1]

                if message.startswith("New pricing round at node"):
                    try:
                        node = int(message.split()[-1])
                        pricing_round += 1
                    except ValueError:
                        print '    ended abruptly'
                        collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_pricing_prob, val_time, val_nVars)
                        print '    leaving', problemFileName
                        done = True
                        continue

                elif message.startswith("Stabilization round "):
                    try:
                        stab_round = int(message.split()[-1])
                    except ValueError:
                        print '    ended abruptly'
                        collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_pricing_prob, val_time, val_nVars)
                        print '    leaving', problemFileName
                        done = True
                        continue

                elif message.startswith("Pricing prob "):
                    try:
                        pricing_prob = int(message.split()[2])

                        # check if the pricing prob should be included in the plot
                        if (params['root_only'] and node > 1) or pricing_round < params['minRound'] or (0 < params['maxRound'] < pricing_round):
                            continue

                        # store all indices
                        ind_node.append(node)
                        ind_pricing_round.append(pricing_round)
                        ind_stab_round.append(stab_round)
                        ind_pricing_prob.append(pricing_prob)

                        # store the data
                        val_time.append(float(message.split()[-1]))
                        val_nVars.append(int(message.split()[5]))
                    except ValueError:
                        print '    ended abruptly'
                        collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_pricing_prob, val_time, val_nVars)
                        print '    leaving', problemFileName
                        done = True
                        continue

def main():
    """Entry point when calling this script"""
    args = sys.argv[1:]
    parsed_args = parse_arguments(args)
    set_params(parsed_args)
    if not os.path.exists(params['outdir']):
        os.makedirs(params['outdir'])
    if params['load']:
        load_data(parsed_args.filenames)
    else:
        generate_files(parsed_args.filenames)

# Calling main script
if __name__ == '__main__':
    main()
