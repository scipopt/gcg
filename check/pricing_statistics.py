#!/usr/bin/env python2

import sys
import os
import argparse
import time
import datetime

import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import lines
from matplotlib import transforms
from matplotlib import ticker

# Define the global parameter list
params = {}

def parse_arguments(args):
    """
    Parse the command-line arguments
    :param args: Command-line arguments provided during execution
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--outdir', type=str,
                        default="plots",
                        help='output directory (default: "plots")')

    parser.add_argument('-r', '--root-only', action='store_true',
                        help='set this flag, to collect data just for the root-node')

    parser.add_argument('-p', '--png', action='store_true',
                        help='set this flag for a non-zoomable png-plot as output')

    parser.add_argument('-s', '--splitrounds', type=int,
                        default=0,
                        help='the complete plot can be split into pieces containing a maximum of SPLITROUNDS rounds each (default is no splitting)')

    parser.add_argument('-m', '--minround', type=int,
                        default=1,
                        help='start the data-collection or the complete plot with pricing-round MINROUND (default is the first round)')

    parser.add_argument('-M', '--maxround', type=int,
                        default=0,
                        help='end the data-collection or the complete plot with pricing-round MAXROUND (default is MAXROUND=0, which will be the last round)')

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

    parser.add_argument('-b', '--bubble-only', action='store_true',
                        help='create only bubble plots')

    parser.add_argument('-B', '--no-bubble', action='store_true',
                        help='create no bubble plots')

    parser.add_argument('-t', '--no-text', action='store_true',
                        help='do not write any text on the plots (such as node or round numbers)')

    parser.add_argument('-S', '--save', action='store_true',
                        help='saves the collected data in a pickle-file, in the OUTDIR, instead of plotting it (see also --load)')

    parser.add_argument('-L', '--load', action='store_true',
                        help='loads earlier collected data from a pickle-file instead of parsing a GCG-outfile and plots it (see also --save)')

    parser.add_argument('filenames', nargs='+',
                        help='names of the files to be used for the plots; should be GCG output with STATISTICS=true, formatted as by the check-scripts for multiple instances or whole testsets')

    parsed_args = parser.parse_args(args)

    # check if the provided arguments are consistent
    if (0 < parsed_args.maxround and parsed_args.maxround < parsed_args.minround) or parsed_args.minround <= 0 or parsed_args.maxround < 0:
        print 'please make sure that 1 <= MINROUND <= MAXROUND (or MAXROUND == 0 for the last possible round)'
        exit()
    if not parsed_args.colors in plt.cm.datad:
        print 'please use a colormap that is supported by pyplot (' + parsed_args.colors + ' is not supported)'
        exit()
    if parsed_args.load and parsed_args.save:
        print 'please load OR save data'
        exit()
    if not parsed_args.save and parsed_args.bubble_only and parsed_args.summary_only:
        print 'you passed --bubble-only and --summary-only, no plot will be drawn'
        exit()

    return parsed_args

def set_params(args):
    """
    Set the global parameters from the parsed command-line arguments
    :param args: parsed command-line arguments
    :return:
    """
    params['outdir'] = args.outdir
    params['root_only'] = args.root_only
    params['details'] = not args.png
    params['nSplit'] = args.splitrounds
    params['minRound'] = args.minround
    params['maxRound'] = args.maxround
    params['colors'] = args.colors
    params['lines'] = args.lines
    params['instances'] = args.instances
    params['no_summary'] = args.no_summary or args.bubble_only
    params['no_bubble'] = args.no_bubble or args.summary_only
    params['no_complete'] = args.summary_only or args.bubble_only
    params['no_text'] = args.no_text
    params['save'] = args.save
    params['load'] = args.load

def get_colmap(pricers):
    """
    Returns a list of colors, with same length as pricers, that can be used for the bar plot
    Each pricing problem has its own color
    :param pricers: a list with pricing_problem ids
    :return: a list of colors as used by pyplot.bar()
    """

    # initialize variables
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
    Removes all texts in figure from the list texts, that overlap, by setting their visibility to False
    :param figure: the figure object to which the texts belong
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
    Make a complete plot from the structured data
    :param data: dataframe with the collected data
    :param name: the problemname
    :return:
    """
    start_time = time.time()

    # set the height of the zero bars
    ymin = -0.15

    # flat out the data again
    data = data.reset_index()

    # workaround for most times being zero (0.01s is SCIPs smallest time-interval)
    # the column pool has no bar and therefore gets no time shift
    data.loc[data.pricing_prob <> -1, 'time'] = data[data.pricing_prob <> -1].time + 0.01

    # calculate the starting time of each round
    data['starting_time'] = data.time.cumsum() - data.time

    # extract the column pool data and delete it from data
    x_colpool = data[data.pricing_prob == -1].starting_time.values
    y_colpool = data[data.pricing_prob == -1].nVars.values
    data = data[data.pricing_prob <> -1].reset_index()

    # define position, width and height of the bars; the first two are defined by time, the last by nVars
    x = data.starting_time.values
    y = (data.nVars - ymin).values
    widths = data.time.values
    colors = get_colmap(data['pricing_prob'].values)

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

    # add the column pool data as a scatter plot
    ax.scatter(x_colpool, y_colpool, color = 'green', marker = 'o', s = 100, zorder = 10)

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
    if len(y_colpool) == 0:
        ymax = 1.01 * max(y)
    # if the maximal nVars of the column pool exceeds nVars of all pricing problems, do not take the cp into account for ymax
    elif 2*max(y) < max(y_colpool) and max(y) >= 2:
        ymax = max(y) * 1.5
    else:
        ymax = 1.01 * max(y.tolist() + y_colpool.tolist())

    # formatting
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([0,x[-1]+widths[-1]])
    # keep only integer values on the y-axis
    # todo: there is a more elegant way
    old_yticks = ax.get_yticks()
    new_yticks = []
    for i,n in enumerate(old_yticks):
        if abs(n - int(n)) < 0.01 and n >= 0:
            new_yticks.append(n)
    ax.set_yticks(new_yticks)
    ax.tick_params(axis = 'both', length = textsize/2, width = textsize/40, labelsize = textsize*0.9)
    ax.set_xlabel('Time / s', size = 1.15*textsize)
    # todo: add padding to the right (to yticklabels)
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
        prev_rnd = data['pricing_round'][0]
        prev_stab = data['stab_round'][0]
        prev_x = 0
        texts = []
        # special cases: no or only (initial) farkas pricing in the plot
        if data.farkas.all():
            ax.text(ax.get_xlim()[1], .8, "\it{Initial Farkas Pricing did not end}", va = 'center', ha = 'left', rotation = 90, color = 'blue', zorder = 11, size = textsize, transform = ax.transAxes)
            farkasLine = True
        elif not data.farkas.any():
            ax.text(ax.get_xlim()[0], .8, "\it{No initial Farkas Pricing}", va = 'center', ha = 'left', rotation = 90, color = 'blue', zorder = 11, size = textsize, transform = ax.transAxes)
            farkasLine = True
        else:
            farkasLine = False

        enfLine = False
        for i in range(len(x)):
            rnd = data['pricing_round'][i]
            stab = data['stab_round'][i]
            if stab > prev_stab or rnd > prev_rnd:
                if rnd > prev_rnd:
                    # bold line for a new pricing round
                    if params['lines'] or (x[i] - prev_x)/totalTime > 0.0005 or enfLine or (not farkasLine and not data['farkas'][i]):
                        line = lines.Line2D([x[i],x[i]],[0,1],color='r',linewidth=1.0, transform = trans)
                        # blue line at the end of farkas pricing
                        if not farkasLine and not data['farkas'][i]:
                            line.set_color('blue')
                            ax.text(x[i], .8, "\it{End of initial Farkas Pricing}", va = 'center', ha = 'left', rotation = 90, color = 'blue', zorder = 11, size = textsize, transform = trans)
                            farkasLine = True
                        ax.add_line(line)
                        if (x[i] - prev_x)/totalTime > 0.0005:
                            enfLine = True
                        else:
                            enfLine = False
                    texts.append(ax.text(prev_x, 1.01, 'Round '+str(prev_rnd), rotation='vertical',va='bottom', ha='left', size = textsize, transform = trans))
                    prev_rnd = rnd
                    prev_stab = stab
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
        prev_node = data['node'][0]
        prev_x = 0
        texts = []
        for i in range(len(x)):
            node = data['node'][i]
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
    For each instance create one summary plot, which shows for each round the
    cumulated running time of the pricers in this round as well as the fraction of
    pricers, that did find a variable.
    :param data: dataframe with the collected (complete) data
    :param name: the problemname
    :return:
    """
    start_time = time.time()

    # extract summary data
    summary = pd.DataFrame()
    summary['time'] = data.groupby(level=['node','pricing_round','stab_round', 'round']).sum().time
    summary['found_frac'] = data.astype(bool).groupby(level=['node','pricing_round','stab_round', 'round']).sum().nVars/data.groupby(level=['node','pricing_round','stab_round', 'round']).count().nVars*100
    summary = summary.reset_index()

    if not data.farkas.all() and data.farkas.any():
        # get the last round of initial farkas pricing
        farkas_end = (data[data.farkas == False].reset_index()['round'].values[0] + data[data.farkas == True].reset_index()['round'].values[-1])/2.

    print '    extracted summary data:', time.time() - start_time
    start_time = time.time()

    fig,ax1 = plt.subplots()
    ax2 = ax1.twinx()

    # get the data for the plot
    x = summary[summary.stab_round <= 0]['round'].values
    y_time = summary[summary.stab_round <= 0].time.values
    y_found_frac = summary[summary.stab_round <= 0].found_frac.values
    x_stab = summary[summary.stab_round > 0]['round'].values
    y_stab_time = summary[summary.stab_round > 0].time.values
    y_stab_found_frac = summary[summary.stab_round > 0].found_frac.values

    # format the plot
    ax1.set_xlabel('Round', size='large')
    ax1.set_ylabel('Time / s', color='k', size='large')
    ax2.set_ylabel('Fraction of successfull pricing problems / \%', color='r', size='large')

    ax1.tick_params(axis='both', labelsize='large')
    ax2.tick_params(axis='both', labelsize='large')

    # todo: do this in the complete plot as well
    ax1.get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    ax1.set_xlim([0, max(x.tolist() + x_stab.tolist()) + 0.9])
    if max(y_time) > 0:
        ax1.set_ylim([-max(y_time)*0.1,max(y_time)*1.1])
    else:
        ax1.set_ylim([-0.001,0.01])
    ax2.set_ylim([-max(y_found_frac)*0.15,max(y_found_frac)*1.15])

    trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)

    p1,p2 = ax1.transData.transform([(1,0),(2,0)])
    perimeter = max([3.5,(p2[0] - p1[0])/5.])

    # plot the data
    ax1.scatter(x,y_time, color='k', s=perimeter**2)
    ax2.scatter(x,y_found_frac, color='r', s=perimeter**2)
    ax1.scatter(x_stab,y_stab_time, color='k', s=perimeter**2, marker='x', alpha=.5)
    ax2.scatter(x_stab,y_stab_found_frac, color='r', s=perimeter**2, marker='x', alpha=.5)

    # add a line after the root-node
    if summary.node.max() > 1:
        x_line = (summary[summary.node > 1].index[0] + 1 + summary[summary.node == 1].index[-1] + 1)/2.
        line = lines.Line2D([x_line,x_line],[0,1.02],color='red',linewidth=.5,linestyle='--', transform = trans)
        line.set_clip_on(False)
        ax1.add_line(line)
        ax1.text(x_line, 1.025, "\it{End of Root}", ha = 'center', size = 'smaller', color = 'red', zorder = 1, transform = trans)
    elif summary.node.max() == 1:
        ax1.text(0.96, 1.025, "\it{End of Root or Root did not end}", ha = 'right', size = 'smaller', color = 'red', zorder = 1, transform = ax1.transAxes)

    # add a line after the initial farkas pricing
    if data.farkas.all():
        ax1.text(0.96, 1.01, "\it{Initial Farkas Pricing did not end}", size = 'smaller', ha = 'right', color = 'blue', zorder = 1, transform = ax1.transAxes)
    elif not data.farkas.any():
        ax1.text(0, 1.01, "\it{No initial Farkas Pricing}", size = 'smaller', ha = 'left', color = 'blue', zorder = 1, transform = ax1.transAxes)
    else:
        x_line = farkas_end
        line = lines.Line2D([x_line,x_line],[0,1],color='blue',linewidth=.5,linestyle='--', transform = trans)
        ax1.add_line(line)
        ax1.text(x_line, 1.01, "\it{End of initial Farkas Pricing}", size = 'smaller', ha = 'center', color = 'blue', zorder = 1, transform = trans)

    print '    plotted summary:', time.time() - start_time
    start_time = time.time()

    # save the plot
    fig.set_size_inches(11.7,8.3)
    plt.tight_layout()
    fig.subplots_adjust(top = 0.96)
    if params['details']:
        plt.savefig(params['outdir'] + '/' + name + '_summary.pdf')
    else:
        plt.savefig(params['outdir'] + '/' + name + '_summary.png')
    plt.close()
    print '    saved summary:', time.time() - start_time

def make_bubble_plot(data, name):
    """
    For each instance create one bubble plot, that shows which pricing problem found a variable in each round
    :param data: dataframe with the collected (complete) data
    :param name: the problemname
    :return:
    """
    start_time = time.time()

    # flat out the data again
    data = data.reset_index()

    pricers = data[data.pricing_prob <> -1].pricing_prob.unique().tolist()
    if not data.farkas.all() and data.farkas.any():
        # get the last round of initial farkas pricing
        farkas_end = (data[data.farkas == False]['round'].values[0] + data[data.farkas == True]['round'].values[-1])/2.

    # add x and y data to plot for every pricer and the column pool
    x = {}
    y = {}
    x_stab = {}
    y_stab = {}
    for p in (pricers + [-1]):
        x[p] = data[(data.pricing_prob == p) & (data.nVars >= 1) & (data.stab_round <= 0)]['round'].values
        y[p] = [p for i in x[p]]
        x_stab[p] = data[(data.pricing_prob == p) & (data.nVars >= 1) & (data.stab_round > 0)]['round'].values
        y_stab[p] = [p for i in x_stab[p]]

    print '    extracted bubble data:', time.time() - start_time
    start_time = time.time()

    colors = get_colmap(pricers)

    fig = plt.gcf()
    ax = plt.gca()

    # format the plot
    ax.set_xlabel('Round', size='large')
    ax.set_ylabel('Pricer ID', size='large')
    ax.set_xlim([0,data['round'].max()+0.9])
    ax.set_ylim([-1.5, max(pricers) + 0.5])
    ax.get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    ax.get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.tick_params(axis='both', labelsize='large')
    # set the size of the markers, such that they do not overlap
    p1,p2,p3 = ax.transData.transform([(1,min(pricers)),(2,min(pricers)),(1,min(pricers)+1)])
    perimeter = max(3.5,min((p2[0] - p1[0])/2.2 , (p3[1] - p1[1])/2.2))

    # plot the data
    for p in pricers:
        ax.scatter(x[p],y[p], color = colors[pricers.index(p)], s=perimeter**2)
        ax.scatter(x_stab[p],y_stab[p], color = colors[pricers.index(p)], s=perimeter**2, marker = 'x', alpha = .5)
    ax.scatter(x[-1],y[-1], facecolors = 'none', edgecolors = 'green', marker = 'o', s=perimeter**2)
    ax.scatter(x_stab[-1],y_stab[-1], facecolors = 'none', edgecolors = 'green', marker = 'o', alpha = .5, s=perimeter**2)
    ax.scatter(x_stab[-1],y_stab[-1], color = 'green', s=perimeter**2, marker = 'x', alpha = .5)

    # add a line after the root-node
    if data.node.max() > 1:
        x_line = (data[data.node > 1]['round'].iloc[0] + data[data.node == 1]['round'].iloc[-1])/2.
        line = lines.Line2D([x_line,x_line],[0,1.02],color='red',linewidth=.5,linestyle='--', transform = trans)
        line.set_clip_on(False)
        ax.add_line(line)
        ax.text(x_line, 1.025, "\it{End of Root}", ha = 'center', size = 'smaller', color = 'red', zorder = 11, transform = trans)
    elif data.node.max() == 1:
        ax.text(0.96, 1.025, "\it{End of Root or Root did not end}", ha = 'right', size = 'smaller', color = 'red', zorder = 11, transform = ax.transAxes)

    # add a line after initial farkas pricing
    if data.farkas.all():
        ax.text(0.96, 1.01, "\it{Initial Farkas Pricing did not end}", size = 'smaller', ha = 'right', color = 'blue', zorder = 11, transform = ax.transAxes)
    elif not data.farkas.any():
        ax.text(0, 1.01, "\it{No initial Farkas Pricing}", size = 'smaller', ha = 'left', color = 'blue', zorder = 11, transform = ax.transAxes)
    else:
        x_line = farkas_end
        line = lines.Line2D([x_line,x_line],[0,1],color='blue',linewidth=.5,linestyle='--', transform = trans)
        ax.add_line(line)
        ax.text(x_line, 1.01, "\it{End of initial Farkas Pricing}", size = 'smaller', ha = 'center', color = 'blue', zorder = 11, transform = trans)

    print '    plotted bubble data:', time.time() - start_time
    start_time = time.time()

    # save the plot
    fig.set_size_inches(11.7,8.3)
    plt.tight_layout()
    fig.subplots_adjust(top = 0.96)
    if params['details']:
        plt.savefig(params['outdir'] + '/' + name + '_bubble.pdf')
    else:
        plt.savefig(params['outdir'] + '/' + name + '_bubble.png')
    plt.close()

    print '    saved bubble plot:', time.time() - start_time

def plots(data, name):
    """
    Master-function for plotting. Splits the data if necessary and calls all three plotting functions (or a subset, according to the params)
    :param data: collected data as dataframe
    :param name: name of the problem
    :return:
    """

    # use tex to render the text output
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    if not params['no_bubble']:
        # build the bubble plot
        make_bubble_plot(data, name)

    if not params['no_summary']:
        # build the summary plot
        make_summary_plot(data, name)

    if params['no_complete']:
        # do not build the complete plot
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
    Plots data, that was parsed and collected earlier from the generate_files() method and saved in pickle-files
    :param files: the pickle files, from which the dataframes are to load (one file per instance)
    :return:
    """
    for file in files:
        if not os.path.exists(file):
            print 'there is no file', file
            continue
        name, ext = os.path.splitext(os.path.basename(file))
        if not (ext == '.pkl'):
            print file, 'is not a pickle file'
            continue
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

def collect_data(name, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas):
    """
    Take lists containing the parsed data and structure them in a multiindex dataframe; then save or plot the data.
    All lists have to be of equal length (representing columns in a table)
    :param name: name of the instance
    :param ind_node: node, used as index
    :param ind_pricing_round: pricing round, used as index
    :param ind_stab_round: stabilization round, used as index
    :param ind_pricing_prob: pricing problem, used as index
    :param val_time: running time of the pricing problem, will be a column in the dataframe
    :param val_nVars: number of found variables, will be a column in the dataframe
    :param val_farkas: is the pricing problem part of the initial farkas pricing? Will be a column in the dataframe
    :return:
    """
    index = pd.MultiIndex.from_arrays([ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob],
                                      names=["node", "pricing_round", "stab_round", "round", "pricing_prob"])
    data = {'time': val_time, 'nVars': val_nVars, 'farkas': val_farkas}
    df = pd.DataFrame(data = data, index = index)

    # save or plot the data
    if params['save']:
        start_time = time.time()
        df.to_pickle(params['outdir'] + '/' + name + '.pkl')
        print '    total saving:', time.time() - start_time
    else:
        start_time = time.time()
        plots(df, name)
        print '    total plotting:', time.time() - start_time

def parse_files(files):
    """
    Parse the (out-)files and structure the pricing-data in a dataframe
    :param files: List of files to be parsed
    :return:
    """
    for file in files:
        with open(file) as _file:
            # initialize all index-lists
            ind_node = []
            ind_pricing_round = []
            ind_stab_round = []
            ind_round = []
            ind_pricing_prob = []

            # initialize the value-lists
            val_time = []
            val_nVars = []
            val_farkas = []

            # initialize all counters
            node = 0
            pricing_round = 0
            stab_round = 0
            pricing_prob = 0
            round_counter = 0

            # initialize all other variables
            problemFileName = None
            farkasDone = False
            done = False

            for line in _file:
                if line.find("@0")<>-1:
                    # if the file is an out-file, generated by the check-script, reset the variables whenever a new instance starts
                    if line.startswith("@01"):
                        # print message, if the previous problem is not done yet
                        if not done and problemFileName:
                            print '    ended abruptly'
                            collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
                            print '    leaving', problemFileName
                            done = True

                        # initialize all index-lists
                        ind_node = []
                        ind_pricing_round = []
                        ind_stab_round = []
                        ind_round = []
                        ind_pricing_prob = []

                        # initialize the value-lists
                        val_time = []
                        val_nVars = []
                        val_farkas = []

                        # initialize all counters
                        node = 0
                        pricing_round = 0
                        stab_round = 0
                        pricing_prob = 0
                        round_counter = 0

                        # initialize all other variables
                        problemFileName = None
                        farkasDone = False
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

                # end of initial farkas pricing
                elif line.startswith("Starting reduced cost pricing..."):
                    farkasDone = True

                # pricer statistics end
                elif line.startswith("SCIP Status        :"):
                    # continue if no data is found
                    if not ind_node or not ind_pricing_round or not ind_stab_round or not ind_pricing_prob or not val_time or not val_nVars:
                        print '    no pricing data found'
                        done = True
                        continue

                    collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)

                    done = True

                    print '    leaving', problemFileName
                    continue

                # ignore all other lines, that do not contain pricer statistics messages
                elif not line.startswith("[src/pricer_gcg.cpp:"):
                    continue

                # extract the pricing-statistics message
                message = line.split("] statistic: ")[-1]

                if message.startswith("New pricing round at node") or message.startswith("New pr, node"):
                    try:
                        node = int(message.split()[-1])
                        pricing_round += 1
                        round_counter += 1
                        stab_round = 0
                    except ValueError:
                        print '    ended abruptly'
                        collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
                        print '    leaving', problemFileName
                        done = True
                        continue

                elif message.startswith("Stabilization round ") or message.startswith("Sr "):
                    try:
                        stab_round = int(message.split()[-1])
                        round_counter += 1
                    except ValueError:
                        print '    ended abruptly'
                        collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
                        print '    leaving', problemFileName
                        done = True
                        continue

                # todo: include non-sparse output
                # assumption: if one or more improving columns have been found in the column pool, no pricing problems were solved. todo: check that
                elif message.startswith("cp: "):
                    try:
                        if int(message.split()[1]) > 0:
                            # check if the column pool output should be included in the data
                            if (params['root_only'] and node > 1) or pricing_round < params['minRound'] or (0 < params['maxRound'] < pricing_round):
                                continue

                            # store all indices
                            ind_node.append(node)
                            ind_pricing_round.append(pricing_round)
                            ind_stab_round.append(stab_round)
                            ind_round.append(round_counter)
                            # the column pool is represented as a pricing problem with ID -1
                            ind_pricing_prob.append(-1)

                            # store the data
                            val_time.append(0.0)
                            val_nVars.append(int(message.split()[1]))
                            val_farkas.append(not farkasDone)
                    except ValueError:
                        print '    ended abruptly'
                        collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
                        print '    leaving', problemFileName
                        done = True
                        continue

                elif message.startswith("Pricing prob ") or message.startswith("P p "):
                    try:
                        pricing_prob = int(message.split()[2])

                        # check if the pricing prob should be included in the data
                        if (params['root_only'] and node > 1) or pricing_round < params['minRound'] or (0 < params['maxRound'] < pricing_round):
                            continue

                        # store all indices
                        ind_node.append(node)
                        ind_pricing_round.append(pricing_round)
                        ind_stab_round.append(stab_round)
                        ind_round.append(round_counter)
                        ind_pricing_prob.append(pricing_prob)

                        # store the data
                        val_time.append(float(message.split()[-1]))
                        if message.startswith("P p "):
                            val_nVars.append(int(message.split()[-3]))
                        else:
                            val_nVars.append(int(message.split()[5]))
                        val_farkas.append(not farkasDone)
                    except ValueError:
                        print '    ended abruptly'
                        collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
                        print '    leaving', problemFileName
                        done = True
                        continue

            if not done:
                collect_data(problemFileName, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)

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
        parse_files(parsed_args.filenames)
        if params['save']:
            logfile = open(params['outdir'] + '/origin.log','w')
            for f in parsed_args.filenames:
                logfile.write(f + '\n')
            logfile.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M'))

# Calling main script
if __name__ == '__main__':
    main()
