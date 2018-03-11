#!/usr/bin/env python2

import matplotlib
matplotlib.use('AGG')

import sys
import os
import argparse
import time
import datetime
from collections import OrderedDict

import pandas as pd
import pickle

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import lines
from matplotlib import transforms
from matplotlib import ticker
from matplotlib import patches as mpatches

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
                        help='names of the instances to be included in the plot/data-collection (can be only part of the name; default is all instances in FILENAMES)')

    parser.add_argument('-a', '--complete-only', action='store_true',
                        help='create only complete plots')

    parser.add_argument('-A', '--no-complete', action='store_true',
                        help='create no complete plots')

    parser.add_argument('-z', '--summary-only', action='store_true',
                        help='create only summary plots')

    parser.add_argument('-Z', '--no-summary', action='store_true',
                        help='create no summary plots')

    parser.add_argument('-b', '--bubble-only', action='store_true',
                        help='create only bubble plots')

    parser.add_argument('-B', '--no-bubble', action='store_true',
                        help='create no bubble plots')

    parser.add_argument('-t', '--time-only', action='store_true',
                        help='create only time plots')

    parser.add_argument('-T', '--no-time', action='store_true',
                        help='create no time plots')

    parser.add_argument('-g', '--gap-only', action='store_true',
                        help='create only gap plots')

    parser.add_argument('-G', '--no-gap', action='store_true',
                        help='create no gap plots')

    parser.add_argument('-n', '--no-text', action='store_true',
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
    if not parsed_args.save and parsed_args.no_bubble and parsed_args.no_summary and parsed_args.no_complete and parsed_args.no_time and parsed_args.no_gap:
        print 'based on the passed parameters, no plot will be drawn'
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
    params['no_summary'] = args.no_summary or args.bubble_only or args.time_only or args.gap_only or args.complete_only
    params['no_bubble'] = args.no_bubble or args.summary_only or args.time_only or args.gap_only or args.complete_only
    params['no_time'] = args.no_time or args.summary_only or args.bubble_only or args.gap_only or args.complete_only
    params['no_gap'] = args.no_gap or args.summary_only or args.bubble_only or args.time_only or args.complete_only
    params['no_complete'] = args.no_complete or args.summary_only or args.bubble_only or args.time_only or args.gap_only
    params['no_text'] = args.no_text
    params['save'] = args.save
    params['load'] = args.load

def get_colmap(pricers):
    """
    Returns a list of colors, with same length as pricers, that can be used for the bar plot
    Also returns a mapping color -> pricer_id, used for legends
    Each pricing problem has its own color
    :param pricers: a list with pricing_problem ids
    :return: a list of colors as used by pyplot.bar()
    """
    # build a list of pricer ids, in which each pricer id appears once and sort it by id
    pricer_ids = []
    for p in pricers:
        if not p in pricer_ids:
            pricer_ids.append(p)
    pricer_ids = sorted(pricer_ids)

    # get a color map of the right length, so that each color-id gets its own color
    cmap = plt.get_cmap(params['colors'],len(pricer_ids))

    # build the mapping
    mapping = OrderedDict()
    for p in pricer_ids:
        mapping[p] = cmap(pricer_ids.index(p))

    # build a list of colors
    colors = [mapping[p] for p in pricers]

    return colors, mapping

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

def get_y1_in_ax(obj, figure):
    """
    Calculates the upper end of the object obj (e.g. texts) in figure in axes coordinates
    :return: upper end of the object as axes coordinate
    """
    return obj.get_window_extent(renderer = figure.canvas.get_renderer()).transformed(figure.gca().transAxes.inverted()).y1

def get_x1_in_data(obj, figure):
    """
    Calculates the right end of the object obj (e.g. texts) in figure in data coordinates
    :return: right end of the object as data coordinate
    """
    return obj.get_window_extent(renderer = figure.canvas.get_renderer()).transformed(figure.gca().transData.inverted()).x1

def get_info_from_filename(filename):
    """
    Parses a filename and extracts the instance name, settings & scip_status
    :param filename: the name of a pickle-file without extension
    :return: info dictionary
    """
    info = {}
    info['status'] = filename.split('.')[-1]
    info['settings'] = filename.split('.')[-2]
    info['instance'] = filename[:-(len(info['settings']) + len(info['status']) + 2)]
    return info

def get_filename_from_info(info):
    """
    Parses a filename and extracts the instance name, settings & scip_status
    :param filename: the name of a pickle-file without extension
    :return: info dictionary
    """
    filename = str(info['instance']) + '.' + info['settings'] + '.' + info['status']
    if 'minRound' in info:
        filename += '.from' + info['minRound']
        if 'maxRound' in info:
            filename += 'to' + info['maxRound']
    return filename

def make_plot(data, info):
    """
    Make a complete plot from the structured data
    :param data: dataframe with the collected data
    :param info: dictionary containing information about the data like the name of the instance, the settings & the scip_status
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

    # set the height of the LP time bars to a maximum value
    data.loc[data.pricing_prob == -2, 'nVars'] = data.nVars.max() * 10

    # extract the column pool data and delete it from data
    x_colpool = data[data.pricing_prob == -1].starting_time.values
    y_colpool = data[data.pricing_prob == -1].nVars.values
    data = data[data.pricing_prob <> -1].reset_index()

    # define position, width and height of the bars; the first two are defined by time, the last by nVars
    x = data.starting_time.values
    y = (data.nVars - ymin).values
    widths = data.time.values
    colors, cmapping = get_colmap(data.pricing_prob.values)

    # sometimes we need just the pricers
    y_pricers = (data[data.pricing_prob >= 0].nVars - ymin).values

    print '    data restructured:', time.time() - start_time
    start_time = time.time()

    # make the bar plot
    if params['details']:
        lw = 0.01
    else:
        lw = 1.0
    plt.bar(x, y, widths, bottom = ymin, align = 'edge', linewidth = lw, edgecolor = 'white', color = colors, label='pricing problems')
    fig = plt.gcf()
    ax = plt.gca()

    # add the column pool data as a scatter plot
    cp_scatter = ax.scatter(x_colpool, y_colpool, color = 'green', marker = 'o', s = 100, zorder = 10, label = 'column pool')

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
        ymax = 1.01 * max(y_pricers)
    # if the maximal nVars of the column pool exceeds nVars of all pricing problems, do not take the cp into account for ymax
    elif 2*max(y_pricers) < max(y_colpool) and max(y_pricers) >= 2:
        ymax = max(y_pricers) * 1.5
    else:
        ymax = 1.01 * max(y_pricers.tolist() + y_colpool.tolist())
    xmin = 0
    xmax = x[-1]+widths[-1]

    # formatting
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([xmin,xmax])
    ax.get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True, nbins = 15))
    ax.tick_params(axis = 'both', length = textsize/2, width = textsize/40, labelsize = textsize*0.9, pad = 15)
    ax.set_xlabel('Time / s', size = 1.15*textsize)
    ax.set_ylabel('\# of variables', size = 1.15*textsize)
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

    print '    data formatted:', time.time() - start_time
    start_time = time.time()

    if params['no_text']:
        # save the figure
        filename = params['outdir'] + '/' + get_filename_from_info(info)
        n = ''
        while os.path.isfile(filename + n):
            if n == '':
                n = 1
            else:
                n += 1
        if params['details']:
            plt.savefig(filename + n + '.pdf')
        else:
            plt.savefig(filename + n + '.png')

        print '    save:', time.time() - start_time
        print '    saved figure'
    else:
        # special cases: no or only (initial) farkas pricing in the plot
        if data.farkas.all():
            ax.text(.991, .99, "\it{Initial Farkas Pricing did not end}", va = 'top', ha = 'right', rotation = 0, color = 'blue', zorder = 11, size = textsize * .95, transform = ax.transAxes, bbox=dict(facecolor = 'white', edgecolor = 'none', alpha = .85, pad = 20))
            farkasLine = True
        elif not data.farkas.any():
            ax.text(.009, .99, "\it{No initial Farkas Pricing}", va = 'top', ha = 'left', rotation = 0, color = 'blue', zorder = 11, size = textsize * .95, transform = ax.transAxes, bbox=dict(facecolor = 'white', edgecolor = 'none', alpha = .85, pad = 20))
            farkasLine = True
        else:
            # add a line at the end of farkas pricing in the loop below
            farkasLine = False

        # add information about the stabilization & pricing rounds
        prev_rnd = data['pricing_round'][0]
        prev_x = 0
        prev_x_drawn = 0
        texts = []

        texts.append(ax.text(-0.001, 1.01, '\\textbf{Round}', rotation=0,va='bottom', ha='right', size = textsize*.75, transform = ax.transAxes))
        for pos,rnd,far in zip(data.drop_duplicates(['pricing_round','stab_round']).starting_time, data.drop_duplicates(['pricing_round','stab_round']).pricing_round, data.drop_duplicates(['pricing_round','stab_round']).farkas):
            if rnd > prev_rnd:
                # bold line for a new pricing round
                if params['lines'] or (pos - prev_x_drawn)/totalTime > 0.002 or (not farkasLine and not far):
                    line = lines.Line2D([pos,pos],[0,1],color='r',linewidth=1.0, transform = trans)
                    # blue line at the end of farkas pricing
                    if not farkasLine and not far:
                        line.set_color('blue')
                        if pos <= (xmax + xmin) / 2:
                            align = 'left'
                        else:
                            align = 'right'
                        ax.text(pos, .99, "\it{End of initial Farkas Pricing}", va = 'top', ha = align, rotation = 0, color = 'blue', zorder = 11, size = textsize * .95, transform = trans, bbox=dict(facecolor = 'white', edgecolor = 'none', alpha = .85, pad = 20))
                        farkasLine = True
                    ax.add_line(line)
                    prev_x_drawn = pos
                # write the round number, if there is space for it
                if len(texts) == 0 or get_x1_in_data(texts[-1], fig) < prev_x:
                    texts.append(ax.text(prev_x, 1.01, str(prev_rnd), rotation='vertical',va='bottom', ha='left', size = textsize, transform = trans))
                prev_rnd = rnd
                prev_x = pos
            else:
                # dashed line for a new stabilization round
                line = lines.Line2D([pos,pos],[0,1],color='orange',linestyle='--',linewidth=0.8, transform = trans)
                ax.add_line(line)
        if len(texts) == 0 or get_x1_in_data(texts[-1], fig) < prev_x:
            texts.append(ax.text(prev_x, 1.01, str(prev_rnd), rotation='vertical',va='bottom', ha='left', size = textsize, transform = trans))
        text_height = max(get_y1_in_ax(texts[0], fig),get_y1_in_ax(texts[-1], fig))

        print '    stab- and pricing-round information:', time.time() - start_time
        start_time = time.time()

        # add information about the nodes
        prev_node = data['node'][0]
        prev_x = 0
        node_header_x = get_x1_in_data(texts[0], fig)
        text_height += 0.0006
        texts = []
        texts.append(ax.text(node_header_x, text_height+0.001, '\\textbf{Node}', ha='right', size = textsize*.75, transform = trans))
        for pos, node in zip(x, data.node):
            if node > prev_node:
                line = lines.Line2D([pos,pos],[1,text_height],color='r',linewidth=1.0, transform = trans)
                line.set_clip_on(False)
                ax.add_line(line)
                # write the node number, if there is space for it
                if len(texts) == 0 or get_x1_in_data(texts[-1], fig) < prev_x:
                    texts.append(ax.text(prev_x, text_height, str(prev_node), ha='left', size = textsize, transform = trans))
                prev_node = node
                prev_x = pos
        if len(texts) == 0 or get_x1_in_data(texts[-1], fig) < prev_x:
            texts.append(ax.text(prev_x, text_height, str(prev_node), ha='left', size = textsize, transform = trans))
        text_height = get_y1_in_ax(texts[-1], fig)

        print '    node information:', time.time() - start_time
        start_time = time.time()

        # draw a legend, but do not include more than 25 pricing problems
        patches = [mpatches.Patch(color = cmapping[p], label = 'pricing problem ' + str(p)) for p in cmapping]
        patches[0].set_label('Master LP Time')
        if len(patches) > 31:
            patches = patches[:31] + [mpatches.Patch(color = 'white', alpha = 0, label = '...')]
        handles = patches + [lines.Line2D([0,0], [0,1], color = 'red', linewidth = 2., label = 'pricing round'), lines.Line2D([0,0], [0,1], color = 'orange', linestyle = '--', linewidth = 1.6, label = 'stabilization round'), cp_scatter]
        plt.legend(handles = handles, bbox_to_anchor = (1.02, .915), loc = 2, fontsize = textsize)

        # add other information
        name = info['instance']
        settings = info['settings']
        status = info['status']
        if len(name) <= 12:
            ax.text(1.093, (text_height + 1.)/2., '\\textbf{\\underline{' + name.replace('_','\_') + '}}', ha = 'center', va = 'center', size = 1.3 * textsize, transform = ax.transAxes)
        else:
            ax.text(1.093, (text_height + 1.)/2., '\\textbf{\\underline{' + name[:11].replace('_','\_') + '\\dots}}', ha = 'center', va = 'center', size = 1.3 * textsize, transform = ax.transAxes)
        ax.text(1.093, .96, '\\textbf{Settings:} \n' + settings.replace('_','\_') + '\n \\textbf{SCIP Status:} \n' + status.replace('_',' '), ha = 'center', va = 'center', size = .9 * textsize, transform = ax.transAxes)
        name = name + '.' + settings

        # save the figure
        plt.tight_layout()
        fig.subplots_adjust(top = 0.98/text_height, right = 0.85, left = 0.03)
        filename = params['outdir'] + '/' + name
        n = ''
        dot = ''
        if params['details']:
            while os.path.isfile(filename + dot + str(n) + '.pdf'):
                if n == '':
                    n = 2
                    dot = '.'
                else:
                    n += 1
            plt.savefig(filename + dot + str(n) + '.pdf')
        else:
            while os.path.isfile(filename + dot + str(n) + '.png'):
                if n == '':
                    n = 2
                    dot = '.'
                else:
                    n += 1
            plt.savefig(filename + dot + str(n) + '.png')
        plt.close()

        print '    save:', time.time() - start_time
        print '    saved figure'

def make_summary_plot(data, info):
    """
    For each instance create one summary plot, which shows for each round the
    cumulated running time of the pricers in this round as well as the fraction of
    pricers, that did find a variable.
    :param data: dataframe with the collected (complete) data
    :param info: dictionary containing information about the data like the name of the instance, the settings & the scip_status
    :return:
    """
    start_time = time.time()

    # extract summary data
    data = data.query('pricing_prob <> -2')
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
    x_mean = summary['round'].values
    y_mean_time = summary.time.rolling(int(.05*(max(x_mean) - min(x_mean))), center = True).mean()
    y_mean_found_frac = summary.found_frac.rolling(int(.05*(max(x_mean) - min(x_mean))), center = True).mean()

    # format the plot
    ax1.set_xlabel('Pricing Round', size='large')
    ax1.set_ylabel('Time / s', color='k', size='large')
    ax2.set_ylabel('Fraction of successfull pricing problems / \%', color='r', size='large')

    ax1.tick_params(axis='both', labelsize='large')
    ax2.tick_params(axis='both', labelsize='large')

    # set the axes limits
    xmin = 0
    xmax = max(x.tolist() + x_stab.tolist()) + 0.9
    ax1.set_xlim([xmin, xmax])
    if max(y_time) > 0:
        ax1.set_ylim([-max(y_time)*0.1,max(y_time)*1.1])
    else:
        ax1.set_ylim([-0.001,0.01])
    ax2.set_ylim([-max(y_found_frac)*0.15,max(y_found_frac)*1.15])

    # make the xticks
    roundsDF = summary.drop_duplicates('pricing_round').reset_index()[['pricing_round','round']]
    xtickpositions = []
    xticklabels = []
    deltaPosMin = int(summary['round'].max() / 20.00001)
    prev_pos = - deltaPosMin
    for label, pos in zip(roundsDF['pricing_round'].tolist(),roundsDF['round'].tolist()):
        if pos - prev_pos > deltaPosMin:
            xtickpositions.append(pos)
            xticklabels.append(str(label))
            prev_pos = pos
    ax1.set_xticks(xtickpositions)
    ax2.set_xticks(xtickpositions)
    ax1.set_xticklabels(xticklabels)
    ax2.set_xticklabels(xticklabels)
    del roundsDF

    trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)

    p1,p2 = ax1.transData.transform([(1,0),(2,0)])
    perimeter = max([3.5,(p2[0] - p1[0])/5.])

    # plot the data
    ax1.scatter(x,y_time, color='k', s=perimeter**2)
    ax2.scatter(x,y_found_frac, color='r', s=perimeter**2)
    ax1.scatter(x_stab,y_stab_time, color='k', s=perimeter**2, marker='x', alpha=.5)
    ax2.scatter(x_stab,y_stab_found_frac, color='r', s=perimeter**2, marker='x', alpha=.5)
    ax1.plot(x_mean,y_mean_time, 'k--')
    ax2.plot(x_mean,y_mean_found_frac, 'r--')

    # add a line after the root-node
    if summary.node.max() > 1:
        x_line = (summary[summary.node > 1].index[0] + 1 + summary[summary.node == 1].index[-1] + 1)/2.
        line = lines.Line2D([x_line,x_line],[0,1.02],color='red',linewidth=.5,linestyle='--', transform = trans)
        line.set_clip_on(False)
        ax1.add_line(line)
        if (x_line - xmin) < 0.1*(xmax - xmin):
            align = 'left'
            x_line = xmin
        elif (xmax - x_line) < 0.1*(xmax - xmin):
            align = 'right'
            x_line = xmax
        else:
            align = 'center'
        ax1.text(x_line, 1.025, "\it{End of Root}", ha = align, size = 'smaller', color = 'red', zorder = 1, transform = trans)
    elif summary.node.max() == 1:
        ax1.text(1., 1.025, "\it{Root did not end}", ha = 'right', size = 'smaller', color = 'red', zorder = 1, transform = ax1.transAxes)

    # add a line after the initial farkas pricing
    if data.farkas.all():
        ax1.text(1., 1.01, "\it{Initial Farkas Pricing did not end}", size = 'smaller', ha = 'right', color = 'blue', zorder = 1, transform = ax1.transAxes)
    elif not data.farkas.any():
        ax1.text(0, 1.01, "\it{No initial Farkas Pricing}", size = 'smaller', ha = 'left', color = 'blue', zorder = 1, transform = ax1.transAxes)
    else:
        x_line = farkas_end
        line = lines.Line2D([x_line,x_line],[0,1],color='blue',linewidth=.5,linestyle='--', transform = trans)
        ax1.add_line(line)
        if (x_line - xmin) < 0.1*(xmax - xmin):
            align = 'left'
            x_line = xmin
        elif (xmax - x_line) < 0.1*(xmax - xmin):
            align = 'right'
            x_line = xmax
        else:
            align = 'center'
        ax1.text(x_line, 1.01, "\it{End of initial Farkas Pricing}", size = 'smaller', ha = align, color = 'blue', zorder = 1, transform = trans)

    # draw a legend
    handles = []
    handles.append(lines.Line2D([0,0], [0,1], color = 'None', marker = 'o', markerfacecolor = 'k', markeredgecolor = 'k', markersize = 5, label = 'Pricing Round'))
    handles.append(lines.Line2D([0,0], [0,1], color = 'None', marker = 'x', markerfacecolor = 'k', markeredgecolor = 'k', markersize = 5, alpha = .7, label = 'Stabilization Round'))
    plt.legend(handles = handles, loc = 3, bbox_to_anchor = (.0, 1.04, .5, 1.04), ncol = 2, mode = 'expand')

    # add other information
    ax1.text(.5, 1.11, '\\textbf{\\underline{' + info['instance'].replace('_','\_') + '}}', ha = 'center', size = 'large', transform = ax1.transAxes)
    ax1.text(.75, 1.055, '\\textbf{Settings:} \\textit{' + info['settings'].replace('_','\_') + '}', ha = 'right', size = 'medium', transform = ax1.transAxes)
    ax1.text(1., 1.055, '\\textbf{SCIP Status:} \\textit{' + info['status'].replace('_',' ') + '}', ha = 'right', size = 'medium', transform = ax1.transAxes)

    print '    plotted summary:', time.time() - start_time
    start_time = time.time()

    # save the plot
    fig.set_size_inches(11.7,8.3)
    plt.tight_layout()
    fig.subplots_adjust(top = 0.87)
    filename = params['outdir'] + '/' + info['instance'] + '.summary.' + info['settings']
    n = ''
    dot = ''
    if params['details']:
        while os.path.isfile(filename + dot + str(n) + '.pdf'):
            if n == '':
                n = 2
                dot = '.'
            else:
                n += 1
        plt.savefig(filename + dot + str(n) + '.pdf')
    else:
        while os.path.isfile(filename + dot + str(n) + '.png'):
            if n == '':
                n = 2
                dot = '.'
            else:
                n += 1
        plt.savefig(filename + dot + str(n) + '.png')
    plt.close()
    print '    saved summary:', time.time() - start_time

def make_bubble_plot(data, info):
    """
    For each instance create one bubble plot, that shows which pricing problem found a variable in each round
    :param data: dataframe with the collected (complete) data
    :param info: dictionary containing information about the data like the name of the instance, the settings & the scip_status
    :return:
    """
    start_time = time.time()

    # flat out the data again
    data = data.query('pricing_prob <> -2').reset_index()

    pricer_min = data.pricing_prob.min()
    pricer_max = data.pricing_prob.max()
    if not data.farkas.all() and data.farkas.any():
        # get the last round of initial farkas pricing
        farkas_end = (data[data.farkas == False]['round'].values[0] + data[data.farkas == True]['round'].values[-1])/2.

    # add x and y data to plot for every pricer and the column pool
    x = []
    y = []
    colors = []
    x_stab = []
    y_stab = []
    colors_stab = []
    cmapping = get_colmap(data[data.pricing_prob <> -1].pricing_prob.unique())[1]
    cmapping[-1] = 'green'
    bubbleDF = data[(data.nVars >= 1) & (data.stab_round <= 0)][['round','pricing_prob','nVars']].reset_index()
    bubbleDF_stab = data[(data.nVars >= 1) & (data.stab_round > 0)][['round','pricing_prob','nVars']].reset_index()
    nVars = {}
    nVars_total = 0
    pricers = data.pricing_prob.unique()
    for p in pricers:
        tmp = bubbleDF[bubbleDF.pricing_prob == p]['round'].tolist()
        x = x + tmp
        y = y + [p for i in tmp]
        colors = colors + [cmapping[p] for i in tmp]
        tmp = bubbleDF_stab[bubbleDF_stab.pricing_prob == p]['round'].tolist()
        x_stab = x_stab + tmp
        y_stab = y_stab + [p for i in tmp]
        colors_stab = colors_stab + [cmapping[p] for i in tmp]
        nVars[p] = bubbleDF[bubbleDF.pricing_prob == p].nVars.sum()
        nVars[p] += bubbleDF_stab[bubbleDF_stab.pricing_prob == p].nVars.sum()
        nVars_total += nVars[p]
    y_bar = sorted(pricers.tolist())
    x_bar = [100*float(nVars[p])/nVars_total for p in y_bar]
    del bubbleDF, bubbleDF_stab

    print '    extracted bubble data:', time.time() - start_time
    start_time = time.time()

    fig = plt.gcf()
    gs = gridspec.GridSpec(1,2,width_ratios=[7,1], wspace = 0.05)
    ax = plt.subplot(gs[0])
    ax_bar = plt.subplot(gs[1], sharey = ax)

    # format the plot
    ax.set_xlabel('Pricing Round', size='large')
    ax.set_ylabel('Pricing Problem ID', size='large')
    xmin = 0
    xmax = data['round'].max()+0.9
    ax.set_xlim([xmin,xmax])
    if pricer_max == pricer_min:
        ax.set_ylim([pricer_min - 1, pricer_max + 1])
    else:
        ax.set_ylim([pricer_min - 0.5, pricer_max + 0.5])
    ax.get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.tick_params(axis='both', labelsize='large')
    ax_bar.set_xlabel('\% of found variables', size='large')
    ax_bar.tick_params(axis='x', labelsize='large')
    ax_bar.get_yaxis().set_visible(False)

    # make the xticks
    roundsDF = data.drop_duplicates('pricing_round').reset_index()[['pricing_round','round']]
    xtickpositions = []
    xticklabels = []
    deltaPosMin = int(data['round'].max() / 20.00001)
    prev_pos = - deltaPosMin
    for label, pos in zip(roundsDF['pricing_round'].tolist(),roundsDF['round'].tolist()):
        if pos - prev_pos > deltaPosMin:
            xtickpositions.append(pos)
            xticklabels.append(str(label))
            prev_pos = pos
    ax.set_xticks(xtickpositions)
    ax.set_xticklabels(xticklabels)
    del roundsDF

    # set the size of the markers, such that they do not overlap
    p1,p2,p3 = ax.transData.transform([(1,pricer_min),(2,pricer_min),(1,pricer_min+1)])
    perimeter = max(3.5,min((p2[0] - p1[0])/2.2 , (p3[1] - p1[1])/2.2))

    # plot the data
    ax.scatter(x,y, color = colors, s=perimeter**2)
    ax.scatter(x_stab,y_stab, color = colors_stab, s=perimeter**2, marker = 'x', alpha = .5)

    # add a line after the root-node
    if data.node.max() > 1:
        x_line = (data[data.node > 1]['round'].iloc[0] + data[data.node == 1]['round'].iloc[-1])/2.
        line = lines.Line2D([x_line,x_line],[0,1.02],color='red',linewidth=.5,linestyle='--', transform = trans)
        line.set_clip_on(False)
        ax.add_line(line)
        if (x_line - xmin) < 0.1*(xmax - xmin):
            align = 'left'
            x_line = xmin
        elif (xmax - x_line) < 0.1*(xmax - xmin):
            align = 'right'
            x_line = xmax
        else:
            align = 'center'
        ax.text(x_line, 1.025, "\it{End of Root}", ha = align, size = 'smaller', color = 'red', zorder = 11, transform = trans)
    elif data.node.max() == 1:
        ax.text(1., 1.025, "\it{Root did not end}", ha = 'right', size = 'smaller', color = 'red', zorder = 11, transform = ax.transAxes)

    # add a line after initial farkas pricing
    if data.farkas.all():
        ax.text(1., 1.01, "\it{Initial Farkas Pricing did not end}", size = 'smaller', ha = 'right', color = 'blue', zorder = 11, transform = ax.transAxes)
    elif not data.farkas.any():
        ax.text(0, 1.01, "\it{No initial Farkas Pricing}", size = 'smaller', ha = 'left', color = 'blue', zorder = 11, transform = ax.transAxes)
    else:
        x_line = farkas_end
        line = lines.Line2D([x_line,x_line],[0,1],color='blue',linewidth=.5,linestyle='--', transform = trans)
        ax.add_line(line)
        if (x_line - xmin) < 0.1*(xmax - xmin):
            align = 'left'
            x_line = xmin
        elif (xmax - x_line) < 0.1*(xmax - xmin):
            align = 'right'
            x_line = xmax
        else:
            align = 'center'
        ax.text(x_line, 1.01, "\it{End of initial Farkas Pricing}", size = 'smaller', ha = align, color = 'blue', zorder = 11, transform = trans)

    # draw the bar plot
    if pricer_max - pricer_min >= 4:
        height = .95
    else:
        height = .1
    ax_bar.barh(y_bar, x_bar, align = 'center', height = height, color = [cmapping[p] for p in y_bar])

    print '    plotted bubble data:', time.time() - start_time
    start_time = time.time()

    # draw a legend
    handles = []
    handles.append(lines.Line2D([0,0], [0,1], color = 'None', marker = 'o', markerfacecolor = 'k', markersize = 5, label = 'Pricer has found at least one variable'))
    handles.append(lines.Line2D([0,0], [0,1], color = 'None', marker = 'o', markerfacecolor = 'green', markersize = 5, label = 'Variables were taken from column pool (ID -1)'))
    handles.append(lines.Line2D([0,0], [0,1], color = 'None', marker = 'x', markerfacecolor = 'k', markeredgecolor = 'k', markersize = 5, alpha = .5, label = 'Pricer has found at least one variable in stab. round'))
    ax.legend(handles = handles, loc = 3, bbox_to_anchor = (.0, 1.04, 1.18, .02), ncol = 3, mode = 'expand')

    # add other information
    ax.text(.3, 1.11, '\\textbf{\\underline{' + info['instance'].replace('_','\_') + '}}', ha = 'center', size = 'large', transform = transforms.blended_transform_factory(fig.transFigure, ax.transAxes))
    ax.text(.75, 1.11, '\\textbf{Settings:} \\textit{' + info['settings'].replace('_','\_') + '}', ha = 'right', size = 'medium', transform = ax.transAxes)
    ax.text(1., 1.11, '\\textbf{SCIP Status:} \\textit{' + info['status'].replace('_',' ') + '}', ha = 'right', size = 'medium', transform = ax.transAxes)

    # save the plot
    fig.set_size_inches(11.7,8.3)
    gs.tight_layout(fig,rect = (0,0,1,.9))
    filename = params['outdir'] + '/' + info['instance'] + '.bubble.' + info['settings']
    n = ''
    dot = ''
    if params['details']:
        while os.path.isfile(filename + dot + str(n) + '.pdf'):
            if n == '':
                n = 2
                dot = '.'
            else:
                n += 1
        plt.savefig(filename + dot + str(n) + '.pdf')
    else:
        while os.path.isfile(filename + dot + str(n) + '.png'):
            if n == '':
                n = 2
                dot = '.'
            else:
                n += 1
        plt.savefig(filename + dot + str(n) + '.png')
    plt.close()

    print '    saved bubble plot:', time.time() - start_time

def make_time_plot(data, info):
    """
    For each instance create a pie chart summarizing the computing time distribution
    :param data: dataframe with the collected (complete) data
    :param info: dictionary containing information about the data like the name of the instance, the settings & the scip_status
    :return:
    """
    start_time = time.time()

    # calculate times for the total summary
    farkas_time = data.query('(pricing_prob >= 0) & (farkas == True)').time.sum() * 100
    redcost_time = data.query('(pricing_prob >= 0) & (farkas == False)').time.sum() * 100
    masterlp_time = data.query('pricing_prob == -2').time.sum() * 100

    min_angle = 15./360.

    # calculate times for the pricer summary
    df = data.query('pricing_prob >= 0').reset_index()[['pricing_prob','time','nVars']].groupby('pricing_prob').sum()
    df['colors'] = get_colmap(df.index.tolist())[0]
    df = df[df.time >= 0.01].sort_values('time', ascending = False)
    pricer_times = df.time[df.time / df.time.sum() >= min_angle].tolist()
    pricer_times_labels = df[df.time / df.time.sum() >= min_angle].index.astype('str').tolist()
    pricer_times_colors = df.colors[df.time / df.time.sum() >= min_angle].tolist()
    if df.time[df.time / df.time.sum() < min_angle].sum() > 0:
        pricer_times.append(df.time[df.time / df.time.sum() < min_angle].sum())
        pricer_times_labels.append('Others')
        pricer_times_colors.append('grey')

    # calculate number of found variables
    pricer_nVars = df.nVars[df.nVars / df.nVars.sum() >= min_angle].tolist()
    pricer_nVars_labels = df[df.nVars / df.nVars.sum() >= min_angle].index.astype('str').tolist()
    pricer_nVars_colors = df.colors[df.nVars / df.nVars.sum() >= min_angle].tolist()
    if df.nVars[df.nVars / df.nVars.sum() < min_angle].sum() > 0:
        pricer_nVars.append(df.nVars[df.nVars / df.nVars.sum() < min_angle].sum())
        pricer_nVars_labels.append('Others')
        pricer_nVars_colors.append('grey')

    # calculate efficiencies
    df['efficiency'] = df.nVars / df.time
    pricer_efficiencies = df.efficiency[df.efficiency / df.efficiency.sum() >= min_angle].tolist()
    pricer_efficiencies_labels = df.efficiency[df.efficiency / df.efficiency.sum() >= min_angle].index.astype('str').tolist()
    pricer_efficiencies_colors = df[df.efficiency / df.efficiency.sum() >= min_angle].colors.tolist()
    if df.efficiency[df.efficiency / df.efficiency.sum() < min_angle].sum() > 0:
        pricer_efficiencies.append(df.nVars[df.efficiency / df.efficiency.sum() < min_angle].sum() / df.time[df.efficiency / df.efficiency.sum() < min_angle].sum())
        pricer_efficiencies_labels.append('Others')
        pricer_efficiencies_colors.append('grey')

    print '    extracted time data:', time.time() - start_time
    start_time = time.time()

    # create the subplots
    fig = plt.gcf()
    gs = gridspec.GridSpec(2,2, wspace = .2)
    ax_total = plt.subplot(gs[0])
    ax_pricers = plt.subplot(gs[1])
    ax_nVars = plt.subplot(gs[2])
    ax_efficiencies = plt.subplot(gs[3])

    # plots
    numeric_labels = []

    numeric_labels += ax_total.pie([redcost_time, farkas_time, masterlp_time], labels = ['redcostpricing', 'initialfarkas', 'masterlp'], startangle = 180, counterclock = False, autopct = (lambda x: str(round(x * sum([redcost_time, farkas_time, masterlp_time]) / 10000.,2))), pctdistance = .75)[2]
    ax_total.axis('equal')

    if len(pricer_times) == 1:
        numeric_labels += ax_pricers.pie([t / sum(pricer_times) for t in pricer_times], labels = pricer_times_labels, colors = pricer_times_colors, startangle = 90, counterclock = False, autopct = str(pricer_times[0]), pctdistance = .75)[2]
    elif len(pricer_times) > 1:
        numeric_labels += ax_pricers.pie([t / sum(pricer_times) for t in pricer_times], labels = pricer_times_labels, colors = pricer_times_colors, startangle = 90, counterclock = False, autopct = (lambda x: str(round(x * sum(pricer_times),2))), pctdistance = .75)[2]
    ax_pricers.axis('equal')

    if len(pricer_nVars) == 1:
        numeric_labels += ax_nVars.pie([t / sum(pricer_nVars) for t in pricer_nVars], labels = pricer_nVars_labels, colors = pricer_nVars_colors, startangle = 90, counterclock = False, autopct = str(pricer_nVars[0]), pctdistance = .75)[2]
    elif len(pricer_nVars) > 1:
        numeric_labels += ax_nVars.pie(pricer_nVars, labels = pricer_nVars_labels, colors = pricer_nVars_colors, startangle = 90, counterclock = False, autopct = (lambda x: str(int(x * sum(pricer_nVars) / 100.))), pctdistance = .75)[2]
    ax_nVars.axis('equal')

    if len(pricer_efficiencies) == 1:
        numeric_labels += ax_efficiencies.pie([e / sum(pricer_efficiencies) for e in pricer_efficiencies], labels = pricer_efficiencies_labels, colors = pricer_efficiencies_colors, startangle = 90, counterclock = False, autopct = str(pricer_efficiencies[0]), pctdistance = .75)[2]
    elif len(pricer_efficiencies) > 1:
        numeric_labels += ax_efficiencies.pie([e / sum(pricer_efficiencies) for e in pricer_efficiencies], labels = pricer_efficiencies_labels, colors = pricer_efficiencies_colors, startangle = 90, counterclock = False, autopct = (lambda x: str(round(x * sum(pricer_efficiencies) / 100.,2))), pctdistance = .75)[2]
    ax_efficiencies.axis('equal')

    # format plots
    fig.suptitle('\\underline{\\textbf{' + info['instance'].replace('_','\_') + '}}', size = 'large')
    ax_total.set_title('Total Timeshares')
    ax_pricers.set_title('Timeshares of the Pricing Problems [s]')
    ax_nVars.set_title('\# of found Variables of the Pricing Problems')
    ax_efficiencies.set_title('Variables per second')
    for label in numeric_labels:
        label.set_color('white')

    print '    plotted time data:', time.time() - start_time
    start_time = time.time()

    # save
    fig.set_size_inches(11.7,8.3)
    filename = params['outdir'] + '/' + info['instance'] + '.times.' + info['settings']
    n = ''
    dot = ''
    if params['details']:
        while os.path.isfile(filename + dot + str(n) + '.pdf'):
            if n == '':
                n = 2
                dot = '.'
            else:
                n += 1
        plt.savefig(filename + dot + str(n) + '.pdf')
    else:
        while os.path.isfile(filename + dot + str(n) + '.png'):
            if n == '':
                n = 2
                dot = '.'
            else:
                n += 1
        plt.savefig(filename + dot + str(n) + '.png')
    plt.close()

    print '    saved time plot:', time.time() - start_time
    start_time = time.time()

def make_gap_plot(data, info, root_bounds):
    """
    For each instance create a plot, comparing the solving time of a prcing problem with the size of the gap in the root node
    :param data: collected pricing data as dataframe
    :param info: dictionary containing information about the data like the name of the instance, the settings & the scip_status
    :param root_bounds: collected root bounds data; necessary for the calculation of the gap at each round
    :return:
    """
    start_time = time.time()

    gap_data = data.query('(farkas == False) & (node == 1) & (pricing_prob >= 0)').groupby(['pricing_round','pricing_prob']).sum()
    rb = root_bounds.copy()
    rb.iter += gap_data.reset_index().pricing_round.min()
    rb = rb.rename(columns = {'iter': 'pricing_round'})[['pricing_round','pb','db','time']].set_index('pricing_round')
    gap_data = rb.join(gap_data, how = 'inner', lsuffix = '_rb')
    gap_data.db = gap_data.db.expanding().max()
    gap_data['gap'] = abs(gap_data.pb - gap_data.db)
    gap_data['gap'] = gap_data.gap / gap_data.gap.values[0]
    gap_data = gap_data.sort_values('time')

    print '    extracted gap data:', time.time() - start_time
    start_time = time.time()

    x = gap_data.time.values
    y = (1 - gap_data.gap).values
    mean_data = gap_data['gap'].groupby(gap_data.time.apply(lambda x: round(x,1))).mean()
    x_mean = mean_data.index.values
    y_mean = (1 - mean_data).values

    fig = plt.gcf()
    ax = fig.gca()

    ax.set_xlabel('Time of one pricing problem', size = 'large')
    ax.set_ylabel('Gap closed', size = 'large')
    ax.tick_params(axis='both', labelsize='large')

    ax.scatter(x, y, color = 'black')
    ax.plot(x_mean, y_mean, 'k--')

    print '    plotted time data:', time.time() - start_time
    start_time = time.time()

    fig.set_size_inches(11.7,8.3)
    filename = params['outdir'] + '/' + info['instance'] + '.gap.' + info['settings']
    n = ''
    dot = ''
    if params['details']:
        while os.path.isfile(filename + dot + str(n) + '.pdf'):
            if n == '':
                n = 2
                dot = '.'
            else:
                n += 1
        plt.savefig(filename + dot + str(n) + '.pdf')
    else:
        while os.path.isfile(filename + dot + str(n) + '.png'):
            if n == '':
                n = 2
                dot = '.'
            else:
                n += 1
        plt.savefig(filename + dot + str(n) + '.png')
    plt.close()

    print '    saved time plot:', time.time() - start_time
    start_time = time.time()

    return


def plots(data, info, root_bounds = None):
    """
    Master-function for plotting. Splits the data if necessary and calls all plotting functions (or a subset, according to the params)
    :param data: collected data as dataframe
    :param info: dictionary containing information about the data like the name of the instance, the settings & the scip_status
    :return:
    """

    # use tex to render the text output
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    # plot everything involving root bounds if possible and requested
    if not (root_bounds is None or params['no_gap']):
        make_gap_plot(data, info, root_bounds)

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
        # do not split the plots, but still check if rounds were neglected
        if params['minRound'] > 1 or params['maxRound'] > 0:
            # todo: include this min/max round information in the filenames (?)
            # todo: write a central method for saving the files
#            name += '_rounds' + str(params['minRound']) + 'to' + str(maxRnd)
            info['rounds_min'] = minRnd
            info['rounds_max'] = maxRnd
        data = data.query('node <= @maxNode & @minRnd <= pricing_round <= @maxRnd')
        if not params['no_bubble']:
            # build the bubble plot
            make_bubble_plot(data, info)
        if not params['no_summary']:
            # build the summary plot
            make_summary_plot(data, info)
        if not params['no_time']:
            # build the time plot
            make_time_plot(data, info)
        if not params['no_complete']:
            # do not build the complete plot
            make_plot(data, info)
    else:
        # split the plots by rounds
        fromRnd = minRnd - 1
        for i in range(1,(maxRnd-minRnd)+1):
            if i % params['nSplit'] <> 0 and i <> (maxRnd-minRnd):
                continue
            toRnd = i+minRnd
            data = data.query('@fromRnd < pricing_round <= @toRnd & node <= @maxNode')
#            name = name + '_rounds' + str(fromRnd + 1) + 'to' + str(toRnd)
            info['rounds_min'] = fromRnd + 1
            info['rounds_max'] = toRnd
            if not params['no_bubble']:
                # build the bubble plot
                make_bubble_plot(data, info)
            if not params['no_summary']:
                # build the summary plot
                make_summary_plot(data, info)
            if not params['no_time']:
                # build the time plot
                make_time_plot(data, info)
            if not params['no_complete']:
                # do not build the complete plot
                make_plot(data, info)
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
        filename, ext = os.path.splitext(os.path.basename(file))
        if not (ext == '.pkl'):
            print file, 'is not a pickle file'
            continue
        if params['instances'] <> '' and not any([(string in filename) for string in params['instances']]):
            print 'skipping', filename
            continue

        start_time = time.time()
        pkl_file = open(file, 'rb')
        objects = pickle.load(pkl_file)
        df = objects['pricing_data']
        info = objects['info']
        if 'root_bounds' in objects and not objects['root_bounds'].empty:
            root_bounds = objects['root_bounds']
        else:
            root_bounds = None
        pkl_file.close()

        print 'entering', info['instance']
        print '    loading data:', time.time() - start_time
        if df.empty:
            print '    no data found'
        else:
            start_time = time.time()
            plots(df, info, root_bounds)
            print '    total plotting:', time.time() - start_time
        print '    leaving', info['instance']

def collect_data(info, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, root_bounds = None):
    """
    Take lists containing the parsed data and structure them in a multiindex dataframe; then save or plot the data.
    All lists have to be of equal length (representing columns in a table)
    :param filename: name of the instance, settings & SCIP_status (format: instance.settings.scip_status)
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

    if not root_bounds is None and root_bounds.empty:
        root_bounds = None

    # save or plot the data
    if params['save']:
        start_time = time.time()
        output = open(params['outdir'] + '/' + info['instance'] + '.' + info['settings'] + '.pkl', 'wb')
        if root_bounds is None:
            pickle.dump({'pricing_data': df, 'info': info}, output, -1)
        else:
            pickle.dump({'pricing_data': df, 'info': info, 'root_bounds': root_bounds}, output, -1)
        print '    total saving:', time.time() - start_time
    else:
        start_time = time.time()
        plots(df, info, root_bounds)
        print '    total plotting:', time.time() - start_time

def parse_files(files):
    """
    Parse the (out-)files and structure the pricing-data in a dataframe
    :param files: List of files to be parsed
    :return:
    """
    for file in files:
        with open(file) as _file:
            first_line_of_file = True
            done = True

            for line in _file:
                if line.find("@0")<>-1 or first_line_of_file:
                    # if the file is an out-file, generated by the check-script, reset the variables whenever a new instance starts
                    if line.startswith("@01") or first_line_of_file:
                        # print message, if the previous problem is not done yet
                        if not done and problemFileName:
                            if not ind_node or not ind_pricing_round or not ind_stab_round or not ind_pricing_prob or not val_time or not val_nVars:
                                print '    no pricing data found'
                                done = True
                                continue
                            print '    ended abruptly'
                            collect_data({'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
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
                        root_bounds_ind = 0

                        # initialize flags
                        first_line_of_file = False
                        farkasDone = False
                        done = False
                        root_bounds = False
                        round_begin = False

                        # initialize all other variables
                        problemFileName = None
                        settings = 'default'
                        scip_status = 'NONE'
                        lptime_begin = 0
                        lptime_end = 0
                        root_bounds_data = pd.DataFrame()
                    else:
                        # ignore lines, where the output ends abrubtly (e.g. when the hard limit of the check-script is reached)
                        continue

                # if the problem is already processed, continue
                elif done:
                    continue

                elif line.startswith("loaded parameter file"):
                    # store current settings
                    settings = line.split()[-1]
                    settings = settings.split("/")[-1]
                    settings = os.path.splitext(settings)[0]

                elif not problemFileName and line.startswith("read problem "):
                    # get the problem name from the file name as in "check.awk"
                    tmparray = line.split("<")[-1].replace(">","").replace("\n","").split("/")[-1].split(".")
                    problemFileName = tmparray[0]
                    if tmparray[-1] == "gz" or tmparray[-1] == "z" or tmparray[-1] == "GZ" or tmparray[-1] == "Z":
                        tmparray.pop()
                    for i in range(1,len(tmparray)-1):
                        problemFileName += "." + tmparray[i]
                    if params['instances'] <> '' and not any([(string in problemFileName) for string in params['instances']]):
                        print 'skipping', problemFileName
                        done = True
                        continue
                    print 'entering', problemFileName

                # end of initial farkas pricing
                elif line.startswith("Starting reduced cost pricing..."):
                    farkasDone = True

                # read the SCIP status; end of pricing statistics
                elif line.startswith("SCIP Status        :"):
                    # continue if no data is found
                    if not ind_node or not ind_pricing_round or not ind_stab_round or not ind_pricing_prob or not val_time or not val_nVars:
                        print '    no pricing data found'
                        done = True
                        continue
                    scip_status = line.split(':')[1].split('[')[1].split(']')[0].replace(' ','_')
                    continue

                # read the root bounds table
                elif line.startswith("Root bounds"):
                    root_bounds = True

                elif root_bounds:
                    if line.startswith('Pricing Summary:'):
                        # end of the root bounds table; nothing more to read for this instance
                        root_bounds_data.iter = root_bounds_data.iter.astype('int')
                        collect_data({'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, root_bounds = root_bounds_data)
                        root_bounds = False
                        done = True
                        print '    leaving', problemFileName
                    elif line.startswith('iter	pb	db	time	dualdiff	dualoptdiff'):
                        # create a dataframe to store the table
                        root_bounds_data = pd.DataFrame(columns = line.split())
                        root_bounds_ind = 0
                    else:
                        root_bounds_data.loc[root_bounds_ind] = [float(s) for s in line.split()]
                        root_bounds_ind += 1

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
                        round_begin = True
                    except ValueError:
                        print '    ended abruptly'
                        collect_data({'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
                        print '    leaving', problemFileName
                        done = True
                        continue

                elif message.startswith("MLP t: "):
                    try:
                        if round_begin:
                            lptime_begin = float(message.split()[-1])
                            if len(ind_round) > 0 and lptime_begin - lptime_end >= 0.01:
                                # store all indices
                                ind_node.append(ind_node[-1])
                                ind_pricing_round.append(ind_pricing_round[-1])
                                ind_stab_round.append(ind_stab_round[-1])
                                ind_round.append(ind_round[-1])
                                # the Master LP Time is represented as a pricing problem with ID -2
                                ind_pricing_prob.append(-2)

                                # store the data
                                val_time.append(lptime_begin - lptime_end)
                                val_nVars.append(0)
                                val_farkas.append(val_farkas[-1])
                            round_begin = False
                        else:
                            lptime_end = float(message.split()[-1])
                            if lptime_end - lptime_begin > 0.005:
                                print 'It seems, that the LP time is not constant during a pricing round. Delta t is', lptime_end - lptime_begin
                    except ValueError:
                        print '    ended abruptly'
                        collect_data({'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
                        print '    leaving', problemFileName
                        done = True
                        continue

                elif message.startswith("Stabilization round ") or message.startswith("Sr "):
                    try:
                        stab_round = int(message.split()[-1])
                        round_counter += 1
                    except ValueError:
                        print '    ended abruptly'
                        collect_data({'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
                        print '    leaving', problemFileName
                        done = True
                        continue

                elif message.startswith("cp: ") or message.startswith("found "):
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
                        collect_data({'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
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
                        collect_data({'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)
                        print '    leaving', problemFileName
                        done = True
                        continue

            if not done:
                collect_data({'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas)

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
            logfile.close()

# Calling main script
if __name__ == '__main__':
    main()
