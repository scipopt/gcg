#!/usr/bin/env python3

import matplotlib
matplotlib.use('AGG')

import sys
import os
import argparse
import time
import datetime
from collections import Counter
from collections import OrderedDict
import gc

import pandas as pd
import pickle
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import lines
from matplotlib import transforms
from matplotlib import ticker
from matplotlib import patches as mpatches

sys.path.append("./")
import misc.vbc_reader as vbc

# Define the global parameter list
params = {}
debug = False
linelimit = 10**6 # modify depending on your system's (python-usable) RAM

# List of all implemented plots; first letter of the name will be the short option string of the command line argument
plotnames = ['complete',
             'summary',
             'bubble',
             'time',
             'gap',
             'depth',
             'nodeID',
             'vartimes']

def parse_arguments(args):
    """
    Parse the command-line arguments
    :param args: Command-line arguments provided during execution
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()

    # plottype arguments
    for plotname in plotnames:
        parser.add_argument('-' + plotname[0], '--' + plotname + '-only', action='store_true',
                            help='create only ' + plotname + ' plots')

        parser.add_argument('-' + plotname[0].upper(), '--no-' + plotname, action='store_true',
                            help='create no ' + plotname + ' plots')

    # arguments that specify which part of the data is to be stored/plotted
    parser.add_argument('--minround', type=int,
                        default=1,
                        help='start the data-collection or the plots with pricing-round MINROUND (default is the first round)')

    parser.add_argument('--maxround', type=int,
                        default=0,
                        help='end the data-collection or the plots with pricing-round MAXROUND (default is MAXROUND=0, which will be the last round)')

    parser.add_argument('--minnode', type=int,
                        default=1,
                        help='start the data-collection or the plots with node MINNODE (default is the root node)')

    parser.add_argument('--maxnode', type=int,
                        default=0,
                        help='end the data-collection or the plots with node MAXNODE (default is MAXNODE=0, which will be the last node)')

    parser.add_argument('--instances', type=str, nargs = '*',
                        default="",
                        help='names of the instances to be included in the plot/data-collection (can be only part of the name; default is all instances in FILENAMES)')

    # arguments that determine the look of the plots
    parser.add_argument('--splitrounds', type=int,
                        default=0,
                        help='the complete plot can be split into pieces containing a maximum of SPLITROUNDS rounds each (default is no splitting)')

    parser.add_argument('--colors', type=str,
                        default="nipy_spectral",
                        help='name of the color-map, that is used for the bars (see matplotlib documentation for maps, default is nipy_spectral)')

    parser.add_argument('--no-farkasline', action='store_true',
                        help='do not draw the blue line which marks the end of Farkas pricing')

    parser.add_argument('--lines', type=int, choices = list(range(0,3)),
                        default=1,
                        help='draw lines between pricing-rounds on the plots (0=never, 1=only for rounds that are not too short, 2=always)')

    parser.add_argument('--aggregate', action = 'store_true',
                        help='for each pricing round, draw only one aggregated bar in the complete plot')

    parser.add_argument('--short-times', action = 'store_true',
                        help='make slices down to really small ones, removing the absolute numbers in the time plot')

    parser.add_argument('--gapincomplete', action = 'store_true',
                        help='draw the gaps in the complete plot')

    parser.add_argument('--no-text', action='store_true',
                        help='do not write any text on the plots (such as node or round numbers)')

    parser.add_argument('--dualoptdiff', action = 'store_true',
                        help='for all plots involving the gap: use the best primal bound instead of the current one')

    parser.add_argument('--lptimeinsummary', action = 'store_true',
                        help='draw the Master LP time in the summary')

    parser.add_argument('--gapperround', action = 'store_true',
                        help='calculate the gap in the root node not for every pricing problem but for whole pricing iterations (affects only the gap plot)')

    # arguments concerning in- and output files
    parser.add_argument('-o', '--outdir', type=str,
                        default="plots",
                        help='output directory (default: "plots")')

    parser.add_argument('-save', '--savepickle', action='store_true',
                        help='saves the collected data in a pickle-file, in the OUTDIR, instead of plotting it (see also --load)')

    parser.add_argument('-load', '--loadpickle', action='store_true',
                        help='loads earlier collected data from a pickle-file instead of parsing a GCG-outfile and plots it (see also --save)')

    parser.add_argument('--png', action='store_true',
                        help='set this flag for a non-zoomable png-plot as output')

    parser.add_argument('--vbconly', action='store_true',
                        help='plot the plots involving only the vbc data (nodeID and depth); no pricing data is needed')

    parser.add_argument('--vbcdir', type=str,
                        default="../check/results/vbc",
                        help='directory of the vbc-files (needed for nodeID & depth plots; default: "../check/results/vbc")')

    parser.add_argument('filenames', nargs='+',
                        help='names of the files to be used for the plots; should be GCG output with STATISTICS=true, formatted as by the check-scripts for multiple instances or whole testsets')

    # parse the given command line arguments
    parsed_args = parser.parse_args(args)

    # check if the provided arguments are consistent
    if (0 < parsed_args.maxround and parsed_args.maxround < parsed_args.minround) or parsed_args.minround <= 0 or parsed_args.maxround < 0:
        print('please make sure that 1 <= MINROUND <= MAXROUND (or MAXROUND == 0 for the last possible round)')
        exit()
    if (0 < parsed_args.maxnode and parsed_args.maxnode < parsed_args.minnode) or parsed_args.minnode <= 0 or parsed_args.maxnode < 0:
        print('please make sure that 1 <= MINNODE <= MAXNODE (or MAXNODE == 0 for the last possible node)')
        exit()
    if not parsed_args.colors in plt.cm.datad:
        print('please use a colormap that is supported by pyplot (' + parsed_args.colors + ' is not supported)')
        exit()
    if parsed_args.loadpickle and parsed_args.savepickle:
        print('please load OR save data')
        exit()
    if not parsed_args.savepickle and all([do_not_draw for key, do_not_draw in vars(parsed_args).items() if key in ['no_' + plotname for plotname in plotnames]]):
        print('based on the passed parameters, no plot will be drawn')
        exit()

    return parsed_args

def set_params(args):
    """
    Set the global parameters from the parsed command-line arguments
    :param args: parsed command-line arguments
    :return:
    """
    # save the arguments in a global dict
    for plotname in plotnames:
        # evaluate which plots to plot
        params['no_' + plotname] = vars(args)['no_' + plotname] or any([vars(args)[other_plotname + '_only'] for other_plotname in plotnames if not other_plotname == plotname])
    for key, val in vars(args).items():
        # set all other params
        if not any([plotname == key.replace('no_','').replace('_only','') for plotname in plotnames]):
            params[key] = val
    return

def get_colmap(pricers, consider_masterlptime = True):
    """
    Returns a list of colors, with same length as pricers, that can be used for the bar plot
    Also returns a mapping color -> pricer_id, used for legends
    Each pricing problem has its own color
    :param pricers: a list with pricing_problem ids
    :param consider_masterlptime: reserve the first color of the colormap for the master lp time in case it is not included in the pricers list
    :return colors: a list of colors as used by pyplot.bar()
    :return mapping: a mapping pricer_id -> color
    """
    if len(pricers) == 0:
        print("   no pricing data found")
        return [(0.8,0.8,0.8,1)]
    # build a list of pricer ids, in which each pricer id appears once and sort it by id
    if consider_masterlptime:
        pricer_ids = [min(pricers) - 2]
    else:
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
        try: mapping[p] = cmap(pricer_ids.index(p))
        except ValueError:
            print("    unable to apply color scheme, using grey for all pricers.")
            for i in range(len(pricer_ids)):
                mapping[i] = (0.8,0.8,0.8,1)
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

def is_fin(x):
    """
    auxilliary function to determine if a value is finite
    """
    abs_inf = 10**20
    return abs_inf != abs(x)

def get_info_from_filename(filename):
    """
    Parses a filename and extracts the instance name, settings & scip_status
    :param filename: the name of a pickle-file without extension
    :return: info dictionary
    """
    info = {}
    info['status'] = filename.split('.')[-1]
    info['settings'] = filename.split('.')[-3]
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

def save_plot(fig, plotname, info):
    """
    Saves the plot
    :param fig: a matplotlib figure, that shall be saved
    :param plotname: name of the plot, e.g. bubble or summary
    :param info: dict containing information about the plot
    :return:
    """
    filename = params['outdir'] + '/' + info['instance'] + '.' + info['settings'] + '.' + plotname
    if 'rounds_min' in info:
        filename += '.fromRound' + info['rounds_min']
        if 'rounds_max' in info:
            filename += 'toRound' + info['rounds_max']
    elif 'rounds_max' in info:
        filename += '.toRound' + info['rounds_max']
    n = ''
    dot = ''
    if not params['png']:
        while os.path.isfile(filename + dot + str(n) + '.pdf'):
            if n == '':
                n = 2
                dot = '.'
            else:
                n += 1
        try:
            plt.savefig(filename + dot + str(n) + '.pdf')
            plt.close()
        except RuntimeError:
            print("Fatal: There was an error with your latex installation. Check that dvipng is installed.\nTerminating.")

    else:
        while os.path.isfile(filename + dot + str(n) + '.png'):
            if n == '':
                n = 2
                dot = '.'
            else:
                n += 1
        plt.savefig(filename + dot + str(n) + '.png')
        plt.close()
    return

def collect_gap_data(data, root_bounds):
    """
    Collect primal/dual bound and gap data per pricing round and problem, in preparation for gap and complete plots
    :param data: collected pricing data as dataframe
    :param root_bounds: collected root bounds data; necessary for the calculation of the gap at each round
    :return: The gap data
    """

    start_time = time.time()

    # calculate the gap data as combination of pricing statistics and root bounds data
    if params['gapperround']:

        gap_data = data.query('(farkas == False) & (node == 1) & (pricing_prob >= 0)').drop(['farkas', 'nVars'], axis = 1).groupby('pricing_round').sum()
    else:
        gap_data = data.query('(farkas == False) & (node == 1) & (pricing_prob >= 0)').drop(['farkas', 'nVars'], axis = 1).groupby(['pricing_round','pricing_prob']).sum()
    rb = root_bounds.copy()
    rb.iter += gap_data.reset_index().pricing_round.min()
    rb = rb.rename(columns = {'iter': 'pricing_round'})[['pricing_round','pb','db']].set_index('pricing_round')
    gap_data = rb.join(gap_data, how = 'inner', lsuffix = '_rb')
    #gap_data.db = gap_data.db.expanding().max()

    if params['dualoptdiff']:
        # check if the primal bound is the upper bound
        primal_is_upper = None
        if not gap_data.pb.dropna().empty and not gap_data.db.dropna().empty:
            if gap_data.pb.dropna().iloc[0] > gap_data.db.dropna().iloc[0]:
                primal_is_upper = True
            else:
                primal_is_upper = False

        # calculate the gap itself (normalized to 1 as the largest gap)
        if primal_is_upper is None:
            print('    cannot calculate dualoptdiff, since it is not clear if the primal bound is the upper or lower')
            gap_data['gap'] = abs(gap_data.pb - gap_data.db)
        elif primal_is_upper:
            gap_data.db = gap_data.db.cummax()
            gap_data['gap'] = gap_data.pb.min() - gap_data.db
        elif not primal_is_upper:
            gap_data.db = gap_data.db.cummin()
            gap_data['gap'] = gap_data.db - gap_data.pb.max()

    else:
        gap_data['gap'] = abs(gap_data.pb - gap_data.db)

    # normalize all finite gaps and set all infinite gaps to one
    max_gap = gap_data.gap.dropna()[is_fin(gap_data.pb) & is_fin(gap_data.db)].max()
    gap_data.gap = gap_data.gap / max_gap
    gap_data.loc[gap_data.gap > 1., 'gap'] = 1.
    gap_data.loc[gap_data.gap < 0., 'gap'] = 0.

    if debug: print('    extracted gap data:', time.time() - start_time)

    return gap_data

def make_complete_plot(data, info, gap_data, incumbent_times, rootlpsol_times):
    """
    Make a complete plot from the structured data
    :param data: dataframe with the collected data
    :param info: dictionary containing information about the data like the name of the instance, the settings & the scip_status
    :param gap data: collected gap data as dataframe
    :return:
    """
    try:
        start_time = time.time()

        # set the height of the zero bars
        ymin = -0.15

        # set some colors
        pricinground_linecolor = 'red'
        stabround_linecolor = 'green'
        farkas_linecolor = 'blue'
        gap_color = 'red'
        rootlpsol_color = 'blue'
        incumbent_color = 'orange'

        # flat out the data again
        data = data.reset_index()

        # workaround for most times being zero (0.01s is SCIPs smallest time-interval)
        # the column pool has no bar and therefore gets no time shift
        data.loc[data.pricing_prob != -1, 'time'] = data[data.pricing_prob != -1].time + 0.01

        # calculate the starting and ending time of each round
        data['starting_time'] = data.time.cumsum() - data.time
        data['ending_time'] = data.time.cumsum()

        if (not gap_data is None) and params['gapincomplete']:
            # calculate data for the gap plot
            if params['gapperround']:
                gap_rounds = [x for x in gap_data.index.values]
            else:
                gap_rounds = [x[0] for x in gap_data.index.values]
            x_gap = data[data['pricing_round'].isin(gap_rounds)].drop_duplicates('pricing_round', 'last').ending_time.values
            y_gap = gap_data.groupby('pricing_round').first()['gap'].values

        # set the height of the LP time bars to a maximum value
        data.loc[data.pricing_prob == -2, 'nVars'] = data.nVars.max() * 10

        # extract the column pool data and delete it from data
        x_colpool = data[data.pricing_prob == -1].starting_time.values
        y_colpool = data[data.pricing_prob == -1].nVars.values
        data = data[data.pricing_prob != -1].reset_index()

        # define position, width and height of the bars; the first two are defined by time, the last by nVars
        x = data.starting_time.values
        y = (data.nVars - ymin).values
        widths = data.time.values
        colors, cmapping = get_colmap(data.pricing_prob.values, consider_masterlptime = False)

        # define positions for incumbent and root lp solution plots
        rootlpsol_times_cnt = Counter(rootlpsol_times)
        incumbent_times_cnt = Counter(incumbent_times)
        incumbent_times_bottoms = list()
        for t in list(incumbent_times_cnt.keys()):
            if t in list(rootlpsol_times_cnt.keys()):
                incumbent_times_bottoms.append(-rootlpsol_times_cnt[t])
            else:
                incumbent_times_bottoms.append(ymin)
                incumbent_times_cnt[t] += ymin
        if len(incumbent_times) > 0:
            ymin_ncols = min([s-t for (s,t) in zip(incumbent_times_bottoms, list(incumbent_times_cnt.values()))])
        else:
            ymin_ncols = ymin

        # sometimes we need just the height of the bars of the pricing problems
        y_pricers = (data[data.pricing_prob >= 0].nVars - ymin).values

        if debug: print('    data restructured:', time.time() - start_time)
        start_time = time.time()

        # make the bar plot
        fig = plt.gcf()
        ax = plt.gca()
        if not params['png']:
            lw = 0.01
        else:
            lw = 1.0
        try:
            ax.bar(x, y, widths, bottom = ymin, align = 'edge', linewidth = lw, edgecolor = 'white', color = colors, label='Pricing problems')
        except MemoryError:
            f = plt.figure()
            f.clear()
            plt.close(f)
            print("Warning: Insufficient memory to make bars in complete plot for instance " + info['instance'] + ". Skipping.")
            return

        # add the column pool data as a scatter plot
        cp_scatter = ax.scatter(x_colpool, y_colpool, color = 'green', marker = 'o', s = 200, zorder = 10, label = 'Column Pool')
        try:
            ax.bar(list(rootlpsol_times_cnt.keys()), width=lw, height=[-t-ymin for t in list(rootlpsol_times_cnt.values())], bottom = ymin, align = 'edge', color = rootlpsol_color, label='Root LP Sol')
            ax.bar(list(incumbent_times_cnt.keys()), width=lw, height=[-t for t in list(incumbent_times_cnt.values())], bottom = incumbent_times_bottoms, align = 'edge', color = incumbent_color, label='Incumbent')
        except MemoryError:
            f = plt.figure()
            f.clear()
            plt.close(f)
            print("Warning: Insufficient memory to make bars in complete plot for instance " + info['instance'] + ". Skipping.")
            return

        if (not gap_data is None) and params['gapincomplete']:
            # add the gap plot
            ax2 = ax.twinx()
            gap_plot = ax2.plot(x_gap, y_gap, 'k--', color = gap_color, linewidth = 10.0, label = 'Gap')

        if debug: print('    data plotted:', time.time() - start_time)
        start_time = time.time()

        # set parameters
        if not params['png']:
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
        ax.set_ylim([ymin_ncols,ymax])
        ax.set_xlim([xmin,xmax])
        ax.get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True, nbins = 15))
        ax.tick_params(axis = 'both', length = textsize/2, width = textsize/40, labelsize = textsize*0.9, pad = 15)
        ax.set_xlabel('Time / s', size = 1.15*textsize)
        ax.set_ylabel('\# of variables', size = 1.15*textsize)
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        if (not gap_data is None) and params['gapincomplete']:
            ax2.set_ylim([ymin_ncols/ymax-0.01, 1.01])
            ax2.tick_params(axis = 'both', length = textsize/2, width = textsize/40, labelsize = textsize*0.9, pad = 15)
            ax2.set_ylabel('Gap', size = 1.15*textsize)

        if debug: print('    data formatted:', time.time() - start_time)
        start_time = time.time()

        if params['no_text']:
            # save the figure
            save_plot(fig, 'complete', info)
            plt.close()

            if debug: print('    save:', time.time() - start_time)
            if debug: print('    saved figure')

            return

        # special cases: no or only (initial) farkas pricing in the plot
        if data.farkas.all():
            ax.text(.991, .99, "\it{Initial Farkas Pricing did not end}", va = 'top', ha = 'right', rotation = 0, color = farkas_linecolor, zorder = 11, size = textsize * .95, transform = ax.transAxes, bbox=dict(facecolor = 'white', edgecolor = 'none', alpha = .85, pad = 20))
            farkasLine = True
        elif not data.farkas.any():
            ax.text(.009, .99, "\it{No initial Farkas Pricing}", va = 'top', ha = 'left', rotation = 0, color = farkas_linecolor, zorder = 11, size = textsize * .95, transform = ax.transAxes, bbox=dict(facecolor = 'white', edgecolor = 'none', alpha = .85, pad = 20))
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
                if params['lines'] == 2 or (params['lines'] == 1 and (pos - prev_x_drawn)/totalTime > 0.002) or (not params['no_farkasline'] and not farkasLine and not far):
                    line = lines.Line2D([pos,pos],[0,1],color=pricinground_linecolor,linewidth=1.0, transform = trans)
                    # blue line at the end of farkas pricing
                    if not farkasLine and not far:
                        line.set_color(farkas_linecolor)
                    ax.add_line(line)
                    prev_x_drawn = pos
                # text for initial Farkas pricing
                if not farkasLine and not far:
                    if pos <= (xmax + xmin) / 2:
                        align = 'left'
                    else:
                        align = 'right'
                    ax.text(pos, .99, "\it{End of initial Farkas Pricing}", va = 'top', ha = align, rotation = 0, color = farkas_linecolor, zorder = 11, size = textsize * .95, transform = trans, bbox=dict(facecolor = 'white', edgecolor = 'none', alpha = .85, pad = 20))
                    farkasLine = True
                # write the round number, if there is space for it
                if len(texts) == 0 or get_x1_in_data(texts[-1], fig) < prev_x:
                    texts.append(ax.text(prev_x, 1.01, str(prev_rnd), rotation='vertical',va='bottom', ha='left', size = textsize, transform = trans))
                prev_rnd = rnd
                prev_x = pos
            else:
                # dashed line for a new stabilization round
                if params['lines'] == 2 or (params['lines'] == 1 and (pos - prev_x_drawn)/totalTime > 0.0005):
                    line = lines.Line2D([pos,pos],[0,1],color=stabround_linecolor,linestyle='--',linewidth=0.8, transform = trans)
                    ax.add_line(line)
                    prev_x_drawn = pos
        if len(texts) == 0 or get_x1_in_data(texts[-1], fig) < prev_x:
            texts.append(ax.text(prev_x, 1.01, str(prev_rnd), rotation='vertical',va='bottom', ha='left', size = textsize, transform = trans))
        text_height = max(get_y1_in_ax(texts[0], fig),get_y1_in_ax(texts[-1], fig))

        if debug: print('    stab- and pricing-round information:', time.time() - start_time)
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

        if debug: print('    node information:', time.time() - start_time)
        start_time = time.time()

        # draw a legend, but do not include more than 25 pricing problems
        if -2 in data.pricing_prob.values:
            patches = [mpatches.Patch(color = cmapping[p], label = 'Pricing Problem ' + str(p)) for p in cmapping]
            patches[0].set_label('Master LP Time')
            if params['aggregate']:
                patches[1].set_label('Pricing Problems')
        else:
            patches = [mpatches.Patch(color = cmapping[p], label = 'Pricing Problem ' + str(p)) for p in cmapping]
            if params['aggregate']:
                patches[0].set_label('Pricing Problems')
        if len(patches) > 31:
            patches = patches[:31] + [mpatches.Patch(color = 'white', alpha = 0, label = '...')]
        handles = patches + [cp_scatter]
        if params['lines'] > 0:
            handles += [lines.Line2D([0,0], [0,1], color = pricinground_linecolor, linewidth = 12., label = 'Pricing Round'), lines.Line2D([0,0], [0,1], color = stabround_linecolor, linestyle = '--', linewidth = 1.6, label = 'Mis-price Iteration')]
        if params['gapincomplete']:
            handles += [lines.Line2D([0,0], [0,1], color = gap_color, linestyle = '--', linewidth = 12., label = 'Gap')]
        handles += [lines.Line2D([0,0], [0,1], color = rootlpsol_color, linewidth = 16., label = 'Root LP Sol'), lines.Line2D([0,0], [0,1], color = incumbent_color, linewidth = 16., label = 'Incumbent')]
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
        save_plot(fig, 'complete', info)
        plt.close()

        if debug: print('    save:', time.time() - start_time)
        if debug: print('    saved figure')

        return
    except MemoryError:
        f = plt.figure()
        f.clear()
        plt.close(f)
        print("Warning: Insufficient memory to make complete plot for instance " + info['instance'] + ". Skipping.")
        gc.collect()
        return

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
    summary = pd.DataFrame()
    summary['time'] = data.query('pricing_prob != -2').groupby(level=['node','pricing_round','stab_round', 'round']).sum().time
    summary['mlp_time'] = data.query('pricing_prob == -2').groupby(level=['node','pricing_round','stab_round', 'round']).sum().time
    summary['mlp_time'] = summary['mlp_time'].fillna(0.00)
    summary['found_frac'] = data.query('pricing_prob != -2').astype(bool).groupby(level=['node','pricing_round','stab_round', 'round']).sum().nVars/data.query('pricing_prob != -2').groupby(level=['node','pricing_round','stab_round', 'round']).count().nVars*100
    summary = summary.reset_index()

    if not data.farkas.all() and data.farkas.any():
        # get the last round of initial farkas pricing
        farkas_end = (data[data.farkas == False].reset_index()['round'].values[0] + data[data.farkas == True].reset_index()['round'].values[-1])/2.

    if debug: print('    extracted summary data:', time.time() - start_time)
    start_time = time.time()

    fig,ax1 = plt.subplots()
    ax2 = ax1.twinx()

    # get the data for the plot
    x = summary[summary.stab_round <= 0]['round'].values
    y_time = summary[summary.stab_round <= 0].time.values
    y_mlp_time = summary[summary.stab_round <= 0].mlp_time.values
    y_found_frac = summary[summary.stab_round <= 0].found_frac.values
    x_stab = summary[summary.stab_round > 0]['round'].values
    y_stab_time = summary[summary.stab_round > 0].time.values
    y_stab_mlp_time = summary[summary.stab_round > 0].mlp_time.values
    y_stab_found_frac = summary[summary.stab_round > 0].found_frac.values
    x_mean = summary['round'].values
    y_mean_time = summary.time.rolling(int(.05*(max(x_mean) - min(x_mean))), center = True).mean()
    y_mean_mlp_time = summary.mlp_time.rolling(int(.05*(max(x_mean) - min(x_mean))), center = True).mean()
    y_mean_found_frac = summary.found_frac.rolling(int(.05*(max(x_mean) - min(x_mean))), center = True).mean()

    # format the plot
    ax1.set_xlabel('Pricing Round', size='large')
    ax1.set_ylabel('Time / s', color='k', size='large')
    ax2.set_ylabel('Fraction of successful pricing problems / \%', color='r', size='large')

    ax1.tick_params(axis='both', labelsize='large')
    ax2.tick_params(axis='both', labelsize='large')

    # set the axes limits
    xmin = 0
    xmax = max(x.tolist() + x_stab.tolist()) + 0.9
    ax1.set_xlim([xmin, xmax])
    if max(y_time) > 0 or max(y_mlp_time):
        ax1.set_ylim([-max(max(y_time), max(y_mlp_time)) * 0.1, max(max(y_time), max(y_mlp_time)) * 1.1])
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
    handles = []
    handles.append(ax1.scatter(x,y_time, color='k', s=perimeter**2, zorder = 2, label = 'Pricing Time'))
    handles.append(ax2.scatter(x,y_found_frac, color='r', s=perimeter**2, label = 'Success (column generated)'))
    ax1.scatter(x_stab,y_stab_time, color='k', s=perimeter**2, marker='x', alpha=.5, zorder = 2, label = 'Pricing Time in Stabilization Round')
    ax2.scatter(x_stab,y_stab_found_frac, color='r', s=perimeter**2, marker='x', alpha=.5, label = 'Success in Stabilization Round')
    ax1.plot(x_mean,y_mean_time, 'k--', zorder = 2, label = None)
    ax2.plot(x_mean,y_mean_found_frac, 'r--', label = None)
    if params['lptimeinsummary']:
        handles.append(ax1.scatter(x,y_mlp_time, color='g', s=perimeter**2, zorder = 1, label = 'Master LP Time'))
        ax1.scatter(x_stab,y_stab_mlp_time, color='g', s=perimeter**2, marker='x', alpha=.5, zorder = 1, label = 'Master LP Time in Stabilization Rounds')
        ax1.plot(x_mean,y_mean_mlp_time, 'g--', zorder = 1, label = None)

    # add a line after the root-node
    if summary.node.max() > 1:
        try:
            x_line = (summary[summary.node > 1].index[0] + 1 + summary[summary.node == 1].index[-1] + 1)/2.
        except:
            print("    Information: Could not add line after root-node for {}".format(info["instance"]))
            x_line = 0
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
        ax1.text(1., 1.025, "\it{Ended within Root}", ha = 'right', size = 'smaller', color = 'red', zorder = 1, transform = ax1.transAxes)

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
    handles.append(lines.Line2D([0,0], [0,0], color = 'None', marker = 'x', markerfacecolor = 'k', markeredgecolor = 'k', markersize = 5, alpha = .7, label = '... in Stabilization Round'))
    plt.legend(handles = handles, loc = 3, bbox_to_anchor = (.0, 1.04, 1., 1.04), ncol = 4, mode = 'expand')

    # add other information
    ax1.text(.3, 1.11, '\\textbf{\\underline{' + info['instance'].replace('_','\_') + '}}', ha = 'center', size = 'large', transform = ax1.transAxes)
    ax1.text(.75, 1.11, '\\textbf{Settings:} \\textit{' + info['settings'].replace('_','\_') + '}', ha = 'right', size = 'medium', transform = ax1.transAxes)
    ax1.text(1., 1.11, '\\textbf{SCIP Status:} \\textit{' + info['status'].replace('_',' ') + '}', ha = 'right', size = 'medium', transform = ax1.transAxes)

    if debug: print('    plotted summary:', time.time() - start_time)
    start_time = time.time()

    # save the plot
    fig.set_size_inches(11.7,8.3)
    plt.tight_layout()
    fig.subplots_adjust(top = 0.88)
    save_plot(fig, 'summary', info)
    plt.close()
    if debug: print('    saved summary:', time.time() - start_time)

def make_bubble_plot(data, info):
    """
    For each instance create one bubble plot, that shows which pricing problem found a variable in each round
    :param data: dataframe with the collected (complete) data
    :param info: dictionary containing information about the data like the name of the instance, the settings & the scip_status
    :return:
    """
    start_time = time.time()

    # flat out the data again and do not consider the master lp time in this plot
    data = data.query('pricing_prob != -2').reset_index()

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
    cmapping = get_colmap(data[data.pricing_prob != -1].pricing_prob.unique())[1]
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

    if debug: print('    extracted bubble data:', time.time() - start_time)
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
        try:
            x_line = (data[data.node > 1]['round'].iloc[0] + data[data.node == 1]['round'].iloc[-1])/2.
        except:
            print("    Information: Could not add line after root-node for {}".format(info["instance"]))
            x_line = 0
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
        ax.text(1., 1.025, "\it{Ended within Root}", ha = 'right', size = 'smaller', color = 'red', zorder = 11, transform = ax.transAxes)

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

    if debug: print('    plotted bubble data:', time.time() - start_time)
    start_time = time.time()

    # draw a legend
    handles = []
    handles.append(lines.Line2D([0,0], [0,1], color = 'None', marker = 'o', markerfacecolor = 'k', markersize = 5, label = 'Pricer has found at least one variable'))
    handles.append(lines.Line2D([0,0], [0,1], color = 'None', marker = 'o', markerfacecolor = 'green', markersize = 5, label = 'Variables were taken from column pool (ID $-1$)'))
    handles.append(lines.Line2D([0,0], [0,1], color = 'None', marker = 'x', markerfacecolor = 'k', markeredgecolor = 'k', markersize = 5, alpha = .5, label = 'Pricer has found at least one variable in stab. round'))
    ax.legend(handles = handles, loc = 3, bbox_to_anchor = (.0, 1.04, 1.18, .02), ncol = 3)

    # add other information
    trans = transforms.blended_transform_factory(fig.transFigure, ax.transAxes)
    ax.text(.3, 1.11, '\\textbf{\\underline{' + info['instance'].replace('_','\_') + '}}', ha = 'center', size = 'large', transform = trans)
    ax.text(.75, 1.11, '\\textbf{Settings:} \\textit{' + info['settings'].replace('_','\_') + '}', ha = 'right', size = 'medium', transform = trans)
    ax.text(.975, 1.11, '\\textbf{SCIP Status:} \\textit{' + info['status'].replace('_',' ') + '}', ha = 'right', size = 'medium', transform = trans)

    # save the plot
    fig.set_size_inches(11.7,8.3)
    gs.tight_layout(fig,rect = (0,0,1,.9))
    save_plot(fig, 'bubble', info)
    plt.close()

    if debug: print('    saved bubble plot:', time.time() - start_time)

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

    # change this to show more pricing problems ###
    short_times = params['short_times']
    if short_times:
        min_angle = 4./360.
        deactivateTimeLabels = True
    else:
        min_angle = 11./360.
        deactivateTimeLabels = False

    # calculate times for the pricer summary
    df = data.query('pricing_prob >= 0').reset_index()[['pricing_prob','time','nVars']].groupby('pricing_prob').sum()
    df['colors'] = get_colmap(df.index.tolist())[0]
    df = df[df.time >= 0.01].sort_values('time', ascending = False)
    pricer_times = df.time[df.time / df.time.sum() >= min_angle].tolist()
    pricer_times_labels = df[df.time / df.time.sum() >= min_angle].index.astype('str').tolist()
    pricer_times_colors = df.colors[df.time / df.time.sum() >= min_angle].tolist()
    if df.time[df.time / df.time.sum() < min_angle].sum() > 0:
        pricer_times.append(df.time[df.time / df.time.sum() < min_angle].sum())
        pricer_times_labels.append('Others ({})'.format(str(len(df.time)-len(pricer_times)+1)))
        pricer_times_colors.append('grey')

    # calculate number of found variables
    pricer_nVars = df.nVars[df.nVars / df.nVars.sum() >= min_angle].tolist()
    pricer_nVars_labels = df[df.nVars / df.nVars.sum() >= min_angle].index.astype('str').tolist()
    pricer_nVars_colors = df.colors[df.nVars / df.nVars.sum() >= min_angle].tolist()
    if df.nVars[df.nVars / df.nVars.sum() < min_angle].sum() > 0:
        pricer_nVars.append(df.nVars[df.nVars / df.nVars.sum() < min_angle].sum())
        pricer_nVars_labels.append('Others ({})'.format(str(len(df.time)-len(pricer_nVars)+1)))
        pricer_nVars_colors.append('grey')

    # calculate efficiencies (var/time)
    df['efficiency'] = df.nVars / df.time
    pricer_efficiencies = df.efficiency[df.efficiency / df.efficiency.sum() >= min_angle].tolist()
    pricer_efficiencies_labels = df.efficiency[df.efficiency / df.efficiency.sum() >= min_angle].index.astype('str').tolist()
    pricer_efficiencies_colors = df[df.efficiency / df.efficiency.sum() >= min_angle].colors.tolist()
    if df.efficiency[df.efficiency / df.efficiency.sum() < min_angle].sum() > 0:
        #pricer_efficiencies.append(df.nVars[df.efficiency / df.efficiency.sum() < min_angle].sum() / df.time[df.efficiency / df.efficiency.sum() < min_angle].sum())
        pricer_efficiencies.append(sum(df['efficiency']) - sum(pricer_efficiencies))
        pricer_efficiencies_labels.append('Others ({})'.format(str(len(df.time)-len(pricer_efficiencies)+1)))
        pricer_efficiencies_colors.append('grey')

    # calculate efficiencies2 (time/var)
    df['efficiency2'] =  df.time / df.nVars
    df.efficiency2 = df.efficiency2.replace([np.inf, -np.inf], 0)
    pricer_efficiencies2 = df.efficiency2[df.efficiency2 / df.efficiency2.sum() >= min_angle].tolist()
    pricer_efficiencies2_labels = df.efficiency2[df.efficiency2 / df.efficiency2.sum() >= min_angle].index.astype('str').tolist()
    pricer_efficiencies2_colors = df[df.efficiency2 / df.efficiency2.sum() >= min_angle].colors.tolist()
    if df.efficiency2[df.efficiency2 / df.efficiency2.sum() < min_angle].sum() > 0:
        #pricer_efficiencies2.append(df.nVars[df.efficiency2 / df.efficiency2.sum() < min_angle].sum() / df.time[df.efficiency2 / df.efficiency2.sum() < min_angle].sum())
        pricer_efficiencies2.append(sum(df['efficiency2']) - sum(pricer_efficiencies2))
        pricer_efficiencies2_labels.append('Others ({})'.format(str(len(df.time)-len(pricer_efficiencies2)+1)))
        pricer_efficiencies2_colors.append('grey')

    if debug: print('    extracted time data:', time.time() - start_time)
    start_time = time.time()

    # create the subplots
    fig = plt.gcf()
    eff = True
    if eff:
        gs = gridspec.GridSpec(2,3, wspace = .2)
        ax_total = plt.subplot(gs[:,0])
        ax_pricers = plt.subplot(gs[0,1])
        ax_nVars = plt.subplot(gs[0,-1])
        ax_efficiencies = plt.subplot(gs[1,1])
        ax_efficiencies2 = plt.subplot(gs[1,-1])
    else:
        gs = gridspec.GridSpec(2,2, wspace = .2)
        ax_total = plt.subplot(gs[0])
        ax_pricers = plt.subplot(gs[1])
        ax_nVars = plt.subplot(gs[2])
        ax_efficiencies = plt.subplot(gs[3])

    plt.figtext(.01,.01,'Number of pricing problems: ' + str(len(df.time)) + '.',size='medium')

    # plots
    numeric_labels = []

    numeric_labels += ax_total.pie([redcost_time, farkas_time, masterlp_time], labels = ['redcostpricing', 'initialfarkas', 'masterlp'], startangle = 180, counterclock = False, autopct = (lambda x: str(round(x * sum([redcost_time, farkas_time, masterlp_time]) / 10000.,2))), pctdistance = .75)[2]
    ax_total.axis('equal')

    if len(pricer_times) == 1:
        if deactivateTimeLabels:
            numeric_labels += ax_pricers.pie([t / sum(pricer_times) for t in pricer_times], labels = pricer_times_labels, colors = pricer_times_colors, startangle = 90, counterclock = False)[1]
        else:
            numeric_labels += ax_pricers.pie([t / sum(pricer_times) for t in pricer_times], labels = pricer_times_labels, colors = pricer_times_colors, startangle = 90, counterclock = False, autopct = str(pricer_times[0]), pctdistance = .75)[2]
    elif len(pricer_times) > 1:
        if deactivateTimeLabels:
            numeric_labels += ax_pricers.pie([t / sum(pricer_times) for t in pricer_times], labels = pricer_times_labels, colors = pricer_times_colors, startangle = 90, counterclock = False)[1]
        else:
            numeric_labels += ax_pricers.pie([t / sum(pricer_times) for t in pricer_times], labels = pricer_times_labels, colors = pricer_times_colors, startangle = 90, counterclock = False, autopct = (lambda x: str(round(x * sum(pricer_times) / 100.,2))), pctdistance = .75)[2]
    ax_pricers.axis('equal')

    if len(pricer_nVars) == 1:
        if deactivateTimeLabels:
            numeric_labels += ax_nVars.pie([t / sum(pricer_nVars) for t in pricer_nVars], labels = pricer_nVars_labels, colors = pricer_nVars_colors, startangle = 90, counterclock = False)[1]
        else:
            numeric_labels += ax_nVars.pie([t / sum(pricer_nVars) for t in pricer_nVars], labels = pricer_nVars_labels, colors = pricer_nVars_colors, startangle = 90, counterclock = False, autopct = str(pricer_nVars[0]), pctdistance = .75)[2]
    elif len(pricer_nVars) > 1:
        if deactivateTimeLabels:
            numeric_labels += ax_nVars.pie(pricer_nVars, labels = pricer_nVars_labels, colors = pricer_nVars_colors, startangle = 90, counterclock = False)[1]
        else:
            numeric_labels += ax_nVars.pie(pricer_nVars, labels = pricer_nVars_labels, colors = pricer_nVars_colors, startangle = 90, counterclock = False, autopct = (lambda x: str(int(round(x * sum(pricer_nVars) / 100.)))), pctdistance = .75)[2]
    ax_nVars.axis('equal')

    if len(pricer_efficiencies) == 1:
        if deactivateTimeLabels:
            numeric_labels += ax_efficiencies.pie([e / sum(pricer_efficiencies) for e in pricer_efficiencies], labels = pricer_efficiencies_labels, colors = pricer_efficiencies_colors, startangle = 90, counterclock = False)[1]
        else:
            numeric_labels += ax_efficiencies.pie([e / sum(pricer_efficiencies) for e in pricer_efficiencies], labels = pricer_efficiencies_labels, colors = pricer_efficiencies_colors, startangle = 90, counterclock = False, autopct = str(pricer_efficiencies[0]), pctdistance = .75)[2]
    elif len(pricer_efficiencies) > 1:
        if deactivateTimeLabels:
            numeric_labels += ax_efficiencies.pie([e / sum(pricer_efficiencies) for e in pricer_efficiencies], labels = pricer_efficiencies_labels, colors = pricer_efficiencies_colors, startangle = 90, counterclock = False)[1]
        else:
            numeric_labels += ax_efficiencies.pie([e / sum(pricer_efficiencies) for e in pricer_efficiencies], labels = pricer_efficiencies_labels, colors = pricer_efficiencies_colors, startangle = 90, counterclock = False, autopct = (lambda x: str(round(x * sum(pricer_efficiencies) / 100.,2))), pctdistance = .75)[2]
    ax_efficiencies.axis('equal')

    if len(pricer_efficiencies2) == 1:
        if deactivateTimeLabels:
            numeric_labels += ax_efficiencies2.pie([e / sum(pricer_efficiencies2) for e in pricer_efficiencies2], labels = pricer_efficiencies2_labels, colors = pricer_efficiencies2_colors, startangle = 90, counterclock = False)[1]
        else:
            numeric_labels += ax_efficiencies2.pie([e / sum(pricer_efficiencies2) for e in pricer_efficiencies2], labels = pricer_efficiencies2_labels, colors = pricer_efficiencies2_colors, startangle = 90, counterclock = False, autopct = str(pricer_efficiencies2[0]), pctdistance = .75)[2]
    elif len(pricer_efficiencies2) > 1:
        if deactivateTimeLabels:
            numeric_labels += ax_efficiencies2.pie([e / sum(pricer_efficiencies2) for e in pricer_efficiencies2], labels = pricer_efficiencies2_labels, colors = pricer_efficiencies2_colors, startangle = 90, counterclock = False)[1]
        else:
            numeric_labels += ax_efficiencies2.pie([e / sum(pricer_efficiencies2) for e in pricer_efficiencies2], labels = pricer_efficiencies2_labels, colors = pricer_efficiencies2_colors, startangle = 90, counterclock = False, autopct = (lambda x: str(round(x * sum(pricer_efficiencies2) / 100.,2))), pctdistance = .75)[2]
    ax_efficiencies2.axis('equal')

    # format plots
    fig.suptitle('\\underline{\\textbf{' + info['instance'].replace('_','\_') + '}} \\qquad \\textbf{Settings:} ' + info['settings'].replace('_','\_') + '\\quad \\textbf{SCIP Status:} ' + info['status'].replace('_',' '), size = 'large')
    ax_total.set_title('Total Timeshares')
    ax_pricers.set_title('Timeshares of the Pricing Problems [s]')
    ax_nVars.set_title('\# of found Variables of the Pricing Problems')
    ax_efficiencies.set_title('Variables per second')
    ax_efficiencies2.set_title('Seconds per variable')
    for label in numeric_labels:
        if deactivateTimeLabels:
            label.set_color('black')
        else:
            label.set_color('white')

    if debug: print('    plotted time data:', time.time() - start_time)
    start_time = time.time()

    # save
    fig.set_size_inches(11.7,8.3)

    if debug: print('    set plot size:', time.time() - start_time)
    start_time = time.time()

    save_plot(fig, 'time', info)

    if debug: print('    saved time plot:', time.time() - start_time)
    start_time = time.time()


    plt.close()

    if debug: print('    close time plot:', time.time() - start_time)
    start_time = time.time()

def make_gap_plot(gap_data, info):
    """
    For each instance create a plot, comparing the solving time of a pricing problem with the size of the gap in the root node
    :param gap data: collected gap data as dataframe
    :param info: dictionary containing information about the data like the name of the instance, the settings & the scip_status
    :return:
    """
    start_time = time.time()
    time_gap_data = gap_data.sort_values('gap')
    mean_data = time_gap_data['time'].groupby(time_gap_data.gap.apply(lambda x: round(x,2))).mean()
    #mean_data = mean_data.groupby(mean_data.apply(lambda x: round(x,2))).mean()
    if len(mean_data) >= 20:
        mean_data = mean_data.rolling(max(int(.08 * len(mean_data)), 5), center = True).mean()

    if debug: print('    sorted gap data:', time.time() - start_time)
    start_time = time.time()

    # set x and y values
    y = time_gap_data.time.values
    x = (1 - time_gap_data.gap).values
    y_mean = (mean_data).values
    x_mean = (1 - mean_data.index).values

    # create the figure
    fig = plt.gcf()
    ax = fig.gca()

    # format the plot
    ax.text(.0, 1.03, '\\textbf{\\underline{' + info['instance'].replace('_','\_') + '}}', ha = 'left', size = 'large', transform = ax.transAxes)
    description = 'Gap between primalbound and'
    if params['dualoptdiff']:
        description +=' \\textit{best}'
    else:
        description += ' \\textit{current}'
    description += ' dual bound \n vs Duration of pricing'
    if params['gapperround']:
        description += ' \\textit{iteration}'
    else:
        description += ' \\textit{problem}'
    ax.text(.24, 1.03, 'In the root node:', ha = 'right', size = 'small', transform = ax.transAxes)
    ax.text(.25, 1.0325, description, ha = 'left', va = 'center', size = 'medium', transform = ax.transAxes)
    ax.text(.75, 1.03, '\\textbf{Settings:} \\textit{' + info['settings'].replace('_','\_') + '}', ha = 'right', size = 'medium', transform = ax.transAxes)
    ax.text(1., 1.03, '\\textbf{SCIP Status:} \\textit{' + info['status'].replace('_',' ') + '}', ha = 'right', size = 'medium', transform = ax.transAxes)
    if params['gapperround']:
        ax.set_ylabel('Time of one pricing round', size = 'large')
        color = 'k'
    else:
        ax.set_ylabel('Time of one pricing problem (resolution: $0.01$s)', size = 'large')
        color = get_colmap(time_gap_data.index.get_level_values('pricing_prob').tolist())[0]
    ax.set_xlabel('Gap closed', size = 'large')
    ax.set_xlim([-0.04,1.04])
    ax.tick_params(axis='both', labelsize='large')
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base = .1))

    # plot the data
    ax.scatter(x, y, color = color, label='a pricing iteration that needed $y$ amount of time\nwhen $x$ gap was already closed')
    ax.plot(x_mean, y_mean, 'k--', label='mean')
    ax.legend(loc='best')

    if debug: print('    plotted gap data:', time.time() - start_time)
    start_time = time.time()

    # save the plot
    fig.set_size_inches(11.7,8.3)
    try:
        fig.tight_layout()
    except RuntimeError:
        print("Fatal: There was an error with your latex installation. Check that type1cm (included in texlive-latex-extra) is installed.\nTerminating.")
    fig.subplots_adjust(top = .925)
    save_plot(fig, 'gap', info)
    plt.close()

    if debug: print('    saved gap plot:', time.time() - start_time)
    start_time = time.time()

    del time_gap_data

    return

def make_nodeID_plot(info, tree_data, tree_info):
    """
    For each instance create a plot, comparing the size of the gap at each node
    :param info: information about the instance
    :return:
    """
    plotname = 'nodeID'
    start_time = time.time()

    if tree_data is None or tree_data.empty or not 'primal_is_upper' in tree_info or (tree_data.primalbound == np.NaN).all():
        print('    no plotable vbc data found')
        return
    tree_data.drop(columns=['node_hex','var'],axis=1,inplace=True) ## dropping unnecessary columns, needed in a following line
    # check if the primal is a lower or upper bound
    if tree_info['primal_is_upper'] == 'Ambiguous':
        print('    the vbc file lists the primal bound as both upper and lower; will skip this plot')
        return
    elif tree_info['primal_is_upper']:
        tree_data['infeasible'] = (tree_data.dualbound > tree_data.primalbound).astype('bool')
    elif not tree_info['primal_is_upper']:
        tree_data['infeasible'] = (tree_data.dualbound < tree_data.primalbound).astype('bool')

    # calculate the gap
    if params['dualoptdiff'] and not tree_data.loc[is_fin(tree_data.primalbound), 'primalbound'].dropna().empty:
        if tree_info['primal_is_upper']:
            tree_data['gap'] = abs(tree_data.loc[is_fin(tree_data.primalbound), 'primalbound'].dropna().min() - tree_data.dualbound).astype('float')
        else:
            tree_data['gap'] = abs(tree_data.loc[is_fin(tree_data.primalbound), 'primalbound'].dropna().max() - tree_data.dualbound).astype('float')
    else:
        tree_data['gap'] = abs(tree_data.primalbound - tree_data.dualbound).astype('float')
    max_gap = tree_data.gap[-tree_data.infeasible & is_fin(tree_data.primalbound) & is_fin(tree_data.dualbound)].max()
    tree_data.loc[-tree_data.infeasible, 'gap'] = tree_data.gap[-tree_data.infeasible] / max_gap
    tree_data.loc[tree_data.infeasible, 'gap'] = 0.
    tree_data.gap = tree_data.gap.fillna(1.)
    tree_data['new_pb'] = (tree_data['primalbound'].sort_values(ascending = not tree_info['primal_is_upper']).diff() != 0).sort_index()
    try:
        mean_data = tree_data.sort_values('node_scip').rolling(5, center = True).mean()
    except TypeError:
        print("Could not calculate means. Skipping.")
        return
    if len(mean_data) >= 20:
        mean_data = mean_data.rolling(max(int(len(mean_data) * .08), 5)).mean()

    if debug: print('    extracted ' + plotname + ' data:', time.time() - start_time)
    start_time = time.time()

    # set x and y values
    x = [tree_data.node_scip[(-tree_data.infeasible) & (-tree_data.new_pb)].values,
        tree_data.node_scip[(-tree_data.infeasible) & (tree_data.new_pb)].values,
        tree_data.node_scip[(tree_data.infeasible) & (-tree_data.new_pb)].values,
        tree_data.node_scip[(tree_data.infeasible) & (tree_data.new_pb)].values]
    y = [(1 - tree_data.gap)[(-tree_data.infeasible) & (-tree_data.new_pb)].values,
        (1 - tree_data.gap)[(-tree_data.infeasible) & (tree_data.new_pb)].values,
        (1 - tree_data.gap)[(tree_data.infeasible) & (-tree_data.new_pb)].values,
        (1 - tree_data.gap)[(tree_data.infeasible) & (tree_data.new_pb)].values]
    x_mean = mean_data.index.values
    y_mean = (1 - mean_data.gap).values

    # create the figure
    fig = plt.gcf()
    ax = fig.gca()

    # format the plot
    ax.text(.0, 1.03, '\\textbf{\\underline{' + info['instance'].replace('_','\_') + '}}', ha = 'left', size = 'large', transform = ax.transAxes)
    description = 'Gap between incumbent solution and'
    if params['dualoptdiff']:
        description +=' \\textit{global (best)}'
    else:
        description +=' \\textit{local (current)}'
    description += ' dual bound of RMP \n vs node in the branch-and-bound tree'
    ax.text(.2, 1.0325, description, ha = 'left', va = 'center', size = 'medium', transform = ax.transAxes)
    ax.text(.75, 1.03, '\\textbf{Settings:} \\textit{' + info['settings'].replace('_','\_') + '}', ha = 'right', size = 'medium', transform = ax.transAxes)
    ax.text(1., 1.03, '\\textbf{SCIP Status:} \\textit{' + info['status'].replace('_',' ') + '}', ha = 'right', size = 'medium', transform = ax.transAxes)
    ax.set_xlabel('Node ID', size = 'large')
    ax.set_ylabel('Gap closed', size = 'large')
    ax.set_ylim([-0.04,1.04])
    ax.tick_params(axis='both', labelsize='large')
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer = True, min_n_ticks = 1))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base = .1))

    # plot the data
    ax.scatter(x[0], y[0], color = 'black', label="primal feasible, not improving the primal bound (but possibly the dual bound)")
    ax.scatter(x[1], y[1], color = 'blue', marker = 'd', label="primal feasible, improving (i.e. lowering) the \\textit{local} primal bound")
    ax.scatter(x[2], y[2], color = 'gray', marker = 'x', label="primal infeasible, not improving the primal bound")
    ax.scatter(x[3], y[3], color = 'green', marker = 'D', label="complete gap was closed (duality gap = 0\%)")
    ax.plot(x_mean, y_mean, 'k--', label='mean')
    ax.legend(loc='best')

    if debug: print('    plotted ' + plotname + ' data:', time.time() - start_time)
    start_time = time.time()

    # save the plot
    fig.set_size_inches(11.7,8.3)
    fig.tight_layout()
    fig.subplots_adjust(top = .925)
    save_plot(fig, plotname, info)
    plt.close()

    if debug: print('    saved ' + plotname + ' plot:', time.time() - start_time)
    start_time = time.time()

    return

def make_depth_plot(info, tree_data, tree_info):
    """
    For each instance create a plot, comparing the size of the gap at each node
    :param info: information about the instance
    :return:
    """
    start_time = time.time()

    if tree_data is None or tree_data.empty or not 'primal_is_upper' in tree_info or (tree_data.primalbound == np.NaN).all():
        print('    no plotable vbc data found')
        return

    # check if the primal is a lower or upper bound
    if tree_info['primal_is_upper'] == 'Ambiguous':
        print('    the vbc file lists the primal bound as both upper and lower; will skip this plot')
        return
    elif tree_info['primal_is_upper']:
        tree_data['infeasible'] = (tree_data.dualbound > tree_data.primalbound).astype('bool')
    elif not tree_info['primal_is_upper']:
        tree_data['infeasible'] = (tree_data.dualbound < tree_data.primalbound).astype('bool')

    # calculate the gap
    if params['dualoptdiff'] and not tree_data.loc[is_fin(tree_data.primalbound), 'primalbound'].dropna().empty:
        if tree_info['primal_is_upper']:
            tree_data['gap'] = abs(tree_data.loc[is_fin(tree_data.primalbound), 'primalbound'].dropna().min() - tree_data.dualbound).astype('float')
        else:
            tree_data['gap'] = abs(tree_data.loc[is_fin(tree_data.primalbound), 'primalbound'].dropna().max() - tree_data.dualbound).astype('float')
    else:
        tree_data['gap'] = abs(tree_data.primalbound - tree_data.dualbound).astype('float')
    max_gap = tree_data.gap[-tree_data.infeasible & is_fin(tree_data.primalbound) & is_fin(tree_data.dualbound)].max()
    tree_data.loc[-tree_data.infeasible, 'gap'] = tree_data.gap[-tree_data.infeasible] / max_gap
    tree_data.loc[tree_data.infeasible, 'gap'] = 0.
    tree_data.gap = tree_data.gap.fillna(1.)
    tree_data['new_pb'] = (tree_data['primalbound'].sort_values(ascending = not tree_info['primal_is_upper']).diff() != 0).sort_index()
    mean_data = tree_data.sort_values('depth').groupby('depth').mean()
    if len(mean_data) >= 20:
        mean_data = mean_data.rolling(max(int(len(mean_data) * .08), 5)).mean()

    if debug: print('    extracted depth data:', time.time() - start_time)
    start_time = time.time()

    # set x and y values
    x = [tree_data.depth[(-tree_data.infeasible) & (-tree_data.new_pb)].values,
        tree_data.depth[(-tree_data.infeasible) & (tree_data.new_pb)].values,
        tree_data.depth[(tree_data.infeasible) & (-tree_data.new_pb)].values,
        tree_data.depth[(tree_data.infeasible) & (tree_data.new_pb)].values]
    y = [(1 - tree_data.gap)[(-tree_data.infeasible) & (-tree_data.new_pb)].values,
        (1 - tree_data.gap)[(-tree_data.infeasible) & (tree_data.new_pb)].values,
        (1 - tree_data.gap)[(tree_data.infeasible) & (-tree_data.new_pb)].values,
        (1 - tree_data.gap)[(tree_data.infeasible) & (tree_data.new_pb)].values]
    x_mean = mean_data.index.values
    y_mean = (1 - mean_data.gap).values

    # create the figure
    fig = plt.gcf()
    ax = fig.gca()

    # format the plot
    ax.text(.0, 1.03, '\\textbf{\\underline{' + info['instance'].replace('_','\_') + '}}', ha = 'left', size = 'large', transform = ax.transAxes)
    description = 'Gap between incumbent solution and'
    if params['dualoptdiff']:
        description +=' \\textit{best}'
    else:
        description +=' \\textit{current}'
    description += ' local dual bound of RMP \n vs node in the branch-and-bound tree'
    ax.text(.25, 1.0325, description, ha = 'left', va = 'center', size = 'medium', transform = ax.transAxes)
    ax.text(.75, 1.03, '\\textbf{Settings:} \\textit{' + info['settings'].replace('_','\_') + '}', ha = 'right', size = 'medium', transform = ax.transAxes)
    ax.text(1., 1.03, '\\textbf{SCIP Status:} \\textit{' + info['status'].replace('_',' ') + '}', ha = 'right', size = 'medium', transform = ax.transAxes)
    ax.set_xlabel('Node Depth', size = 'large')
    ax.set_ylabel('Gap closed', size = 'large')
    ax.set_ylim([-0.04,1.04])
    ax.tick_params(axis='both', labelsize='large')
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer = True, min_n_ticks = 1))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base = .1))

    # plot the data
    ax.scatter(x[0], y[0], color = 'black', label="primal feasible, not improving the primal bound (but possibly the dual bound)")
    ax.scatter(x[1], y[1], color = 'blue', marker = 'd', label="primal feasible, improving (i.e. lowering) the \\textit{local} primal bound")
    ax.scatter(x[2], y[2], color = 'gray', marker = 'x', label="primal infeasible, not improving the primal bound")
    ax.scatter(x[3], y[3], color = 'green', marker = 'D', label="complete gap was closed (duality gap = 0\%)")
    ax.plot(x_mean, y_mean, 'k--', label='mean')
    ax.legend(loc='best')
    if debug: print('    plotted depth data:', time.time() - start_time)
    start_time = time.time()

    # save the plot
    fig.set_size_inches(11.7,8.3)
    fig.tight_layout()
    fig.subplots_adjust(top = .925)
    save_plot(fig, 'depth', info)
    plt.close()

    if debug: print('    saved depth plot:', time.time() - start_time)
    start_time = time.time()

    return

def make_vartimes_plot(settings, incumbent_times_tot, rootlpsol_times_tot):
    """
    Plot a distribution of relative creation times of master variables that are present in incumbent and root lp solutions
    """

    # set some colors
    rootlpsol_color = 'blue'
    incumbent_color = 'green'

    print('summarize variable creation times')

    start_time = time.time()

    # make the histogram
    fig = plt.gcf()
    ax = plt.gca()
    if not params['png']:
        lw = 0.01
    else:
        lw = 1.0

    ax.hist([rootlpsol_times_tot, incumbent_times_tot], 10, stacked=True, color = [rootlpsol_color, incumbent_color], edgecolor = 'white', alpha=0.75, label = ['Root LP Sol', 'Incumbent'])

    # set parameters
    if not params['png']:
        fig.set_size_inches(8*11.7,8*8.3)
        textsize = ax.get_window_extent().height * 0.013
    else:
        fig.set_size_inches(11.7,8.3)
        textsize = 12

    # formatting
#    ax.set_ylim([0.0, 1.0])
    ax.set_xlim([0.0, 1.0])
#    ax.get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True, nbins = 15))
    ax.tick_params(axis = 'both', length = textsize/2, width = textsize/40, labelsize = textsize*0.9, pad = 15)
    ax.set_xlabel('Relative solution time', size = 1.15*textsize)
    ax.set_ylabel('\# of variables', size = 1.15*textsize)
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

    # draw a legend
    handles = [lines.Line2D([0,0], [0,1], color = rootlpsol_color, linewidth = 4., label = 'Root LP Sol'), lines.Line2D([0,0], [0,1], color = incumbent_color, linewidth = 4., label = 'Incumbent')]
    plt.legend(handles = handles, bbox_to_anchor = (1.02, .915), loc = 2, fontsize = textsize)

    # save the figure
    plt.tight_layout()
    fig.subplots_adjust(top = 0.98, right = 0.85, left = 0.04)
    save_plot(fig, 'vartimes', {'instance': "all", 'settings': settings, 'status': "None"})
    plt.close()

    if debug: print('    save:', time.time() - start_time)
    if debug: print('    saved figure')

    return

def plots(data, info, incumbent_times, rootlpsol_times, root_bounds = None):
    """
    Master-function for plotting. Splits the data if necessary and calls all plotting functions (or a subset, according to the params)
    :param data: collected data as dataframe
    :param info: dictionary containing information about the data like the name of the instance, the settings & the scip_status
    :return:
    """
    print("    starting plotting...")

    # use tex to render the text output
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    # prepare gap data, if necessary
    if (not root_bounds is None) and (not params['no_gap'] or (not params['no_complete'] and params['gapincomplete'])):
        if debug: print("Collecting gap data")
        gap_data = collect_gap_data(data, root_bounds)
    else:
        gap_data = None


    # plot everything involving root bounds and the b&b-tree if possible and requested
    if not (root_bounds is None or params['no_gap']):
        make_gap_plot(gap_data, info)

    if not params['no_depth'] or not params['no_nodeID']:
        try:
            #print("Reading .vbc file...")
            tree_data, tree_info = vbc.read(params['vbcdir'] + '/' + info['instance'] + '.' + info['settings'] + '.vbc')
        except TypeError:
            print("    Warning: no vbc file found for instance {}. Skipping.".format(info['instance']))
            tree_data = None
            tree_info = None
            return
    if not params['no_depth']:
        make_depth_plot(info, tree_data, tree_info)
    if not params['no_nodeID']:
        make_nodeID_plot(info, tree_data, tree_info)

    # set parameters determining which part of the data shall be plotted
    if params['maxround'] <= 0:
        maxRnd = data.index.get_level_values('pricing_round').max()
    else:
        maxRnd = min(params['maxround'], data.index.get_level_values('pricing_round').max())
    if params['minround'] <= 1:
        minRnd = data.index.get_level_values('pricing_round').min()
    else:
        minRnd = max(params['minround'], data.index.get_level_values('pricing_round').min())

    if params['maxnode'] <= 0:
        maxNode = data.index.get_level_values('node').max()
    else:
        maxNode = min(params['maxnode'], data.index.get_level_values('node').max())
    if params['minnode'] <= 1:
        minNode = data.index.get_level_values('node').min()
    else:
        minNode = max(params['minnode'], data.index.get_level_values('node').min())
    if params['minnode'] > 1 or params['maxnode'] > 0:
        info['nodes_min'] = minNode
        info['nodes_max'] = maxNode

    if params['splitrounds'] <= 0:
        # do not split the plots, but still check if rounds or nodes were neglected
        if params['minround'] > 1 or params['maxround'] > 0:
            info['rounds_min'] = minRnd
            info['rounds_max'] = maxRnd
        data = data.query('@minNode <= node <= @maxNode & @minRnd <= pricing_round <= @maxRnd')
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
            gc.collect()
            # do not build the complete plot
            make_complete_plot(data, info, gap_data, incumbent_times, rootlpsol_times)

    else:
        # split the plots by rounds
        fromRnd = minRnd - 1
        for i in range(1,(maxRnd-minRnd)+1):
            if i % params['splitrounds'] != 0 and i != (maxRnd-minRnd):
                continue
            toRnd = i+minRnd
            data = data.query('@fromRnd < pricing_round <= @toRnd & @minNode <= node <= @maxNode')
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
                gc.collect()
                # do not build the complete plot
                make_complete_plot(data, info, gap_data, incumbent_times, rootlpsol_times)

            fromRnd = toRnd

    return

def collect_vartimes(data, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot):
    """
    Collects variable creation times of incumbent and root lp solutions over all instances
    """
    total_time = data.time.sum()

    incumbent_times_tot += [min((x / total_time if total_time != 0 else np.inf), 1.0) for x in incumbent_times]
    rootlpsol_times_tot += [min((x / total_time if total_time != 0 else np.inf), 1.0) for x in rootlpsol_times]

def load_data(files):
    """
    Plots data, that was parsed and collected earlier from the generate_files() method and saved in pickle-files
    :param files: the pickle files, from which the dataframes are to load (one file per instance)
    :return:
    """

    # global statistics over all files
    settings_global = 'default'
    incumbent_times_tot = list()
    rootlpsol_times_tot = list()

    for file in files:
        # check if the file exists
        if not os.path.exists(file):
            print('there is no file', file)
            continue
        filename, ext = os.path.splitext(os.path.basename(file))
        # extension has to be pkl
        if not (ext == '.pkl'):
            print(file, 'is not a pickle file')
            continue
        # check if the user wants to skip the instance
        if params['instances'] != '' and not any([(string in filename) for string in params['instances']]):
            print('skipping', filename)
            continue

        start_time = time.time()

        # restore the objects from the pickle file
        pkl_file = open(file, 'rb')
        objects = pickle.load(pkl_file)
        df = objects['pricing_data']
        info = objects['info']
        incumbent_times = objects['incumbent_times']
        rootlpsol_times = objects['rootlpsol_times']
        if 'root_bounds' in objects and not objects['root_bounds'].empty:
            root_bounds = objects['root_bounds']
        else:
            root_bounds = None
        pkl_file.close()

        print('entering', info['instance'])
        if debug: print('    loading data:', time.time() - start_time)
        if df.empty:
            print('    no data found')
        else:
            settings_global = info['settings']
            collect_vartimes(df, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot)
            start_time = time.time()
            # call the plotting master method
            plots(df, info, incumbent_times, rootlpsol_times, root_bounds)
            print('    total plotting:', time.time() - start_time)
        print('    leaving', info['instance'])

    if not params['no_vartimes']:
        make_vartimes_plot(settings_global, incumbent_times_tot, rootlpsol_times_tot)

    return

def collect_data(info, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot, root_bounds = None):
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
    # temporary workaround for damaged output at start of outfiles during cluster runs (8/2021)
    if get_info_from_filename(info["file"])["settings"].startswith("default_"):
        info["settings"] = get_info_from_filename(info["file"])["settings"]
    # build a dataframe from the index and value lists
    index = pd.MultiIndex.from_arrays([ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob],
                                      names=["node", "pricing_round", "stab_round", "round", "pricing_prob"])
    data = {'time': val_time, 'nVars': val_nVars, 'farkas': val_farkas}
    df = pd.DataFrame(data = data, index = index)

    # treat an empty root bounds table as no root bounds
    if not root_bounds is None and root_bounds.empty:
        root_bounds = None

    # save or plot the data
    if params['savepickle']:
        start_time = time.time()
        output = open(params['outdir'] + '/' + info['instance'] + '.' + info['settings'] + '.pricing' + '.pkl', 'wb')
        if root_bounds is None:
            pickle.dump({'pricing_data': df, 'info': info, 'incumbent_times': incumbent_times, 'rootlpsol_times': rootlpsol_times}, output, -1)
        else:
            pickle.dump({'pricing_data': df, 'info': info, 'incumbent_times': incumbent_times, 'rootlpsol_times': rootlpsol_times, 'root_bounds': root_bounds}, output, -1)
        output.close()
        print('    total saving:', time.time() - start_time)
    else:
        collect_vartimes(df, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot)
        start_time = time.time()
        plots(df, info, incumbent_times, rootlpsol_times, root_bounds)
        print('    total plotting:', time.time() - start_time)
    return

def parse_files(files):
    """
    Parse the (out-)files and structure the pricing-data in a dataframe
    :param files: List of files to be parsed
    :return:
    """

    # global statistics over all files
    settings_global = 'default'
    incumbent_times_tot = list()
    rootlpsol_times_tot = list()

    for file in files:
        with open(file) as _file:
            first_line_of_file = True
            done = True
            line_count_for_instance = 0

            for line in _file:
                line_count_for_instance += 1
                if line.find("@0")!=-1 or first_line_of_file:
                    # if the file is an out-file, generated by the check-script, reset the variables whenever a new instance starts
                    if line.startswith("@01") or first_line_of_file:
                        # print message, if the previous problem is not done yet
                        if not done and problemFileName:
                            if not ind_node or not ind_pricing_round or not ind_stab_round or not ind_pricing_prob or not val_time or not val_nVars:
                                print('    no pricing data found')
                                done = True
                                continue
                            print('    ended abruptly')
                            collect_data({'file': file, 'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot)
                            print('    leaving', problemFileName)
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
                        line_count_for_instance = 0

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
                        incumbent_times = list()
                        rootlpsol_times = list()
                    else:
                        # ignore lines, where the output ends abrubtly (e.g. when the hard limit of the check-script is reached)
                        continue

                # if the problem is already processed, continue
                elif done:
                    continue

                # if the line counter exceeds the limit (i.e. the instance is to big)
                elif line_count_for_instance > linelimit:
                        print("Warning: Line count for instance " + str(problemFileName) + " was exceeded. Stopping to collect data for it.\nWarning: Results of this plot may have limited expressive power.")
                        done = True
                        continue

                elif line.startswith("loaded parameter file"):
                    # store current settings
                    settings = line.split()[-1]
                    settings = settings.split("/")[-1]
                    settings = os.path.splitext(settings)[0]
                    settings_global = settings

                elif not problemFileName and line.startswith("read problem "):
                    # get the problem name from the file name as in "check.awk"
                    tmparray = line.split("<")[-1].replace(">","").replace("\n","").split("/")[-1].split(".")
                    problemFileName = tmparray[0]
                    if tmparray[-1] == "gz" or tmparray[-1] == "z" or tmparray[-1] == "GZ" or tmparray[-1] == "Z":
                        tmparray.pop()
                    for i in range(1,len(tmparray)-1):
                        problemFileName += "." + tmparray[i]
                    if params['instances'] != '' and not any([(string in problemFileName) for string in params['instances']]):
                        print('skipping', problemFileName)
                        done = True
                        continue
                    print('entering', problemFileName)

                # end of initial farkas pricing
                elif line.startswith("Starting reduced cost pricing..."):
                    farkasDone = True

                # read the SCIP status; end of pricing statistics
                elif line.startswith("SCIP Status        :"):
                    # continue if no data is found
                    if not ind_node or not ind_pricing_round or not ind_stab_round or not ind_pricing_prob or not val_time or not val_nVars:
                        print('    no pricing data found')
                        done = True
                        continue
                    scip_status = line.split(':')[1].split('[')[1].split(']')[0].replace(' ','_')
                    continue

                # read variable creation times of incumbent and root solution
                elif line.startswith("VAR:"):
                    tmparray = line.split()
                    time = tmparray[3]
                    if time == 'time':
                        continue
                    solval = float(tmparray[8])
                    rootlpsolval = float(tmparray[9])
                    if time == 'time':
                        continue
                    if solval > 0.0:
                        incumbent_times.append(float(time))
                    if rootlpsolval > 0.0:
                        rootlpsol_times.append(float(time))

                # read the root bounds table
                elif line.startswith("Root bounds"):
                    root_bounds = True

                elif root_bounds:
                    if line.startswith('Pricing Summary:'):
                        root_bounds = False
                    elif line.startswith('iter') or line.startswith(' iter'):
                        # create a dataframe to store the table
                        root_bounds_data = pd.DataFrame(columns = line.split())
                        root_bounds_ind = 0
                    else:
                        root_bounds_data.loc[root_bounds_ind] = [float(s) for s in line.split()]
                        root_bounds_ind += 1

                elif line.startswith("Pricing time in colpool"):
                    # nothing more to read for this instance
                    try: root_bounds_data.iter = root_bounds_data.iter.astype('int')
                    except AttributeError: print("Fatal: Could not read data for instance {}. Have you tested with STATISTICS=true?\nTerminating.".format(problemFileName))
                    collect_data({'file': file, 'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot, root_bounds = root_bounds_data)
                    done = True
                    print('    leaving', problemFileName)
                    continue

                # ignore all other lines, that do not contain pricer statistics messages
                elif not line.startswith("[pricer_gcg.cpp:"):
                    continue

                # extract the pricing-statistics message
                message = line.split("] statistic: ")[-1]

                if message.startswith("New pricing round at node") or message.startswith("New pr, node"):
                    try:
                        node = int(message.split()[-1])
                        pricing_round += 1
                        round_counter += 1
                        stab_round = 0
                        aggr_time = 0.0
                        aggr_nVars = 0
                        round_begin = True
                    except ValueError:
                        print('    ended abruptly')
                        collect_data({'file': file, 'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot)
                        print('    leaving', problemFileName)
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
                            if params['aggregate']:
                                # store all indices
                                ind_node.append(node)
                                ind_pricing_round.append(pricing_round)
                                ind_stab_round.append(stab_round)
                                ind_round.append(round_counter)
                                ind_pricing_prob.append(pricing_prob)
                                # store the data
                                val_time.append(aggr_time)
                                val_nVars.append(aggr_nVars)
                                val_farkas.append(not farkasDone)

                            lptime_end = float(message.split()[-1])
                            if lptime_end - lptime_begin > 0.005:
                                print('It seems, that the LP time is not constant during a pricing round. Delta t is', lptime_end - lptime_begin)
                    except ValueError:
                        print('    ended abruptly')
                        collect_data({'file': file, 'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot)
                        print('    leaving', problemFileName)
                        done = True
                        continue

                elif message.startswith("Stabilization round ") or message.startswith("Sr "):
                    try:
                        if params['aggregate']:
                            # store all indices
                            ind_node.append(node)
                            ind_pricing_round.append(pricing_round)
                            ind_stab_round.append(stab_round)
                            ind_round.append(round_counter)
                            ind_pricing_prob.append(pricing_prob)
                            # store the data
                            val_time.append(aggr_time)
                            val_nVars.append(aggr_nVars)
                            val_farkas.append(not farkasDone)
                            aggr_time = 0.0
                            aggr_nVars = 0

                        stab_round = int(message.split()[-1])
                        round_counter += 1
                    except ValueError:
                        print('    ended abruptly')
                        collect_data({'file': file, 'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot)
                        print('    leaving', problemFileName)
                        done = True
                        continue

                elif message.startswith("cp: ") or message.startswith("found "):
                    try:
                        if int(message.split()[1]) > 0:
                            # check if the column pool output should be included in the data
                            if node < params['minnode'] or (0 < params['maxnode'] < node) or pricing_round < params['minround'] or (0 < params['maxround'] < pricing_round):
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
                        print('    ended abruptly')
                        collect_data({'file': file, 'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot)
                        print('    leaving', problemFileName)
                        done = True
                        continue

                elif message.startswith("Pricing prob ") or message.startswith("P p "):
                    try:
                        if params['aggregate']:
                            pricing_prob = 0
                        else:
                            pricing_prob = int(message.split()[2])

                        # check if the pricing prob should be included in the data
                        if node < params['minnode'] or (0 < params['maxnode'] < node) or pricing_round < params['minround'] or (0 < params['maxround'] < pricing_round):
                            continue

                        if params['aggregate']:
                            # sum up the data over all pricing problems
                            aggr_time += float(message.split()[-1])
                            if message.startswith("P p "):
                                aggr_nVars += int(message.split()[-3])
                            else:
                                aggr_nVars += int(message.split()[5])
                        else:
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
                        print('    ended abruptly')
                        collect_data({'file': file, 'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot)
                        print('    leaving', problemFileName)
                        done = True
                        continue

            if not done:
                collect_data({'file': file, 'instance': problemFileName, 'settings': settings, 'status': scip_status}, ind_node, ind_pricing_round, ind_stab_round, ind_round, ind_pricing_prob, val_time, val_nVars, val_farkas, incumbent_times, rootlpsol_times, incumbent_times_tot, rootlpsol_times_tot)

    if not params['no_vartimes']:
        make_vartimes_plot(settings_global, incumbent_times_tot, rootlpsol_times_tot)

def main():
    """Entry point when calling this script"""
    args = sys.argv[1:]
    parsed_args = parse_arguments(args)
    set_params(parsed_args)
    if not os.path.exists(params['outdir']):
        os.makedirs(params['outdir'])
    if not os.path.exists(params['vbcdir']):
        os.makedirs(params['vbcdir'])
    if params['loadpickle'] or params['filenames'][0].endswith(".pkl"):
        load_data(parsed_args.filenames)
    elif params['vbconly']:
        # use tex to render the text output
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        for file in parsed_args.filenames:
            # read the tree data from a vbc file
            try:
                tree_data, tree_info = vbc.read(params['vbcdir'] + '/' + info['instance'] + '.' + info['settings'] + '.vbc')
            except TypeError:
                print("    Warning: no vbc file found for instance {}. Skipping.".format(info['instance']))
                tree_data = None
                tree_info = None
                continue
            try:
                make_nodeID_plot({'instance': os.path.splitext(os.path.basename(file))[0], 'settings': 'unknown', 'status': 'unknown'}, tree_data, tree_info)
                make_depth_plot({'instance': os.path.splitext(os.path.basename(file))[0], 'settings': 'unknown', 'status': 'unknown'}, tree_data, tree_info)
            except:
                print("    Warning: Unknown error. Skipping {}".format(file))
    else:
        parse_files(parsed_args.filenames)
        if params['savepickle']:
            logfile = open(params['outdir'] + '/origin.log','w')
            for f in parsed_args.filenames:
                logfile.write(f + '\n')
            logfile.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M'))
            logfile.close()

# Calling main script
if __name__ == '__main__':
#    try:
        main()
#    except KeyboardInterrupt:
#        print('KeyboardInterrupt.\nTerminating.')
#        try:
#            sys.exit(0)
#        except SystemExit:
#            os._exit(0)
