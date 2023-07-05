#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Zachary Sethna
"""

from __future__ import print_function, division
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection


from tcr_utils import ClonesAndCounts, ClonesAndPvals
from TcellRepertoire import TcellRepertoire
from CloneTrack import CloneTrack

def plot_detection_limit(ax, xmin = None, ymin = None, detection_limit_type = 'patch'):
    
    ax_xlims = ax.get_xlim()
    ax_ylims = ax.get_ylim()
    
    if detection_limit_type == 'patch':
        if xmin is not None and ymin is not None:    
            y_undetected_zone = mpatches.Rectangle(xy = [ax_xlims[0], ax_ylims[0]], width = ax_xlims[1] - ax_xlims[0], height = ymin - ax_ylims[0], edgecolor = 'None', facecolor = 'k', alpha = 0.2, zorder = -1)
            ax.add_patch(y_undetected_zone)
    
            x_undetected_zone = mpatches.Rectangle(xy = [ax_xlims[0], ymin], width = xmin - ax_xlims[0], height = ax_ylims[1] - ax_ylims[0], edgecolor = 'None', facecolor = 'k', alpha = 0.2, zorder = -1)
            ax.add_patch(x_undetected_zone)
    
        elif ymin is not None:
            y_undetected_zone = mpatches.Rectangle(xy = [ax_xlims[0], ax_ylims[0]], width = ax_xlims[1] - ax_xlims[0], height = ymin - ax_ylims[0], edgecolor = 'None', facecolor = 'k', alpha = 0.2, zorder = -1)
            ax.add_patch(y_undetected_zone)
        elif xmin is not None:
            x_undetected_zone = mpatches.Rectangle(xy = [ax_xlims[0], ax_ylims[0]], width = xmin - ax_xlims[0], height = ax_ylims[1] - ax_ylims[0], edgecolor = 'None', facecolor = 'k', alpha = 0.2, zorder = -1)
            ax.add_patch(x_undetected_zone)
            
    elif detection_limit_type == 'line':
    
        if xmin is not None and ymin is not None:
            
            ax.plot([xmin, xmin, ax_xlims[1]], [ax_ylims[1], ymin, ymin], 'k--', lw =1, zorder = -1)
            
        elif ymin is not None:
            ax.plot(ax_xlims, [ymin, ymin], 'k--', lw =1, zorder = -1)
        elif xmin is not None:
            ax.plot([xmin, xmin], ax_ylims, 'k--', lw =1, zorder = -1)

#%%
def plot_scatter(x_rep, y_rep, c_dict = None, kwargs_for_scatter = {}, **kwargs):
    """Plots clone trajectories

    Parameters
    ----------
    x_rep : TcellRepertoire or ClonesAndCounts
    
    kwargs_for_plots : dict
        kwargs for the plotting of clone trajectories. Keywords will be added to
        or replace default values. If a c_dict is included the 'c' field will
        be overwritten to allow for coloring by the c_dict.
        Default is: dict(c = 'k', alpha = 0.5, lw = 2)
    kwargs_for_av_line : dict
        kwargs for the plotting of clone (geometric) mean line. Keywords will be
        added to or replace default values. Default is:
            dict(c = 'k',
                 alpha = 1,
                 lw = 2,
                 capsize = 2,
                 elinewidth = 2,
                 label = 'Geometric mean over n = %s clones'%(len(clones)))

    Optional Parameters (**kwargs)
    ------------------------------
    sample_order : iterable
        Order of samples for clone trajectories.



    Returns
    -------
    None : NoneType
        If keyword return_axes == False (Default), only None is returned
    ax, cbar, fig : tuple
        If keyword return_axes == True, the axes, colorbar, and figure are
        returned. If no c_dict is provided the cbar returned will be None,
        if ax is provided as a kwargs the fig returned will be None.

    """
    
    shared_clones_only = kwargs.get('shared_clones_only', False)

    fontsize = kwargs.get('fontsize', 16)
    tick_fontsize = kwargs.get('tick_fontsize', 14)

    fontsize_dict = dict(xlabel_fontsize = kwargs.get('label_fontsize', fontsize),
                         ylabel_fontsize = kwargs.get('label_fontsize', fontsize),
                         title_fontsize = kwargs.get('title_fontsize', fontsize),
                         cbar_label_fontsize = kwargs.get('label_fontsize', fontsize),
                         xtick_fontsize = kwargs.get('xtick_fontsize', tick_fontsize),
                         ytick_fontsize = kwargs.get('ytick_fontsize', tick_fontsize),
                     )
    for kw, val in kwargs.items():
        if kw in fontsize_dict: fontsize_dict[kw] = val

    d_kwargs_for_scatter = dict(
                                c = 'k',
                                s = 7,
                                alpha = 0.2,
                                )
    if not kwargs.get('use_edgecolor', False):
        d_kwargs_for_scatter['edgecolor'] = 'None'

    for kw, val in kwargs_for_scatter.items():
        d_kwargs_for_scatter[kw] = val
#    if 'c' not in kwargs_for_scatter and 'color' not in kwargs_for_scatter: kwargs_for_scatter['c'] = 'k'
#    if 's' not in kwargs_for_scatter and 'size' not in kwargs_for_scatter: kwargs_for_scatter['s'] = 7
#    if 'alpha' not in kwargs_for_scatter: kwargs_for_scatter['alpha'] = 0.2
#    if 'edgecolor' not in kwargs_for_scatter and not kwargs.get('use_edgecolor', False): kwargs_for_scatter['edgecolor'] = 'None'


    if type(x_rep) == type(TcellRepertoire()):
        x_clones_and_counts = x_rep.get_clones_and_counts(**kwargs)
    else:
        x_clones_and_counts = x_rep
    
    if type(y_rep) == type(TcellRepertoire()):
        y_clones_and_counts = y_rep.get_clones_and_counts(**kwargs)
    else:
        y_clones_and_counts = y_rep
    x_norm = x_clones_and_counts.norm
    y_norm = y_clones_and_counts.norm
    if shared_clones_only:
        clones = x_clones_and_counts.clone_intersection(y_clones_and_counts)
        x_to_plot = np.array([x_clones_and_counts[c] for c in clones])
        y_to_plot = np.array([y_clones_and_counts[c] for c in clones])
    else:
        clones = x_clones_and_counts.clone_union(y_clones_and_counts)
        x_to_plot = np.array([x_clones_and_counts[c] for c in clones])
        y_to_plot = np.array([y_clones_and_counts[c] for c in clones])
        x_to_plot[x_to_plot == 0] = 1/3.
        y_to_plot[y_to_plot == 0] = 1/3.


    if kwargs.get('normalize', False):
        x_to_plot /= x_norm
        y_to_plot /= y_norm


    if 'ax' in kwargs:
        ax = kwargs['ax']
    elif 'figsize' in kwargs:
        fig, ax = plt.subplots(figsize = kwargs['figsize'])
    elif c_dict is not None:
        fig, ax = plt.subplots(figsize = (7, 5))
    else:
        fig, ax = plt.subplots(figsize = (5.45, 5))


    if c_dict is not None:
        colors_to_plot = np.array([c_dict[c] for c in clones])
        missing_clones = [i for i, c in enumerate(clones) if c not in c_dict]
        kwargs_missing = {kw: val for kw, val in d_kwargs_for_scatter.items() if kw not in ['cmap', 'norm']}
        kwargs_color = {kw: val for kw, val in d_kwargs_for_scatter.items() if kw not in ['c', 'alpha']}
        if type(c_dict) == ClonesAndPvals:
            if 'cmap' not in kwargs_color: kwargs_color['cmap'] = 'cool_r'
            c_order = sorted([i for i, c in enumerate(clones) if c in c_dict], key = c_dict.get, reverse = True)
            if 'norm' not in kwargs_color:
                vmin = kwargs.get('vmin', 1e-6)
                vmax = kwargs.get('vmax', 1)
                colors_to_plot = np.clip(colors_to_plot, vmin, 1)
                kwargs_color['norm'] = mpl.colors.LogNorm(vmin = vmin, vmax = vmax)
        else:
            if 'cmap' not in kwargs_color: kwargs_color['cmap'] = 'plasma'
            c_order = sorted([i for i, c in enumerate(clones) if c in c_dict], key = c_dict.get)
            if 'norm' not in kwargs_color:
                if 'vmin' in kwargs and 'vmax' in kwargs:
                    vmin = kwargs['vmin']
                    vmax = kwargs['vmax']
                    colors_to_plot = np.clip(colors_to_plot, vmin, vmax)
                    if kwargs.get('colorscale_normalize', 'log') == 'log':
                        kwargs_color['norm'] = mpl.colors.LogNorm(vmin = vmin, vmax = vmax)
                    else:
                        kwargs_color['norm'] = mpl.colors.Normalize(vmin = vmin, vmax = vmax)
                else:
                    if kwargs.get('colorscale_normalize', 'log') == 'log':
                        kwargs_color['norm'] = mpl.colors.LogNorm()
                    else:
                        kwargs_color['norm'] = mpl.colors.Normalize()


        ax.scatter(x_to_plot[missing_clones], y_to_plot[missing_clones], **kwargs_missing)
        #ax.scatter(x_to_plot[missing_clones], y_to_plot[missing_clones], c = np.ones(len(missing_clones)), **kwargs_color)
        
        tcr_sc = ax.scatter(x_to_plot[c_order], y_to_plot[c_order], c = colors_to_plot[c_order], **kwargs_color)
        cbar = plt.colorbar(tcr_sc)
        if 'cbar_label' in kwargs:
            cbar.set_label(kwargs['cbar_label'], fontsize = fontsize_dict['cbar_label_fontsize'])

    else:
        ax.scatter(x_to_plot, y_to_plot, **d_kwargs_for_scatter)

    if (kwargs.get('split_clones', False) and kwargs.get('plot_clone_linkage', True)) == True: #No idea why the True == True test needed here... but buggy otherwise
        print('got here!')
        r_kwargs = kwargs
        r_kwargs['split_clones'] = False
        r_x_clones_and_counts = x_rep.get_clones_and_counts(**r_kwargs)
        r_y_clones_and_counts = y_rep.get_clones_and_counts(**r_kwargs)
        if len(r_x_clones_and_counts) < len(x_clones_and_counts):
            c_clone_def = x_rep.get_clone_def()
            for kw, kw_val in kwargs.items():
                if kw in c_clone_def: c_clone_def[kw] = kw_val
            for c_clone in x_rep.full_clone_list:
                c_split_clone_inds = [c.index(c) for c in c_clone.split_clone_rep(**c_clone_def) if c in clones]
                if len(c_split_clone_inds) > 1:
                    ax.plot([x_to_plot[c_ind] for c_ind in c_split_clone_inds], [y_to_plot[c_ind] for c_ind in c_split_clone_inds], 'k', lw = 1, zorder = 0)

        if len(r_y_clones_and_counts) < len(y_clones_and_counts):
            c_clone_def = y_rep.get_clone_def()
            for kw, kw_val in kwargs.items():
                if kw in c_clone_def: c_clone_def[kw] = kw_val
            for c_clone in y_rep.full_clone_list:
                c_split_clone_inds = [clones.index(c) for c in c_clone.split_clone_rep(**c_clone_def) if c in clones]
                if len(c_split_clone_inds) > 1:
                    ax.plot([x_to_plot[c_ind] for c_ind in c_split_clone_inds], [y_to_plot[c_ind] for c_ind in c_split_clone_inds], 'k', lw = 1, zorder = 0)

    dmin_val = min(np.min(x_to_plot), np.min(y_to_plot))
    dmax_val = max(np.max(x_to_plot), np.max(y_to_plot))
    if 'xlim' in kwargs:
        ax.set_xlim(kwargs['xlim'])
    else:
        ax.set_xlim([dmin_val/1.5, dmax_val*2])
    if 'ylim' in kwargs:
        ax.set_ylim(kwargs['ylim'])
    else:
        ax.set_ylim([dmin_val/1.5, dmax_val*2])


    if kwargs.get('plot_diagonal', False):
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        min_val = min(xlim[0], ylim[0])
        if kwargs.get('normalize', False):
            ax.plot([min_val, 1], [min_val, 1], 'k--', lw = 1)
        else:
            ax.plot([min_val*x_norm/(x_norm + y_norm), x_norm], [min_val*y_norm/(x_norm + y_norm), y_norm], 'k--', lw =1, zorder = 0)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.tick_params(axis='x', labelsize= fontsize_dict['xtick_fontsize'])
    ax.tick_params(axis='y', labelsize= fontsize_dict['ytick_fontsize'])
    if not shared_clones_only and kwargs.get('plot_detection_limit', True):
        # x_max = ax.get_xlim()[1]
        # y_max = ax.get_ylim()[1]
        if kwargs.get('normalize', False):
            x_min = 0.6/x_norm
            y_min = 0.6/y_norm
        else:
            x_min = 0.6
            y_min = 0.6
        #ax.plot([x_min, x_min,x_max], [y_max, y_min, y_min], 'k--', lw = 1)
        plot_detection_limit(ax = ax, xmin = x_min, ymin=y_min, detection_limit_type = kwargs.get('detection_limit_type', 'patch'))

    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'], fontsize = fontsize_dict['xlabel_fontsize'])
    elif len(x_rep.name) > 0:
        if kwargs.get('normalize', False):
            ax.set_xlabel('Frequency in %s'%(x_rep.name), fontsize = fontsize_dict['xlabel_fontsize'])
        else:
            ax.set_xlabel('Counts in %s'%(x_rep.name), fontsize = fontsize_dict['xlabel_fontsize'])

    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'], fontsize = fontsize_dict['ylabel_fontsize'])
    elif len(y_rep.name) > 0:
        if kwargs.get('normalize', False):
            ax.set_ylabel('Frequency in %s'%(y_rep.name), fontsize = fontsize_dict['ylabel_fontsize'])
        else:
            ax.set_ylabel('Counts in %s'%(y_rep.name), fontsize = fontsize_dict['ylabel_fontsize'])

    if 'title' in kwargs:
        ax.set_title(kwargs['title'], fontsize = fontsize_dict['title_fontsize'])

    if 'ax' not in kwargs:
        fig.tight_layout()

    if 'savefig_filename' in kwargs and 'ax' not in kwargs:
        fig.savefig(kwargs['savefig_filename'])

    if kwargs.get('return_axes', False):
        if c_dict is not None:
            return ax, cbar
        else:
            return ax
    else:
        return None

#%%
def plot_clone_trajectories(rep_time_series, clones = None, c_dict = None, kwargs_for_plots = {}, kwargs_for_av_line = {}, kwargs_for_tot_line = {}, **kwargs):
    """Plots clone trajectories

    Parameters
    ----------
    rep_time_series : RepertoireTimeSeries
        Repertoire time series object that clones trajectories will be pulled
        from
    clones : iterable or None (default)
        List of clones of trajectories to plot. If None the keys of c_dict will
        be used. If neither clones or c_dict are provided, all of the clones
        from rep_time_series will be used as default.
    
    kwargs_for_plots : dict
        kwargs for the plotting of clone trajectories. Keywords will be added to
        or replace default values. If a c_dict is included the 'c' field will
        be overwritten to allow for coloring by the c_dict.
        Default is: dict(c = 'k', alpha = 0.5, lw = 2)
    kwargs_for_av_line : dict
        kwargs for the plotting of clone (geometric) mean line. Keywords will be
        added to or replace default values. Default is:
            dict(c = 'k',
                 alpha = 1,
                 lw = 2,
                 capsize = 2,
                 elinewidth = 2,
                 label = 'Geometric mean over n = %s clones'%(len(clones)))

    Optional Parameters (**kwargs)
    ------------------------------
    sample_order : iterable
        Order of samples for clone trajectories.
    clone_ordering : 'ascending', 'descending'
        Sorts the clones according to the c_dict values. For pval colorings 
        descending is recommended (most significant clones will be plotted on 
        top), for frequency colorings ascending is recommended (largest on top).


    Returns
    -------
    None : NoneType
        If keyword return_axes == False (Default), only None is returned
    fig, ax, cbar : tuple
        If keyword return_axes == True, the axes, colorbar, and figure are
        returned. If no c_dict is provided the cbar returned will be None,
        if ax is provided as a kwargs the fig returned will be None.

    """


    cbar = None
    fig = None
    d_kwargs_for_plots = dict(c = 'k',
                             alpha = 0.5,
                             lw = 2
                            )

    for kw, val in kwargs_for_plots.items():
        d_kwargs_for_plots[kw] = val

    #set fontsizes
    fontsize = kwargs.get('fontsize', 12)
    tick_fontsize = kwargs.get('tick_fontsize', 10)
    fontsize_dict = dict(xlabel_fontsize = kwargs.get('label_fontsize', fontsize),
                         ylabel_fontsize = kwargs.get('label_fontsize', fontsize),
                         title_fontsize = fontsize,
                         cbar_label_fontsize = kwargs.get('label_fontsize', fontsize),
                         xtick_fontsize = kwargs.get('xtick_fontsize', tick_fontsize),
                         legend_fontsize = kwargs.get('legend_fontsize', fontsize)
                         )
    for kw, val in kwargs.items():
        if kw in fontsize_dict: fontsize_dict[kw] = val

    show_av = kwargs.get('plot_average', True)
    show_tot = kwargs.get('plot_total', True)
    sample_order = kwargs.get('sample_order', rep_time_series.sample_order)

    if kwargs.get('real_time_ax', True):
        time_ax = np.array([rep_time_series.sample_timepoints[sample]*kwargs.get('time_ax_rescale', 1) for sample in sample_order])
    else:
        time_ax = np.array(range(len(sample_order)))

    #set clone for plotting
    if clones is not None:
        pass
    elif c_dict is not None:
        clones = list(c_dict.keys())
    else:
        clones = rep_time_series.clones
        
#    if len(clones) < 4: show_av = False

    if 'ax' in kwargs:
        ax = kwargs['ax']
    elif 'figsize' in kwargs:
        fig, ax = plt.subplots(figsize = kwargs['figsize'])
    elif c_dict is not None and kwargs.get('plot_cbar', True):
        fig, ax = plt.subplots(figsize = (7, 5))
    else:
        fig, ax = plt.subplots(figsize = (5.45, 5))

    if c_dict is not None:
        if kwargs.get('clone_ordering', '') == 'descending':
            clones = sorted(clones, key = c_dict.get, reverse = True)
        elif kwargs.get('clone_ordering', '') =='ascending':
            clones = sorted(clones, key = c_dict.get)
        line_color_vals = np.array([c_dict[c] for c in clones])
        
        #set color info
        d_kwargs_for_plots.__delitem__('c')
        #d_kwargs_for_plots['c']  = line_color_vals
        if 'cmap' not in d_kwargs_for_plots: d_kwargs_for_plots['cmap'] = 'cool_r'
        if 'norm' not in d_kwargs_for_plots: d_kwargs_for_plots['norm'] = mpl.colors.LogNorm(vmin = kwargs.get('vmin', 1e-6), vmax = kwargs.get('vmax', 1))
        
        if kwargs.get('plot_cbar', True):
            c_sc_for_cbar = ax.scatter([], [], c = [], cmap = d_kwargs_for_plots['cmap'], norm = d_kwargs_for_plots['norm'])
            cbar = plt.colorbar(c_sc_for_cbar)
            cbar.set_label(kwargs.get('cbar_label', ''), fontsize = 14)
    
    if len(clones) > 0:
        c_freqs_to_plot = rep_time_series.get_timecourse(clones, trajectory_type = 'frequency_for_plot', sample_order = sample_order)

        if c_dict is None:
            for c_freqs in c_freqs_to_plot:
                ax.plot(time_ax, c_freqs, **d_kwargs_for_plots)
        else:
            lc = LineCollection([np.column_stack([x, y]) for x, y in zip([time_ax for _ in range(len(clones))], c_freqs_to_plot)],
                                 **d_kwargs_for_plots)
            lc.set_array(line_color_vals)
            ax.add_collection(lc)
        #cbar = plt.colorbar(lc)
    
    if len(clones) > 1:
        if show_av:
            c_mean_freq_to_plot = np.exp(np.mean(np.log(c_freqs_to_plot), axis = 0))
            c_std_freq_to_plot = np.exp(np.std(np.log(c_freqs_to_plot), axis = 0))
            c_error_arr = [c_mean_freq_to_plot - c_mean_freq_to_plot/c_std_freq_to_plot, c_mean_freq_to_plot*c_std_freq_to_plot -c_mean_freq_to_plot]

            d_kwargs_for_av_line = dict(c = 'k',
                               alpha = 1,
                               lw = 2,
                               capsize = 2,
                               elinewidth = 2,
                               label = 'Geometric mean (n = %s)'%(len(c_freqs_to_plot))
                               )
            if c_dict is None: d_kwargs_for_av_line['c'] = 'r'
            for kw, val in kwargs_for_av_line.items():
                d_kwargs_for_av_line[kw] = val
            ax.errorbar(time_ax, c_mean_freq_to_plot, yerr = c_error_arr,  **d_kwargs_for_av_line)
            
        if show_tot:
            
            tot_freq_timecourse = rep_time_series.get_timecourse(clones, trajectory_type = 'frequency_for_plot', sample_order = sample_order, aggregate = True).squeeze()
            
            d_kwargs_for_tot_line = dict(c = 'r',
                               alpha = 1,
                               lw = 2,
                               label = 'Total frequency (n = %s)'%(len(c_freqs_to_plot))
                               )
            for kw, val in kwargs_for_tot_line.items():
                d_kwargs_for_av_line[kw] = val
            ax.plot(time_ax, tot_freq_timecourse, **d_kwargs_for_tot_line)
        
        if kwargs.get('show_legend', True):
            ax.legend(fontsize = fontsize_dict['legend_fontsize'], frameon = False)

    max_sample_norm = max([rep_time_series[sample].norm for sample in sample_order])
    #ax.plot(range(c_freqs_to_plot.shape[1]), np.ones(c_freqs_to_plot.shape[1])* 1/(2*max_sample_norm), 'k--', lw = 1)
        
    # if len(clones) > 0:
    #     if show_tot and len(clones) > 1:
    #         ax.set_ylim([1/(5*max_sample_norm), max(1e-1, 1.2*np.max(tot_freq_timecourse))])
    #     else:
    #         ax.set_ylim([1/(5*max_sample_norm), max(1e-1, 1.2*np.max(c_freqs_to_plot))])
    # else:
    #     ax.set_ylim([1/(5*max_sample_norm), 1e-1])
    
    if 'xlim' in kwargs:
        ax.set_xlim(kwargs['xlim'])
    if 'ylim' in kwargs:
        ax.set_ylim(kwargs['ylim'])
    else:
        ax.set_ylim([1/(5*max_sample_norm), 1])
    
    if kwargs.get('plot_detection_limit', True):
        #ax.plot(time_ax, np.ones(len(sample_order))* 1/(2*max_sample_norm), 'k--', lw = 1)
        plot_detection_limit(ax = ax, ymin = 1/(2*max_sample_norm), detection_limit_type = kwargs.get('detection_limit_type', 'patch'))
        
    
    ax.set_yscale('log')

    if kwargs.get('use_rep_names_as_ticks', False):
        c_xtick_labels = [rep_time_series.sample_metadata[sample]['full_name'] for sample in sample_order]
        plt.xticks(time_ax, c_xtick_labels, rotation = -30, fontsize = fontsize_dict['xtick_fontsize'])
        plt.setp(ax.xaxis.get_majorticklabels(), ha="left", rotation_mode="anchor")

    #plt.xlabel('Time point', fontsize = 14)
    ax.set_ylabel('Clone frequency', fontsize = 14)

    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'], fontsize = fontsize_dict['xlabel_fontsize'])

    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'], fontsize = fontsize_dict['ylabel_fontsize'])
    else:
        ax.set_ylabel('Clone frequency', fontsize = fontsize_dict['ylabel_fontsize'])

    if 'title' in kwargs:
        ax.set_title(kwargs['title'], fontsize = fontsize_dict['title_fontsize'])

    if 'ax' not in kwargs:
        fig.tight_layout()
        if 'savefig_filename' in kwargs:
            fig.savefig(kwargs['savefig_filename'])

    if kwargs.get('return_axes', False):
        return  fig, ax, cbar

    else:
        return None

#%%
def plot_rank_frequency(clones_and_counts, kwargs_for_plots = {}, **kwargs):
    """Plots rank frequency

    Parameters
    ----------
    clones_and_counts : ClonesAndCounts or list of ClonesAndCounts
        Data to determine the ranks and frequencies
        
    kwargs_for_plots : dict or list of dicts
        kwargs for the plotting of individual curves. Keywords will be added to
        or replace default values. Colors will be determined by cycle by 
        default (cycle can be determined by color_cycle keyword), but this can 
        be overrided by including 'c' keyword.
        Default is: dict(alpha = 1, lw = 2)

    Optional Parameters (**kwargs)
    ------------------------------
    
    color_cycle : list
        Used for color cycle of plots

    Returns
    -------
    None : NoneType
        If keyword return_axes == False (Default), only None is returned
    fig, ax : tuple
        If keyword return_axes == True, the axes and figure are
        returned. If ax is provided as a kwargs the fig returned will be None.

    """

    fig = None

    if type(clones_and_counts) is not list:
        clones_and_counts = [clones_and_counts]
    
    d_kwargs_for_plots = [dict(alpha = 1,
                              lw = 2
                              ) for c_and_c in clones_and_counts]
    
    if type(kwargs_for_plots) is dict:
        kwargs_for_plots = [kwargs_for_plots for c_and_c in clones_and_counts]
        
    for i, c_kwargs_for_plots in enumerate(kwargs_for_plots):    
        for kw, val in c_kwargs_for_plots.items():
            d_kwargs_for_plots[i][kw] = val

    #set fontsizes
    fontsize = kwargs.get('fontsize', 12)
    tick_fontsize = kwargs.get('tick_fontsize', 10)
    fontsize_dict = dict(xlabel_fontsize = kwargs.get('label_fontsize', fontsize),
                         ylabel_fontsize = kwargs.get('label_fontsize', fontsize),
                         title_fontsize = fontsize,
                         cbar_label_fontsize = kwargs.get('label_fontsize', fontsize),
                         xtick_fontsize = kwargs.get('xtick_fontsize', tick_fontsize),
                         legend_fontsize = kwargs.get('legend_fontsize', fontsize)
                         )
    for kw, val in kwargs.items():
        if kw in fontsize_dict: fontsize_dict[kw] = val

    if 'ax' in kwargs:
        ax = kwargs['ax']
    elif 'figsize' in kwargs:
        fig, ax = plt.subplots(figsize = kwargs['figsize'])
    else:
        fig, ax = plt.subplots(figsize = (5.45, 5))

    if 'color_cyle' in kwargs:
        ax.set_prop_cycle(color=kwargs['color_cyle'])
    
    for c_and_c, c_d_kwargs_for_plots in zip(clones_and_counts, d_kwargs_for_plots):
        c_freqs = np.array(sorted(c_and_c.values(), reverse = True))/c_and_c.norm
        c_ranks = range(1, len(c_freqs) + 1)
        if kwargs.get('show_fit', True):
            slope_fit, intercept_fit, r_value_fit, p_value_fit, std_err_fit = linregress(np.array([np.log(c_ranks), np.log(c_freqs)]).T)
            #c_d_kwargs_for_plots['label'] = c_d_kwargs_for_plots.get('label', '') + ' (%.2e%sr^%.2e, %s = %.2e)'%(np.exp(intercept_fit), r'$\times$', slope_fit, r'$R^2$', r_value_fit)
            c_d_kwargs_for_plots['label'] = c_d_kwargs_for_plots.get('label', '') + ' (r^%.2f)'%(slope_fit)

            ax.plot(c_ranks, c_freqs, **c_d_kwargs_for_plots)
        else:
            ax.plot(c_ranks, c_freqs, **c_d_kwargs_for_plots)
        
            
        
    max_rank = max([len(c_and_c.clones()) for c_and_c in clones_and_counts])
    max_sample_norm = max([c_and_c.norm for c_and_c in clones_and_counts])
        
    if 'xlim' in kwargs:
        ax.set_xlim(kwargs['xlim'])
    else:
        ax.set_xlim([1, 10**(np.ceil(np.log10(max_rank)))])
    if 'ylim' in kwargs:
        ax.set_ylim(kwargs['ylim'])
    else:
        ax.set_ylim([10**(-np.ceil(np.log10(max_sample_norm))), 1])
        
    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_xlabel(kwargs.get('xlabel', 'Rank'), fontsize = fontsize_dict['xlabel_fontsize'])
    ax.set_ylabel(kwargs.get('ylabel', 'Frequency'), fontsize = fontsize_dict['ylabel_fontsize'])
    
    if 'title' in kwargs:
        ax.set_title(kwargs['title'], fontsize = fontsize_dict['title_fontsize'])

    if kwargs.get('show_legend', True):
        ax.legend(fontsize = fontsize_dict['legend_fontsize'], frameon = False)


    if 'ax' not in kwargs:
        fig.tight_layout()
        if 'savefig_filename' in kwargs:
            fig.savefig(kwargs['savefig_filename'])

    if kwargs.get('return_axes', False):
        return  fig, ax
    else:
        return None
    
#%%
