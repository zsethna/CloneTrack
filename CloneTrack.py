#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Zachary Sethna
"""

from __future__ import print_function, division
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

import json
import numpy as np
#from tcr_errors import InputError, ChainMismatch, UnknownGene, SequenceMismatch
from tcr_utils import ClonesAndCounts, ClonesAndPvals, TCRClonePvalue
#from tcr_seq_and_clone_definitions import CloneDefinition, TCRseq, TcellClone
from TcellRepertoire import TcellRepertoire

#%%
class CloneTrack(TCRClonePvalue):
    
    pseudo_count = 1/3
    def __init__(self, data_dir = None, **kwargs):
        self.samples = {}
        self.clones = []
        self.sample_metadata = {}
        self.sample_timepoints = {}
        self.sample_order = []
        self.timepoint_sample_dict = {}
        if data_dir is not None:
            with open(os.path.join(data_dir, 'metadata.json'), 'r') as metadata_infile:
                self.sample_metadata = json.load(metadata_infile)
            for sample, c_metadata in self.sample_metadata.items():
                #print(sample, c_metadata['filename'])
                self.load_adaptive_sample(sample, os.path.join(data_dir, c_metadata['filename']), timepoint = c_metadata.get('timepoint', None), **kwargs)
        self.pvalues = {}


    def __getitem__(self, sample):
        return self.samples[sample]
    
    def __repr__(self):
        
        
        pval_str = ', '.join(['%s: %s'%(pval_name, len(clones_and_pvals)) for pval_name, clones_and_pvals in self.pvalues.items()])
        
        ordered_sample_str = ', '.join(['%s: %s'%(sample_name, self.sample_timepoints[sample_name]) for sample_name in self.sample_order])
        unordered_sample_str = ', '.join([sample_name for sample_name in self.samples.keys() if sample_name not in self.sample_order])
        
        out_list = ['Total clones: %s'%(len(self.clones)),
                    'Significant clones: {%s}'%(pval_str),
                    'Timepoint samples: {%s}'%(ordered_sample_str),
                    'Additional samples: [%s]'%(unordered_sample_str)
                    ]
        return "RepertoireTimeSeries(%s)"%(',\n'.join(out_list))

    def load_adaptive_sample(self, name, infile_name, timepoint = None, **kwargs):

        self.samples[name] = TcellRepertoire(name = name, adaptive_filename = infile_name, **kwargs)
        self.clones = self[name].clone_union(self.clones)
        if timepoint is not None:
            self.timepoint_sample_dict[timepoint] = name
            self.sample_timepoints[name] = timepoint
            self.sample_order = sorted(self.sample_timepoints.keys(), key = self.sample_timepoints.get)

    def add_expansion_p_value(self, name, baseline_sample, last_sample = None, foldchange_thresh = 2, **kwargs):

        if last_sample is not None:
            reps_to_sweep = [self[sample].get_clones_and_counts(**kwargs) for sample in self.sample_order[self.sample_order.index(baseline_sample)+1:self.sample_order.index(last_sample)+1]]
        else:
            reps_to_sweep = [self[sample].get_clones_and_counts(**kwargs) for sample in self.sample_order[self.sample_order.index(baseline_sample)+1:]]

        baseline_rep = self[baseline_sample].get_clones_and_counts(**kwargs)

        clones_all = sum(reps_to_sweep, ClonesAndCounts()).clones()

        c_pvals_all = [self.compute_fisher_pvalues(c_rep, baseline_rep, foldchange_thresh = foldchange_thresh) for c_rep in reps_to_sweep]

        if kwargs.get('use_qc', False):
            clones_all = [c for c in clones_all if self.timecourse_qc(c, **kwargs)]

        if 'de_novo_count_thresh' in kwargs:
            #exclude clones with counts higher than de_novo_count_thresh in baseline or earlier timepoints
            prev_reps = [self[sample].get_clones_and_counts(**kwargs) for sample in self.sample_order[:self.sample_order.index(baseline_sample)]] + [baseline_rep]
            clones_all = [c for c in clones_all if all([c_rep[c] <= kwargs['de_novo_count_thresh'] for c_rep in prev_reps])]

        if 'de_novo_freq_thresh' in kwargs:
            #exclude clones with freqs higher than de_novo_freq_thresh in baseline or earlier timepoints
            prev_reps = [self[sample].get_clones_and_counts(**kwargs) for sample in self.sample_order[:self.sample_order.index(baseline_sample)]] + [baseline_rep]
            clones_all = [c for c in clones_all if all([c_rep.get_frequency(c) <= kwargs['de_novo_freq_thresh'] for c_rep in prev_reps])]

        exp_pvals = ClonesAndPvals()
        for clone in clones_all:
            exp_pvals[clone] = np.clip(min([c_rep_pval[clone] for c_rep_pval in c_pvals_all])*len(c_pvals_all), 0, 1)

        self.pvalues[name] = exp_pvals
            
    def get_timecourse(self, clones, trajectory_type = 'count', sample_order = None, aggregate = False):

        if sample_order is None: sample_order = self.sample_order

        if type(clones) == str: clones = [clones]

        if trajectory_type == 'count':
            if aggregate:
                return np.sum(np.array([[self[sample][clone] for sample in sample_order] for clone in clones]), axis = 0).reshape(1, -1)
            else:
                return np.array([[self[sample][clone] for sample in sample_order] for clone in clones])
        elif trajectory_type == 'pseudo_count':
            c_counts_array = self.get_timecourse(clones, trajectory_type = 'count', sample_order = sample_order, aggregate = aggregate)
            c_counts_array[c_counts_array == 0] = self.pseudo_count
            return c_counts_array
        else:
            timepoint_norms = np.array([self[sample].norm for sample in sample_order])
            if trajectory_type == 'frequency':
                return self.get_timecourse(clones, trajectory_type = 'count', sample_order = sample_order, aggregate = aggregate)/timepoint_norms
            elif trajectory_type == 'pseudo_frequency':
                return self.get_timecourse(clones, trajectory_type = 'pseudo_count', sample_order = sample_order, aggregate = aggregate)/timepoint_norms
            elif trajectory_type == 'frequency_for_plot':
                c_freq_array = self.get_timecourse(clones, trajectory_type = 'frequency', sample_order = sample_order, aggregate = aggregate)
                c_freq_array[c_freq_array == 0] = self.pseudo_count/max(timepoint_norms)
                return c_freq_array

    def timecourse_qc(self, clones, cell_count_thresh = 3, timepoint_thresh = 2, sample_order = None, **kwargs):

        return np.sum(self.get_timecourse(clones, sample_order= sample_order) >= cell_count_thresh, axis = 1) >= timepoint_thresh