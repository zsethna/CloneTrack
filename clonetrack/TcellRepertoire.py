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

from tcr_errors import InputError, ChainMismatch, UnknownGene, SequenceMismatch
from tcr_utils import ClonesAndCounts, TCRClonePvalue
from tcr_seq_and_clone_definitions import CloneDefinition, TCRseq, TcellClone


class TcellRepertoire(CloneDefinition, TCRClonePvalue, ClonesAndCounts):
    """Class for tracking the clones of a T cell repertoire sample.

    Inherits from ClonesAndCounts and thus will behave similarly to a defaultdict
    with default_value of 0. However, it will remove all clones with counts<=0
    and has support added for the magic method __add__ which can add two
    ClonesAndCount classes.

    Class Attributes
    ----------------
    default_value : 0
        Default value returned if a clone not in the class is accessed.

    Attributes
    ----------
    clones_and_data : dict
        Dictionary keyed by clones, values are counts for each clone.
    norm : int or float
        Normalization of the counts. Generally this is sum of the counts of all
        clones in the class, however it can be set independently if the
        normalization needs to be different. This value will be automatically
        updated if items are added/deleted or if the counts are updated for
        a given clone.

    Methods
    -------
    get_frequency(clone, add_pseudo_count = False)
        Returns the frequency of the clone in a repertoire defined by the
        clones/counts in the class (self[clone]/self.norm). If pseudocounts are
        added, it will change the default count from 0 to 1/3. (So a the
        frequency returned from a clone not in the class will be
        1/(3*self.norm)).
    set_norm(norm = None)
        Sets self.norm. If no norm is provided it will set the norm to the sum
        of the counts over all clones (i.e. self.norm = sum(self.values()))

    """
    def __init__(self, **kwargs):

        CloneDefinition.__init__(self, **kwargs)
        TCRClonePvalue.__init__(self)
        ClonesAndCounts.__init__(self)

        self.name = ''
        self.filenames = []
        self.full_clone_list = []
        self.pvalues = {}

        for kw, kw_val in kwargs.items():
            if kw in self.__dict__:
                self.__dict__[kw] = kw_val
            elif 'filename' in kw:
                if 'adaptive' in kw:
                    self.load_adaptive_file(kw_val)
                elif '10x_clonotype' in kw:
                    self.load_10x_clonotypes(kw_val)

    def __repr__(self):
        clonedef_str = '{%s}'%(', '.join([kw + '=' + str(kw_val) for kw, kw_val in self.get_clone_def().items()]))
        out_list = ['Name: %s'%(self.name),
                    'Number of clones: %s'%(len(self.clones_and_data)),
                    'Norm: %s'%(self.norm),
                    'Clone definition: %s'%(clonedef_str)]
        
        if len(self.pvalues) >  0:
            pval_str = 'Significant clones: {%s}'%(', '.join(['%s: %s'%(pval_name, len(clones_and_pvals)) for pval_name, clones_and_pvals in self.pvalues.items()]))
            out_list.append(pval_str)
        
        
        
        data_head_str = ', '.join(['%s: %s'%(clone, datum) for clone, datum in list(self.clones_and_data.items())[:10]])
        if len(self.clones_and_data) > 10:
            data_head_str += ', ...'
            
        out_list.append('Data: {%s}'%(data_head_str))
        
        return "TcellRepertoire(%s)"%(',\n'.join(out_list))
#            return "TcellRepertoire(%s)"%('Name: %s, Number of clones: %s, Norm: %s, Clone def: %s%s, data: {%s}'%(self.name, len(self.clones_and_data), self.norm, clonedef_str, pval_str, data_head_str))
        
        # else:
        #     return 'TR(Clones: %i, Norm: %i)'%(len(self.clones_and_data), self.norm)
    
    def get_clones_and_counts(self, clones = None, **kwargs):
        c_clone_def = self.get_clone_def()
        for kw, kw_val in kwargs.items():
            if kw in c_clone_def: c_clone_def[kw] = kw_val
        clone_and_count_dict = ClonesAndCounts()
        if clones is None:
            if kwargs.get('split_clones', False):
                #If a clone has multiple TRB or TRA seqs this will return all possible
                #clone representations for SINGLE TRB/TRA seq clones. We keep the normalization constant.
                c_norm = 0
                for clone in self.full_clone_list:
                    c_norm += clone.count
                    c_clone_reps = clone.split_clone_rep(**c_clone_def)
                    for c_clone_rep in c_clone_reps:
                        clone_and_count_dict[c_clone_rep] += clone.count
                clone_and_count_dict.set_norm(norm = c_norm)
            else:
                for clone in self.full_clone_list:
                    clone_and_count_dict[clone.clone_rep(**c_clone_def)] += clone.count
        else:
            if kwargs.get('split_clones', False):
                #If a clone has multiple TRB or TRA seqs this will return all possible
                #clone representations for SINGLE TRB/TRA seq clones. We keep the normalization constant.
                c_norm = 0
                for clone in self.full_clone_list:
                    c_norm += clone.count
                    c_clone_reps = clone.split_clone_rep(**c_clone_def)
                    for c_clone_rep in c_clone_reps:
                        if c_clone_rep in clones:
                            clone_and_count_dict[c_clone_rep] += clone.count
                clone_and_count_dict.set_norm(norm = c_norm)
            else:
                c_norm = 0
                for clone in self.full_clone_list:
                    c_clone_rep = clone.clone_rep(**c_clone_def)
                    c_norm += clone.count
                    if c_clone_rep in clones:
                        clone_and_count_dict[c_clone_rep] += clone.count
                clone_and_count_dict.set_norm(norm = c_norm)
        return clone_and_count_dict


    def add_comparison_pvalue(self, name, comparison_repertoire, pval_type = 'fisher', **kwargs):

        self_clone_and_counts = self.get_clones_and_counts(**kwargs)

        comparison_clone_and_counts = comparison_repertoire.get_clones_and_counts(**kwargs)

        if 'clone_and_pvals' in kwargs:
            self.pvalues[name] = kwargs['clone_and_pvals'].copy()
        elif pval_type == 'fisher':
            self.pvalues[name] = self.compute_fisher_pvalues(self_clone_and_counts, comparison_clone_and_counts, **kwargs)
        elif pval_type == 'binomial':
            self.pvalues[name] = self.compute_binomial_pvalues(self_clone_and_counts, comparison_clone_and_counts, **kwargs)
        else:
            InputError('pval_type', pval_type)

    def get_significant_clones(self, name, pval_thresh):
#        if name not in self.pvalues:
#            InputError('name', name)
        return self.pvalues[name].get_significant_clones(pval_thresh)

    def load_adaptive_file(self, infile_names, primer_correction = 6, ntseq_index = 0, aaseq_index = 1, count_index = 2, use_single_gene_assign = True):
        if type(infile_names) is str:
            infile_names = [infile_names]
        for infile_name in infile_names:
            #print('Loading adaptive file: %s'%(infile_name))
            with open(infile_name, 'r') as infile:
                all_L = [l.split('\t') for l in infile.read().split('\n') if len(l) > 0]
            self.filenames.append(infile_name)
            c_header = [h.lower().strip() for h in all_L[0]]

            v_gene_index = c_header.index('vmaxresolved')
            j_gene_index = c_header.index('jmaxresolved')

            for line in all_L[1:]:
                try:
                    if len(line[aaseq_index].strip()) == 0 or '*' in line[aaseq_index]: continue

                    nt_read = line[ntseq_index]
                    ntseq = nt_read[:len(nt_read)-primer_correction][-3*len(line[aaseq_index]):]

                    count = float(line[count_index])

                    v_gene = line[v_gene_index].strip()
                    j_gene = line[j_gene_index].strip()

                    max_kwargs = CloneDefinition(use_alleles = True).get_clone_def()
                    max_kwargs['ntseq'] = ntseq

                    for i, v in enumerate(v_gene.split('/')):
                        max_kwargs['v' + str(i)] = v
                        if use_single_gene_assign: break

                    for i, j in enumerate(j_gene.split('/')):
                        max_kwargs['j' + str(i)] = j
                        if use_single_gene_assign: break

                    c_tcr_seq = TCRseq(**max_kwargs)

                    if c_tcr_seq.chain is None: c_tcr_seq.chain = 'TRB'

                    self.full_clone_list.append(TcellClone(**{c_tcr_seq.chain: max_kwargs, 'count': count}))
                    #self.clones_and_counts[c_clone_rep] += count
                except:
                    #print('Sequence load failure: %s, %s, %s'%(ntseq, v_gene, j_gene))
                    pass

        self.set_clones_and_counts_attr()

    def load_10x_clonotypes(self, infile_names, ntseq_index = 4, count_index = 1, id_index = 0):
        if type(infile_names) is str:
            infile_names = [infile_names]
        for infile_name in infile_names:
            #print('Loading 10x clonotypes file: %s'%(infile_name))
            with open(infile_name, 'r') as infile:
                all_L = [l.split(',') for l in infile.read().strip().split('\n') if len(l) > 0]
            self.filenames.append(infile_name)

            for i, l in enumerate(all_L[1:]):
                try:
                    c_seqs = l[ntseq_index].split(';')

                    c_clone = TcellClone(count = float(l[count_index]))
                    for c_seq in c_seqs:
                        c_chain, c_ntseq = tuple(c_seq.split(':'))
                        if c_chain == 'TRA':
                            c_clone.load_TRA(ntseq = c_ntseq)
                        elif c_chain == 'TRB':
                            c_clone.load_TRB(ntseq = c_ntseq)

                    self.full_clone_list.append(c_clone)
                except:
                    print(i, l)

        self.set_clones_and_counts_attr()

    def set_clones_and_counts_attr(self, **kwargs):
        for kw, kw_val in kwargs.items():
            if kw in self.__dict__:
                self.__dict__[kw] = kw_val
        self.clones_and_data = self.get_clones_and_counts().clones_and_data
        self.set_norm()

#%%
class IndividualTcellReps(TCRClonePvalue):
    
    pseudo_count = 1/3
    def __init__(self, data_dir = None, **kwargs):
        self.samples = {}
        self.clones = []
        self.sample_metadata = {}
        if data_dir is not None:
            with open(os.path.join(data_dir, 'metadata.json'), 'r') as metadata_infile:
                self.sample_metadata = json.load(metadata_infile)
            for sample, c_metadata in self.sample_metadata.items():
                #print(sample, c_metadata['filename'])
                self.load_adaptive_sample(sample, os.path.join(data_dir, c_metadata['filename']), **kwargs)
        self.pvalues = {}

    def __getitem__(self, sample):
        return self.samples[sample]
    
    def __repr__(self):
        pval_str = ', '.join(['%s: %s'%(pval_name, len(clones_and_pvals)) for pval_name, clones_and_pvals in self.pvalues.items()])
        sample_str = ', '.join(self.samples.keys())
        out_list = ['Total clones: %s'%(len(self.clones)),
                    'Significant clones: {%s}'%(pval_str),
                    'Samples (n = %s): [%s]'%(len(self.samples), sample_str)]
        return "IndividualTcellReps(%s)"%(',\n'.join(out_list))


    def load_adaptive_sample(self, name, infile_name,  **kwargs):

        self.samples[name] = TcellRepertoire(name = name, adaptive_filename = infile_name, **kwargs)
        self.clones = self[name].clone_union(self.clones)
        
    def add_comparison_pvalue(self, specific_repertoire, baseline_repertoire, name = None, pval_type = 'fisher', **kwargs):
        
        if type(specific_repertoire) == str:
            if specific_repertoire in self.samples:
                specific_repertoire = self.samples[specific_repertoire]
            else:
                InputError('specific_repertoire', specific_repertoire)
                
        if type(baseline_repertoire) == str:
            if baseline_repertoire in self.samples:
                baseline_repertoire = self.samples[baseline_repertoire]
            else:
                InputError('baseline_repertoire', baseline_repertoire)
                
        if name is None:
            name = '%s_vs_%s_%s'%(specific_repertoire.name, baseline_repertoire.name, pval_type)
        
        specific_clones_and_counts = specific_repertoire.get_clones_and_counts(**kwargs)

        baseline_clones_and_counts = baseline_repertoire.get_clones_and_counts(**kwargs)

        if 'clone_and_pvals' in kwargs:
            self.pvalues[name] = kwargs['clone_and_pvals'].copy()
        elif pval_type == 'fisher':
            self.pvalues[name] = self.compute_fisher_pvalues(specific_clones_and_counts, baseline_clones_and_counts, **kwargs)
        elif pval_type == 'binomial':
            self.pvalues[name] = self.compute_binomial_pvalues(specific_clones_and_counts, baseline_clones_and_counts, **kwargs)
        else:
            InputError('pval_type', pval_type)