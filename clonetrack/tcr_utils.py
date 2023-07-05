#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Zachary Sethna
"""

from __future__ import print_function, division
import numpy as np
from scipy.stats import fisher_exact, binom_test


def gene_to_num_str(gene_name, gene_type):
    """Strips excess gene name info to number string.

    Parameters
    ----------
    gene_name : str
        Gene or allele name
    gene_type : char
        Genomic cassette type. (i.e. V, D, or J)

    Returns
    -------
    num_str : str
        Reduced gene or allele name with leading zeros and excess
        characters removed.

    """

    num_str = gene_name.lower().split(gene_type.lower())[-1]
    num_str = '-'.join([g.lstrip('0') for g in num_str.split('-')])
    num_str = '*'.join([g.lstrip('0') for g in num_str.split('*')])

    return gene_type.lower() + num_str

def nt2aa(ntseq):
    """Translate a nucleotide sequence into an amino acid sequence.

    Parameters
    ----------
    ntseq : str
        Nucleotide sequence composed of A, C, G, or T (uppercase or lowercase)

    Returns
    -------
    aaseq : str
        Amino acid sequence

    Example
    --------
    >>> nt2aa('TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC')
    'CAWSVAPDRGGYTF'

    """
    nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}
    aa_dict ='KQE*TPASRRG*ILVLNHDYTPASSRGCILVFKQE*TPASRRGWMLVLNHDYTPASSRGCILVF'

    return ''.join([aa_dict[nt2num[ntseq[i]] + 4*nt2num[ntseq[i+1]] + 16*nt2num[ntseq[i+2]]] for i in range(0, len(ntseq), 3) if i+2 < len(ntseq)])

#%% Custom dictionaries for handling T cell clones
class ClonesAndData(object):
    """Base class for data associated with clones.

    This class will behave very similarly to defaultdict from collections. The
    standard dictionary methods (items, keys, values, get) function
    identically, as do the magic methods of __len__, __setitem__, __iter__,
    __delitem__, and __eq__. Accessing items, __getitem__, behaves similarly to
    that of a default dictionary in that it returns the class default_value if
    the key isn't in the class, however unlike a defaultdict it will not add
    the key and default_value to the class.

    Class Attributes
    ----------------
    default_value : None
        Default value returned if a clone not in the class is accessed.

    Attributes
    ----------
    clones_and_data : dict
        Dictionary keyed by clones, values are the data associated with the
        clone.

    Methods
    -------
    clones()
        Returns list of clones in the class (very similar to self.keys())
    clone_intersection(clones)
        Returns the intersection of self.clones() and clones
    clone_union(clones)
        Returns the union of self.clones() and clones

    """

    default_value = None

    def __init__(self, clones_and_data = {}):
        self.clones_and_data = clones_and_data.copy()

    def __len__(self):
        return len(self.clones_and_data)

    def __getitem__(self, clone):
        return self.clones_and_data.get(clone, self.default_value)

    def __setitem__(self, clone, datum):
        self.clones_and_data[clone] = datum

    def __delitem__(self, clone):
        del self.clones_and_data[clone]

    def __iter__(self):
        return iter(self.clones_and_data)

    def __eq__(self, clones_and_data):
        if type(clones_and_data) != type(self): return False
        return self.clones_and_data == clones_and_data.clones_and_data

    def items(self):
        return self.clones_and_data.items()

    def keys(self):
        return self.clones_and_data.keys()

    def values(self):
        return self.clones_and_data.values()

    def get(self, clone, default = default_value):
        return self.clones_and_data.get(clone, default)

    def clones(self):
        return list(self.clones_and_data.keys())

    def clone_intersection(self, clones):
        try:
            clones = set(clones.keys())
        except AttributeError:
            clones = set(clones)

        return list(clones.intersection(self.clones_and_data))

    def clone_union(self, clones):
        try:
            clones = set(clones.keys())
        except AttributeError:
            clones = set(clones)

        return list(clones.union(self.clones_and_data))
    
    def __repr__(self):
        head_str = ', '.join(['%s: %s'%(clone, datum) for clone, datum in list(self.clones_and_data.items())[:10]])
        if len(self.clones_and_data) > 10:
            head_str += ', ...'
        return "ClonesAndData(%s)"%('Number of clones: %s, Data: {%s}'%(len(self.clones_and_data), head_str))

class ClonesAndCounts(ClonesAndData):
    """Class for cell/read counts associated with a clone.

    Inherits from ClonesAndData and thus will behave similarly to a defaultdict
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


    default_value = 0

    def __init__(self, clones_and_counts = {}):
        ClonesAndData.__init__(self)
        self.clones_and_data = clones_and_counts.copy()
        self.set_norm()

    def __setitem__(self, clone, count):
        if count > 0:
            self.norm += count - self.__getitem__(clone)
            self.clones_and_data[clone] = count
        else:
            self.__delitem__(clone)

    def __delitem__(self, clone):
        self.norm -= self.__getitem__(clone)
        del self.clones_and_data[clone]

    def __add__(self, clones_and_counts):
        tot_clones_and_counts = ClonesAndCounts(self.clones_and_data)
        for clone, count in clones_and_counts.clones_and_data.items():
            tot_clones_and_counts[clone] += count
        return tot_clones_and_counts

    def get(self, clone, default = default_value):
        return self.clones_and_data.get(clone, default)

    def set_norm(self, norm = None):
        """Set normalization constant.

        Sets the attribute norm.

        Parameters
        ----------
        norm : float or None (default is None)
            Either sets the norm to the provided value or to the sum of the
            counts of all clones in the class (i.e. sum(self.values())).
        """
        if norm is None:
            self.norm = sum(self.clones_and_data.values())
        else:
            self.norm = norm

    def get_frequency(self, clone, add_pseudo_count = False):
        """Get a clone's frequency.

        Parameters
        ----------
        clone : str
            Clone whose frequency is to be returned
        add_pseudo_count : bool (default is False)
            If True, adds a pseudocount of 1/3 when computing the frequency of
            a clone not in the class.

        Returns
        -------
        freq : float
            Frequency of the clone in the class.

        """
        return self.get(clone, add_pseudo_count/3)/self.norm
    
    def __repr__(self):
        head_str = ', '.join(['%s: %s'%(clone, datum) for clone, datum in list(self.clones_and_data.items())[:10]])
        if len(self.clones_and_data) > 10:
            head_str += ', ...'
        return "ClonesAndCounts(%s)"%('Number of clones: %s, Norm: %s, Data: {%s}'%(len(self.clones_and_data), self.norm, head_str))


class ClonesAndPvals(ClonesAndData):
    """Class for pvalues associated with a clone.

    Inherits from ClonesAndData and thus will behave similarly to a defaultdict
    with default_value of 1. Most functions will act only on the clones that 
    pass the significance threshold (e.g. len, iter, etc). 
    Non-significant clone pvals are preserved in self.clone_and_data and will
    still be returned if queried by get or __getitem__.

    Class Attributes
    ----------------
    default_value : 1
        Default value returned if a clone not in the class is accessed.

    Attributes
    ----------
    clones_and_data : dict
        Dictionary keyed by clones, values are pvalues for each clone.
    pval_thresh : float
        Significance threshold
    significant_clones_and_pvals : dict
        Dictionary keyed by significant clones, values are pvalues for each
        clone.

    Methods
    -------
    get_significant_clones(pval_thresh)
        Returns a ClonesAndPvals object with only the clones/pvals for clones
        that are more significant than the pval_thresh provided
        (so self[clone] < pval_thresh)

    """

    default_value = 1

    def __init__(self, clones_and_pvals = {}, pval_thresh = 1e-2):
        self.clones_and_data = clones_and_pvals.copy()
        self.pval_thresh = pval_thresh
        
        self.significant_clones_and_pvals = {}
        for clone, pval in self.clones_and_data.items():
            if pval<self.pval_thresh:
                self.significant_clones_and_pvals[clone] = pval
                
    def __len__(self):
        return len(self.significant_clones_and_pvals)
    
    def __setitem__(self, clone, pval):
        self.clones_and_data[clone] = pval
        if pval < self.pval_thresh:
            self.significant_clones_and_pvals[clone] = pval
        else:
            try:
                del self.significant_clones_and_pvals[clone]
            except KeyError:
                pass
        
    def __delitem__(self, clone):
        del self.clones_and_data[clone]
        try:
            del self.significant_clones_and_pvals[clone]
        except KeyError:
            pass
        
    def __iter__(self):
        return iter(self.significant_clones_and_pvals)

    def items(self):
        return self.significant_clones_and_pvals.items()

    def keys(self):
        return self.significant_clones_and_pvals.keys()

    def values(self):
        return self.significant_clones_and_pvals.values()

    def get(self, clone, default = default_value):
        return self.clones_and_data.get(clone, default)

    def get_significant_clones(self, pval_thresh = None):
        """Determines which clones are significant.

        Parameters
        ----------
        pval_thresh : float
            Pvalue threshold to determine the significance of clones. Default
            will be self.pval_thresh

        Returns
        -------
        sig_clones_and_pvals : ClonesAndPvals
            ClonesAndPval object with the clones/pvals for clones that pass
            the significance test based on pval_thresh.

        """
        if pval_thresh is None:
            pval_thresh = self.pval_thresh
        sig_clones_and_pvals = ClonesAndPvals(pval_thresh=pval_thresh)

        for clone, pval in self.clones_and_data.items():
            if pval<pval_thresh:
                sig_clones_and_pvals[clone] = pval

        return sig_clones_and_pvals
    
    
    def clones(self):
        return list(self.significant_clones_and_pvals.keys())

    def clone_intersection(self, clones):
        try:
            clones = set(clones.keys())
        except AttributeError:
            clones = set(clones)

        return list(clones.intersection(self.significant_clones_and_pvals))

    def clone_union(self, clones):
        try:
            clones = set(clones.keys())
        except AttributeError:
            clones = set(clones)

        return list(clones.union(self.significant_clones_and_pvals))

    def __repr__(self):
        head_str = ', '.join(['%s: %s'%(clone, datum) for clone, datum in list(self.significant_clones_and_pvals.items())[:10]])
        if len(self.significant_clones_and_pvals) > 10:
            head_str += ', ...'
        return "ClonesAndPvals(%s)"%('Number of significant clones: %s, Pvalue thresh: %.1e, Significant clones: {%s}'%(len(self.significant_clones_and_pvals), self.pval_thresh, head_str))


#%% Pval class
class TCRClonePvalue(object):
    """Class for computing pvalues for a TCR clone.

    Attributes
    ----------

    Methods
    -------
    compute_fisher_pvalues(specific_clones_and_counts, baseline_clones_and_counts, foldchange_thresh = 1, **kwargs)
        Compute pvalues for fisher exact test.
    compute_binomial_pvalues(specific_clones_and_counts, baseline_clones_and_counts, frac_thresh = 0.5, **kwargs)
        Compute pvalues for binomial test.
    """
    def __init__(self):
        pass
    
    def __repr__(self):
        return "TCRClonePvalue"

    def compute_fisher_pvalues(self, specific_clones_and_counts, baseline_clones_and_counts, foldchange_thresh = 1, multiple_hypothesis_correction = 'bonferroni', **kwargs):
        """Compute pvalues for fisher test.

        Parameters
        ----------
        specific_clones_and_counts : ClonesAndCounts
            ClonesAndCounts object for the specific repertoire.
        baseline_clones_and_counts : ClonesAndCounts
            ClonesAndCounts object for the baseline/comparison repertoire.
        foldchange_thresh : float (default is 1)
            The fisher test will be computed with regard to this foldchange.
            For foldchange_thresh == 1, this will be the standard fisher exact
            test. For foldchange_thresh != 1, we implement this by changing the
            normalization value of the baseline/comparison repertoire to
            effectively change the frequencies of the clones in that repertoire.

        Optional Parameters (**kwargs)
        ------------------------------
        pval_cutoff : float (default is 1)
            Significance cutoff for the clones included in the returned
            clones_and_pvals.
        x_norm : int or float
            Manual normalization for baseline_clones_and_counts. Will still
            be adjusted by foldchange_thresh.
        y_norm : int or float
            Manual normalization for specific_clones_and_counts.


        Returns
        -------
        clones_and_pvals : ClonesAndPvals
            ClonesAndPval object with multiple hypothesis adjusted pvals for
            clones based off of the fisher test.

        """


        y_clones_and_counts = specific_clones_and_counts
        x_clones_and_counts = baseline_clones_and_counts

        pval_cutoff = kwargs.get('pval_cutoff', 1)
        
        if 'x_norm' in kwargs:
            x_norm = int(kwargs['x_norm']/foldchange_thresh)
        else:
            x_norm = int(x_clones_and_counts.norm/foldchange_thresh)
        if 'y_norm' in kwargs:
            y_norm = int(kwargs['y_norm'])
        else:
            y_norm = int(y_clones_and_counts.norm)

        clones = y_clones_and_counts.clone_union(x_clones_and_counts)

        c_count_combo_pval_dict = {}

        c_kwargs = {kw: kw_val for kw, kw_val in kwargs.items() if kw == 'pval_thresh'}
        clones_and_pvals = ClonesAndPvals(**c_kwargs)
        for clone in clones:
            x = x_clones_and_counts[clone]
            y = y_clones_and_counts[clone]

            if x/x_norm > y/y_norm:
                c_count_combo_pval_dict[(x, y)] = 1
            if (x,y) not in c_count_combo_pval_dict:
                if multiple_hypothesis_correction == 'bonferroni':
                    c_count_combo_pval_dict[(x, y)] = np.clip(fisher_exact(np.array([[x, x_norm - x], [y, y_norm - y]]))[1]*len(clones), 0, 1)
                elif multiple_hypothesis_correction is None or multiple_hypothesis_correction == 'None':
                    c_count_combo_pval_dict[(x, y)] = np.clip(fisher_exact(np.array([[x, x_norm - x], [y, y_norm - y]]))[1], 0, 1)

            if c_count_combo_pval_dict[(x, y)] <= pval_cutoff:
                clones_and_pvals[clone] = c_count_combo_pval_dict[(x, y)]

        return clones_and_pvals

    def compute_binomial_pvalues(self, specific_clones_and_counts, baseline_clones_and_counts, frac_thresh = 0.5, multiple_hypothesis_correction = 'bonferroni', **kwargs):
        """Compute pvalues for binomial test.

        Parameters
        ----------
        specific_clones_and_counts : ClonesAndCounts
            ClonesAndCounts object for the specific repertoire.
        baseline_clones_and_counts : ClonesAndCounts
            ClonesAndCounts object for the baseline/comparison repertoire.
        frac_thresh : float (default is 0.5)
            The binomial test will be computed with regard to this relative
            fraction.

        Optional Parameters (**kwargs)
        ------------------------------
        pval_cutoff : float (default is 1)
            Significance cutoff for the clones included in the returned
            clones_and_pvals.

        Returns
        -------
        clones_and_pvals : ClonesAndPvals
            ClonesAndPval object with multiple hypothesis adjusted pvals for
            clones based off of the binomial test.

        """

        y_clones_and_counts = specific_clones_and_counts
        x_clones_and_counts = baseline_clones_and_counts

        pval_cutoff = kwargs.get('pval_cutoff', 1)

        clones = y_clones_and_counts.clone_union(x_clones_and_counts)

        c_count_combo_pval_dict = {}

        c_kwargs = {kw: kw_val for kw, kw_val in kwargs.items() if kw == 'pval_thresh'}
        clones_and_pvals = ClonesAndPvals(**c_kwargs)

        for clone in clones:
            x = x_clones_and_counts[clone]
            y = y_clones_and_counts[clone]

            if y/(x+y) < frac_thresh:
                c_count_combo_pval_dict[(x, y)] = 1
            if (x,y) not in c_count_combo_pval_dict:
                if multiple_hypothesis_correction == 'bonferroni':
                    c_count_combo_pval_dict[(x, y)] = np.clip(binom_test(np.array([y, x]), alternative = 'greater', p = frac_thresh)*len(clones), 0, 1)
                elif multiple_hypothesis_correction is None or multiple_hypothesis_correction == 'None':
                    c_count_combo_pval_dict[(x, y)] = np.clip(binom_test(np.array([y, x]), alternative = 'greater', p = frac_thresh), 0, 1)
            clones_and_pvals[clone] = c_count_combo_pval_dict[(x, y)]

            if c_count_combo_pval_dict[(x, y)] <= pval_cutoff:
                clones_and_pvals[clone] = c_count_combo_pval_dict[(x, y)]

        return clones_and_pvals
