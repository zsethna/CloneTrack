#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Zachary Sethna
"""

from __future__ import print_function, division
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

from tcr_errors import ChainMismatch, UnknownGene, SequenceMismatch
from tcr_utils import nt2aa, gene_to_num_str

class SequenceDefinition(object):
    """Class for TCR sequence type definition

    Determines what level of specificity of sequence definition will be used.
    Also includes magic methods for the operators =, <=, >=. This will indicate
    if a definition is

    Attributes/Keyword Arguments
    ----------------------------
    seq_type : 'ntseq' or 'aaseq'
        Defines whether the CDR3 sequence is defined on the nucleotide sequence
        or amino acid sequence level. (Default: 'ntseq')
    use_genes : bool
        If True the cassette gene information will be included. (Default: True)
    use_alleles : bool
        If True the allele information will be included in the sequence.
        (Default: False)
    use_chain : bool
        If True will specify whether it is a TRB or TRA chain. (Default: False)

    Methods
    -------
    get_seq_def()
        Returns dictionary version of sequence definition (useful for **kwargs).
    """

    def __init__(self, **kwargs):

        self.seq_type = kwargs.get('seq_type', 'ntseq')
        self.use_genes = kwargs.get('use_genes', True)
        self.use_alleles = kwargs.get('use_alleles', False)
        self.use_chain = kwargs.get('use_chain', False)

        
    def get_seq_def(self):
        """Returns the sequence definition dictionary
        
        Returns
        -------
        dict
            Sequence definition dictionary

        """
        return dict(seq_type = self.seq_type,
                    use_genes = self.use_genes,
                    use_alleles = self.use_alleles,
                    use_chain = self.use_chain)
    
    def __repr__(self):
        return "SequenceDefinition(%s)"%(', '.join([kw + '=' + str(kw_val) for kw, kw_val in self.get_seq_def().items()]))
        
        
    def __eq__(self, seq_def):
        if type(seq_def) == SequenceDefinition:
            seq_def = seq_def.get_seq_def()
        return self.get_seq_def() == seq_def

    def __le__(self, seq_def):
        if type(seq_def) == SequenceDefinition:
            seq_def = seq_def.get_seq_def()

        if self.seq_type == 'ntseq' and seq_def['seq_type'] == 'aaseq':
            return False

        return not any([self.__dict__[cond] for cond in ['use_genes', 'use_alleles', 'use_chain'] if not seq_def[cond]])

    def __ge__(self, seq_def):
        if type(seq_def) == SequenceDefinition:
            seq_def = seq_def.get_seq_def()

        if self.seq_type == 'aaseq' and seq_def['seq_type'] == 'ntseq':
            return False

        return all([self.__dict__[cond] for cond in ['use_genes', 'use_alleles', 'use_chain'] if seq_def[cond]])

class CloneDefinition(SequenceDefinition):
    """Class for T cell clone type definition

    Inherits from SequenceDefinition and adds the specification of what
    combination of TRB and TRA chains to use.

    Attributes/Keyword Arguments
    ----------------------------
    seq_type : 'ntseq' or 'aaseq'
        Defines whether the CDR3 sequence is defined on the nucleotide sequence
        or amino acid sequence level. (Default: 'ntseq')
    use_genes : bool
        If True the cassette gene information will be included. (Default: True)
    use_alleles : bool
        If True the allele information will be included in the sequence.
        (Default: False)
    use_chain : bool
        If True will specify whether it is a TRB or TRA chain. (Default: False)
    chains : list
        Specifies which chains (of TRB and TRA) are to be included.
        (Default: ['TRB', 'TRA'])

    Methods
    -------
    get_clone_def()
        Returns dictionary version of clone definition (useful for **kwargs).
    """

    def __init__(self, **kwargs):

        SequenceDefinition.__init__(self)

        self.chains = ['TRB', 'TRA']

        self.use_chain = True

        for kw, kw_val in kwargs.items():
            if kw in self.__dict__: self.__dict__[kw] = kw_val

    def get_clone_def(self):
        """Returns the clone definition dictionary
        
        Returns
        -------
        dict
            Clone definition dictionary
            
        """
        return dict(seq_type = self.seq_type,
                    use_genes = self.use_genes,
                    use_alleles = self.use_alleles,
                    use_chain = self.use_chain,
                    chains = self.chains)

    def __repr__(self):
        return "CloneDefinition(%s)"%(', '.join([kw + '=' + str(kw_val) for kw, kw_val in self.get_clone_def().items()]))

    def __eq__(self, clone_def):

        if type(clone_def) == CloneDefinition:
            clone_def = clone_def.get_clone_def()
        return self.get_clone_def() == clone_def

    def __le__(self, clone_def):
        if type(clone_def) == CloneDefinition:
            clone_def = clone_def.get_clone_def()

        if not set(self.chains).issubset(clone_def['chains']):
            return False

        if self.seq_type == 'ntseq' and clone_def['seq_type'] == 'aaseq':
            return False

        return not any([self.__dict__[cond] for cond in ['use_genes', 'use_alleles', 'use_chain'] if not clone_def[cond]])


    def __ge__(self, clone_def):
        if type(clone_def) == CloneDefinition:
            clone_def = clone_def.get_clone_def()

        if not set(clone_def['chains']).issubset(self.chains):
            return False

        if self.seq_type == 'aaseq' and clone_def['seq_type'] == 'ntseq':
            return False

        return all([self.__dict__[cond] for cond in ['use_genes', 'use_alleles', 'use_chain'] if clone_def[cond]])

class TCRseq(SequenceDefinition):
    """Class for TCR sequences (single chain)

    The class stores all provided sequence information for a TCR sequence
    (single chain -- only the TRB or TRA) and will return a standardized syntax
    string based on the SequenceDefinition that can be used as a hash for the
    sequence. Inherits all of the attributes of SequenceDefinition.
    
    Attributes
    ----------
    ntseq : str
        Nucleotide sequence of TCR CDR3 sequence (including conserved C and F)
    aaseq : str
        Amino acid sequence of TCR CDR3 sequence (including conserved C and F)
    chain : 'TRB' or 'TRA'
        Specifies what chain the sequence corresponds to
    v_genes : list
        List of V genes associated with the sequence
    j_genes : list
        List of J genes associated with the sequence
    v_alleles : list
        List of V alleles associated with the sequence
    j_alleles : list
        List of J alleles associated with the sequence
    
    Methods
    -------
    seq_rep(**kwargs)
        Returns standardized sequence string based on the sequence definition
        (modified by any **kwargs)
    set_chain(chain)
        Sets self.chain
    add_gene(gene)
        Adds genes/alleles to gene/allele lists
    """

    def __init__(self, ntseq = None, aaseq = None, **kwargs):


        SequenceDefinition.__init__(self, **kwargs)

        self.ntseq = ntseq
        if ntseq is not None:
            self.aaseq = nt2aa(ntseq)
            if aaseq is not None and aaseq != self.aaseq:
                raise SequenceMismatch(self.aaseq, aaseq)
        else:
            self.aaseq = aaseq
            self.seq_type = 'aaseq'

        self.chain = None
        self.v_genes = []
        self.j_genes = []
        self.v_alleles = []
        self.j_alleles = []

        for kw, kw_val in kwargs.items():
            if kw.lower() == 'chain':
                self.set_chain(kw_val.upper())
            elif kw.lower()[0] in 'vdj':
                self.add_gene(kw_val, gene_type = kw.lower())

    def set_chain(self, chain):
        """Sets chain of the sequence

        Parameters
        ----------
        chain : 'TRB' or 'TRA'
            Chain of sequence

        Raises
        ------
        ChainMismatch
            Error if chain is already specified and different than input

        Attributes modified
        -------------------
        self.chain

        """

        if self.chain is None:
            self.chain = chain
        elif chain != self.chain:
            raise ChainMismatch(self.chain, chain)

    def add_gene(self, gene_name, gene_type = '', reset = False):
        """Method to add genes/alleles to sequence info.
        
        This function will also infer/check the chain of the sequence based on
        provided gene names.
        
        Parameters
        ----------
        gene_name : str
            Gene/allele name
        gene_type : 'V' or 'J', optional
            Specified type of gene. The default is ''.
        reset : bool, optional
            Resets the gene list to only include provided gene. 
            The default is False.

        Raises
        ------
        UnknownGene
        ChainMismatch

        Attributes Modified
        -------------------
        self.chain
        self.v_genes
        self.v_alleles
        self.j_genes
        self.j_alleles

        """

        if 'A' in gene_name.upper():
            self.set_chain('TRA')
        elif 'B' in gene_name.upper():
            self.set_chain('TRB')

        if gene_type.lower()[0] == 'v' or 'v' in gene_name.lower():
            if reset:
                self.v_genes = [gene_to_num_str(gene_name.lower(), 'v').split('*')[0]]
                self.v_alleles = [gene_to_num_str(gene_name.lower(), 'v')]
            else:
                self.v_genes.append(gene_to_num_str(gene_name.lower(), 'v').split('*')[0])
                self.v_alleles.append(gene_to_num_str(gene_name.lower(), 'v'))
        elif gene_type.lower()[0] == 'j' or 'j' in gene_name.lower():
            if reset:
                self.j_genes = [gene_to_num_str(gene_name.lower(), 'j').split('*')[0]]
                self.j_alleles = [gene_to_num_str(gene_name.lower(), 'j')]
            else:
                self.j_genes.append(gene_to_num_str(gene_name.lower(), 'j').split('*')[0])
                self.j_alleles.append(gene_to_num_str(gene_name.lower(), 'j'))
        else:
            raise UnknownGene(gene_name, gene_type)

    def seq_rep(self, **kwargs):
        """Returns standardized sequence string based on sequence definition.

        Returns
        -------
        str
            Sequence string

        """

        seq_def = self.get_seq_def()
        for kw, kw_val in kwargs.items():
            if kw in seq_def: seq_def[kw] = kw_val


        if seq_def['use_chain'] and self.chain is not None:
            s = [self.chain]
        else:
            s = []

        if seq_def['seq_type'] == 'ntseq' and self.ntseq is not None:
            s.append(self.ntseq)
        elif seq_def['seq_type'] == 'aaseq' and self.aaseq is not None:
            s.append(self.aaseq)

        if seq_def['use_alleles']:
            if len(self.v_alleles) > 0: s.append('/'.join(self.v_alleles))
            if len(self.j_alleles) > 0: s.append('/'.join(self.j_alleles))
        elif seq_def['use_genes']:
            if len(self.v_genes) > 0: s.append('/'.join(self.v_genes))
            if len(self.j_genes) > 0: s.append('/'.join(self.j_genes))



        #self.seq = self.seq_rep()
        return ','.join(s)

    def compare_full_sequence(self, tcr_seq):
        if type(self) != type(tcr_seq):
            raise ValueError
        return self.__dict__ == tcr_seq.__dict__

    def __repr__(self):
        return "TCRseq(%s)"%(self.seq_rep())

#%

class TcellClone(CloneDefinition):
    """Class for TCR clones (single or paired chains)

    The class stores all provided sequence information for TRA and TRB 
    sequences and will return a standardized syntax string based on the 
    CloneDefinition that can be used as a hash for the clone. Inherits all of 
    the attributes of CloneDefinition.
    
    Attributes
    ----------
    TRA : list of TCRseq
        List of all TCRseq associated with the clone's TRA
    TRB : list of TCRseq
        List of all TCRseq associated with clone's TRB
    count : float
        Count associated with the clone. Generally cell count.
    
    Methods
    -------
    clone_rep(**kwargs)
        Returns standardized sequence string based on the clone definition
        (modified by any **kwargs)
    split_clone_rep(**kwargs)
        Returns all possible standardized sequence string with 1 TRA and 1 TRB
        based on the clone definition (modified by any **kwargs)
    load_TRB(**kwargs)
        Load TRB sequence
    load_TRA(**kwargs)
        Load TRA sequence
    """
    
    def __init__(self, **kwargs):

        CloneDefinition.__init__(self, **kwargs)

        self.TRB = []
        self.TRA = []
        self.count = 0

        for kw, kw_val in kwargs.items():
            if kw.upper() == 'TRB':
                self.load_TRB(**kw_val)
            elif kw.upper() == 'TRA':
                self.load_TRA(**kw_val)
            elif kw in self.__dict__:
                self.__dict__[kw] = kw_val
        
        if len(self.TRB) > 0 and len(self.TRA) == 0:
            self.chains = [c for c in self.chains if c!='TRA']
        elif len(self.TRA) > 0 and len(self.TRB) == 0:
            self.chains = [c for c in self.chains if c!='TRB']
            
    def load_TRB(self, **kwargs):
        """Load TRB sequence

        Parameters
        ----------
        **kwargs : keywords
            Input for TCRSeq.__init__ to specify TRB sequence

        Attributes modified
        -------------------
        self.TRB

        """
        seq_def = self.get_clone_def()
        kwargs['chain'] = 'TRB'
        for kw, kw_val in kwargs.items():
            seq_def[kw] = kw_val
        self.TRB.append(TCRseq(**seq_def))

    def load_TRA(self, **kwargs):
        """Load TRA sequence

        Parameters
        ----------
        **kwargs : keywords
            Input for TCRSeq.__init__ to specify TRA sequence

        Attributes modified
        -------------------
        self.TRA

        """
        seq_def = self.get_clone_def()
        kwargs['chain'] = 'TRA'
        for kw, kw_val in kwargs.items():
            seq_def[kw] = kw_val
        self.TRA.append(TCRseq(**seq_def))

    def clone_rep(self, **kwargs):
        """Returns standardized clone string based on clone definition.

        Returns
        -------
        str
            Clone string

        """

        clone_def = self.get_clone_def()
        for kw, kw_val in kwargs.items():
            if kw in clone_def: clone_def[kw] = kw_val

        return ';'.join([c_seq.seq_rep(**clone_def) for c_seq in sum([self.__dict__[c_chain] for c_chain in clone_def['chains']], [])])

    def split_clone_rep(self, **kwargs):
        """Returns all possible clone strings based on clone definition.
        
        If a clone has multiple TRA or TRB sequences this will enumerate the
        possible combinations.

        Returns
        -------
        tuple
            Tuple of clone strings

        """
        clone_def = self.get_clone_def()
        for kw, kw_val in kwargs.items():
            if kw in clone_def: clone_def[kw] = kw_val

        if len(clone_def['chains']) <= 1:
            return tuple([c_seq.seq_rep(**clone_def) for c_seq in self.__dict__[clone_def['chains'][0]]])
        else:
            return tuple([';'.join([c_seq1.seq_rep(**clone_def), c_seq2.seq_rep(**clone_def)]) for c_seq1 in self.__dict__[clone_def['chains'][0]] for c_seq2 in self.__dict__[clone_def['chains'][1]]])

    def __repr__(self):
        return "TcellClone(%s)"%(self.clone_rep())
    