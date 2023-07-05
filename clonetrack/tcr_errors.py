#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Zachary Sethna
"""

from __future__ import print_function


class Error(Exception):
    pass


class InputError(Error):
    """Exception for inconsistent input

    Attributes
    ----------
    input_kw : str
    input_arg

    """

    def __init__(self, input_kw, input_arg):
        self.input_kw = input_kw
        self.input_arg = input_arg

class ChainMismatch(Error):
    """Exception for mismatched chains

    Attributes
    ----------
    chainA : str
        TCR chain (i.e. TRA or TRB)
    chainB : str
        TCR chain (i.e. TRA or TRB)

    """

    def __init__(self, chainA, chainB):
        self.chainA = chainA
        self.chainB = chainB

class UnknownGene(Error):
    """Exception for unknown or mischaracterized gene

    Attributes
    ----------
    gene_name : str
        Gene or allele name
    gene_type : char
        Genomic cassette type. (i.e. V, D, or J)

    """

    def __init__(self, gene_name, gene_type = None):
        self.gene_name = gene_name
        self.gene_type = gene_type

class SequenceMismatch(Error):
    """Exception for unknown or mischaracterized gene

    Attributes
    ----------
    gene_name : str
        Gene or allele name
    gene_type : char
        Genomic cassette type. (i.e. V, D, or J)

    """

    def __init__(self, seqA, seqB):
        self.seqA = seqA
        self.seqB = seqB

class PhenotypeMismatch(Error):
    """Exception for mismatched phenotype definitions

    Attributes
    ----------
    phenotypesA : list or int
        list of phenotypes or len(phenotypes)
    phenotypesB : list or int
        list of phenotypes or len(phenotypes)

    """

    def __init__(self, phenotypesA, phenotypesB):
        self.phenotypesA = phenotypesA
        self.phenotypesB = phenotypesB