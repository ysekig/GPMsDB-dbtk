#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'GPL3.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import os
import sys

class DefaultValues():
    try:
        GENERIC_PATH = os.environ['GPMsDB_PATH']
    except KeyError:
        print('  ERROR ')
        print("The 'GPMsDB_PATH' environment variable is not defined.")
        print('Please set this variable to your reference data package.' + '\n')
        sys.exit(1)

    NO_THREAD = 4           #number of default threads

    GPMsDB_PATH = GENERIC_PATH

    HMMER_TABLE_OUT = 'hmmer.analyze.txt'
    HMMER_OUT = 'hmmer.analyze.ali.txt'
    MARKER_FILE = os.path.join(GPMsDB_PATH, 'hmm', 'ribosomal.hmm')

    PRODIGAL_AA = 'genes.faa'
    PRODIGAL_NT = 'genes.fna'
    PRODIGAL_GFF = 'genes.gff'

    E_VAL = 1e-10
    LENGTH = 0.7
    PSEUDOGENE_LENGTH = 0.3

    MARKER_GENE_STATS = 'peak_list_genomes.tsv'
    PFAM_CLAN_FILE = os.path.join(GPMsDB_PATH, 'hmm', 'ribosomal.hmm')

    CUSTOM_LIST_R = os.path.join(GPMsDB_PATH, 'custom', 'custom_ribosomals.db')
    CUSTOM_LIST_O = os.path.join(GPMsDB_PATH, 'custom', 'custom_others.db')
    CUSTOM_LIST_GENES = os.path.join(GPMsDB_PATH, 'custom', 'custom_genes.db')
    CUSTOM_LIST_NAME = os.path.join(GPMsDB_PATH, 'custom', 'custom_names.db')
    CUSTOM_LIST_TAX = os.path.join(GPMsDB_PATH, 'custom', 'custom_taxonomy.db')
    
