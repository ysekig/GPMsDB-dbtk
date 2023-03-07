#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'CC BY-NC-SA 4.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import sys
import os
import logging
import pickle
import time
import ntpath

from GPMsDB_dbtk.common import (makeSurePathExists,checkDirExists,checkFileExists)
from GPMsDB_dbtk.defaultValues import DefaultValues
from GPMsDB_dbtk.util.resultsParser import ResultsParser
from GPMsDB_dbtk.db import Db
from GPMsDB_dbtk.util.markerGeneFinder import MarkerGeneFinder
from GPMsDB_dbtk.common import StopWatch,logger_init



class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger('GPMsDB_tk')
        self.stopwatch = StopWatch(self.logger)

    def binFiles(self, binFolder, binExtension):
        binFiles = []
        if binFolder is not None:
            all_files = os.listdir(binFolder)
            for f in all_files:
                if f.endswith(binExtension):
                    binFile = os.path.join(binFolder, f)
                    if os.stat(binFile).st_size == 0:
                        self.logger.warning("Skipping genome %s as it has a size of 0 bytes." % f)
                    else:
                        binFiles.append(binFile)

        if not binFiles:
            self.logger.error("No genomes found. Check the extension (-x) used to identify bins.")
            sys.exit(1)

        return sorted(binFiles)

    def genome_wf(self, options):
        logger_init(self.logger, options.out_dir, silent = options.silent)
        self.logger.info('[genome_wf] Generate peak peaks from a set of genome fasta files.')

        genFiles = self.binFiles(options.gen_dir, options.extension)

        makeSurePathExists(options.out_dir)
        checkFileExists(DefaultValues.MARKER_FILE)

        mgf = MarkerGeneFinder(options.threads)
        binIdToModels = mgf.find(genFiles,
                                 options.out_dir,
                                 DefaultValues.HMMER_TABLE_OUT,
                                 DefaultValues.HMMER_OUT,
                                 DefaultValues.MARKER_FILE)

        self.logger.info('[genome_wf] Summarizing genome statistics.')

        checkDirExists(options.out_dir)

        RP = ResultsParser(binIdToModels)
        RP.analyseResults(options.out_dir,
                          DefaultValues.HMMER_TABLE_OUT,
                          bIgnoreThresholds=False,
                          evalueThreshold=DefaultValues.E_VAL,
                          lengthThreshold=DefaultValues.LENGTH,
                          bSkipPseudoGeneCorrection=False,
                          bSkipAdjCorrection=False
                          )

        RP.printSummary(anaFolder=options.out_dir)
        markerGenesFile = RP.cacheResults(options.out_dir)

        self.logger.info('Genome peak lists written to: ' + str(markerGenesFile))

        self.stopwatch.lap()

    def list_db(self, options):
        logger_init(self.logger, None, silent = options.silent)
        self.logger.info('[db_list] List all custom database entries in db')

        run = Db()
        run.list()

        self.stopwatch.lap()

    def update_db(self, options):
        logger_init(self.logger, None, silent = options.silent)
        self.logger.info('[update_db] Add custom genomes into the custom database')

        run = Db()
        run.add(options.file)

        self.stopwatch.lap()

    def remove_genome(self, options):
        logger_init(self.logger, None, silent = options.silent)
        self.logger.info('[remove_genome] Remove custom genomes from the custom database')

        run = Db()
        run.remove(options.accessions)

        self.stopwatch.lap()

    def parse_options(self, options):
        if options.subparser_name == 'data':
            self.update_db(options)
        elif options.subparser_name == 'genome_wf':
            self.genome_wf(options)
        elif options.subparser_name == 'list_db':
            self.list_db(options)
        elif options.subparser_name == 'update_db':
            self.update_db(options)
        elif options.subparser_name == 'remove_genome':
            self.remove_genome(options)
        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        return 0

