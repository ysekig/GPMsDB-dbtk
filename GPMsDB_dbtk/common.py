#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'CC BY-NC-SA 4.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import os
import errno
import sys
import logging
import time
import ntpath

import GPMsDB_dbtk
from GPMsDB_dbtk.defaultValues import DefaultValues


def checkFileExists(inputFile):
    if not os.path.exists(inputFile):
        logger = logging.getLogger('GPMsDB_tk')
        logger.error('Input file does not exists: ' + inputFile)
        sys.exit(1)


def checkDirExists(inputDir):
    if not os.path.exists(inputDir):
        logger = logging.getLogger('GPMsDB_tk')
        logger.error('Input directory does not exists: ' + inputDir)
        sys.exit(1)


def makeSurePathExists(path):
    if not path:
        return

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logger = logging.getLogger('GPMsDB_tk')
            logger.error('Specified path does not exist: ' + path)
            sys.exit(1)

def genomeIdFromFilename(filename):
    genId = os.path.basename(filename)
    genId = os.path.splitext(genId)[0]

    return genId


class StopWatch():
    def __init__(self, logger):
        self.time_start = time.time()
        self.time_latest = self.time_start
        self.logger = logger

    def clear(self):
        self.time_start = time.time()
        self.time_latest = self.time_start

    def lap(self):
        now = time.time()
        self.time_latest = now
        lap = now - self.time_latest
        total = now - self.time_start
        lap2 = int(lap + 0.5)
        h = lap2 // 3600 
        m = (lap2 - h * 3600) // 60
        s = lap2 - h * 3600 - m * 60
        total2 = int(total + 0.5)
        h2 = total2 // 3600 
        m2 = (total2 - h2 * 3600) // 60
        s2 = total2 - h2 * 3600 - m2 * 60
        self.logger.info(
            f"\n {{lap time: {h:02}:{m:02}:{s:02}, total time: {h2:02}:{m2:02}:{s2:02}}}")


def version():
    versionFile = open(os.path.join(GPMsDB_dbtk.__path__[0], 'VERSION'))
    return versionFile.readline().strip()


def logger_init(logger, output_dir=None, filename="GPMsDB-tk.log", silent=False):
    GPMsDB_tk_logger = logger
    GPMsDB_tk_logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                   datefmt="%Y-%m-%d %H:%M:%S")
    stream_logger = logging.StreamHandler(sys.stdout)
    stream_logger.setFormatter(log_format)
    GPMsDB_tk_logger.addHandler(stream_logger)
    if silent:
        GPMsDB_tk_logger.is_silent = True
        stream_logger.setLevel(logging.ERROR)

    if output_dir != None:
        os.makedirs(output_dir, exist_ok=True)
        timestamp_file_logger = logging.FileHandler(
            os.path.join(output_dir, filename), 'a')
        timestamp_file_logger.setFormatter(log_format)
        GPMsDB_tk_logger.addHandler(timestamp_file_logger)

    GPMsDB_tk_logger.info('%s v%s' % ("GPMsDB-dbtk", version()))
    GPMsDB_tk_logger.info(ntpath.basename(
        sys.argv[0]) + ' ' + ' '.join(sys.argv[1:]))
