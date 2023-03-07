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
import shutil
import multiprocessing as mp
import logging
import uuid
import tempfile

from biolib.external.hmmer import HMMER,HmmModelParser

from GPMsDB_dbtk.util.prodigal import Prodigal
from GPMsDB_dbtk.common import genomeIdFromFilename, makeSurePathExists
from GPMsDB_dbtk.defaultValues import DefaultValues
from GPMsDB_dbtk.mw import Mw


class MarkerGeneFinder():
    def __init__(self, threads):
        self.logger = logging.getLogger('GPMsDB_tk')
        self.totalThreads = threads

    def find(self, genFiles, outDir, tableOut, hmmerOut, markerFile):
        HMMER()

        self.threadsPerSearch = max(1, int(self.totalThreads / len(genFiles)))
        self.logger.info("Identifying genes in %d seqs with %d threads:" % (len(genFiles), self.totalThreads))

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for genFile in genFiles:
            workerQueue.put(genFile)

        for _ in range(self.totalThreads):
            workerQueue.put(None)

        seqIdToModels = mp.Manager().dict()

        try:
            calcProc = [mp.Process(target=self.__processGenome, args=(outDir, tableOut, hmmerOut, markerFile, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
            writeProc = mp.Process(target=self.__reportProcess, args=(len(genFiles), seqIdToModels, writerQueue))

            writeProc.start()

            for p in calcProc:
                p.start()

            for p in calcProc:
                p.join()

            writerQueue.put((None, None))
            writeProc.join()
        except:
            for p in calcProc:
                p.terminate()

            writeProc.terminate()

        d = {}
        for binId in seqIdToModels.keys():
            d[binId] = seqIdToModels[binId]

        return d

    def __processGenome(self, outDir, tableOut, hmmerOut, markerFile, queueIn, queueOut):
        markerSetParser = MarkerSetParser(self.threadsPerSearch)

        while True:
            binFile = queueIn.get(block=True, timeout=None)
            if binFile == None:
                break

            binId = genomeIdFromFilename(binFile)
            binDir = os.path.join(outDir, 'bins', binId)
            makeSurePathExists(binDir)

            prodigal = Prodigal(binDir)
            prodigal.run(binFile)
            aaGeneFile = prodigal.aaGeneFile

            hmmModelFile = markerSetParser.createHmmModelFile(binId, markerFile)

            hmmer = HMMER()
            tableOutPath = os.path.join(binDir, tableOut)
            hmmerOutPath = os.path.join(binDir, hmmerOut)

            keepAlignStr = '--noali'
            hmmer.search(hmmModelFile, aaGeneFile, tableOutPath, hmmerOutPath,
                         '--cpu ' + str(self.threadsPerSearch) + ' --notextw -E 0.1 --domE 0.1 ' + keepAlignStr,
                         False)

            M = Mw()
            ms_dic = M.run(aaGeneFile)

            queueOut.put((binId, hmmModelFile))

    def __reportProcess(self, numGenomes, seqIdToModels, queueIn):
        numProcessedGenomes = 0
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished processing %d of %d (%.2f%%) seqs.' % (numProcessedGenomes, numGenomes, float(numProcessedGenomes) * 100 / numGenomes)
            sys.stderr.write('%s\r' % statusStr)
            sys.stderr.flush()

        while True:
            binId, hmmModelFile = queueIn.get(block=True, timeout=None)
            if binId == None:
                break

            modelParser = HmmModelParser(hmmModelFile)
            models = modelParser.models()

            seqIdToModels[binId] = models

            if os.path.exists(hmmModelFile):
                os.remove(hmmModelFile)

            indexFile = hmmModelFile + '.ssi'
            if os.path.exists(indexFile):
                os.remove(indexFile)

            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedGenomes += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) seqs.' % (numProcessedGenomes, numGenomes, float(numProcessedGenomes) * 100 / numGenomes)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')
            
            
            
class MarkerSetParser():
    def __init__(self, threads=1):
        self.logger = logging.getLogger('GPMsDB_tk')
        self.numThreads = threads

    def createHmmModelFile(self, binId, markerFile):
        hmmModelFile = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        shutil.copyfile(markerFile, hmmModelFile)

        return hmmModelFile

