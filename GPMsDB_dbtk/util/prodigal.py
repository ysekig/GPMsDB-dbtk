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
import stat
import subprocess
import logging
import shutil
import numpy as np

from biolib.seq_io import read_fasta, write_fasta

from GPMsDB_dbtk.defaultValues import DefaultValues
from GPMsDB_dbtk.common import checkFileExists


class ProdigalError(BaseException):
    pass


class Prodigal():
    def __init__(self, outDir):
        self.logger = logging.getLogger('GPMsDB_tk')
        self.checkForProdigal()
        self.aaGeneFile = os.path.join(outDir, DefaultValues.PRODIGAL_AA)
        self.ntGeneFile = os.path.join(outDir, DefaultValues.PRODIGAL_NT)
        self.gffFile = os.path.join(outDir, DefaultValues.PRODIGAL_GFF)

    def run(self, query, bNucORFs=True):

        prodigal_input = query

        if prodigal_input.endswith('.gz'):
            tmp_dir = tempfile.mkdtemp()
            prodigal_input = os.path.join(tmp_dir, os.path.basename(prodigal_input[0:-3]) + '.fna')
            writeFasta(seqs, prodigal_input)

        seqs = read_fasta(prodigal_input)
        totalBases = 0
        for seqId, seq in seqs.items():
            totalBases += len(seq)

        tableCodingDensity = {}
        for translationTable in [4, 11]:
            aaGeneFile = self.aaGeneFile + '.' + str(translationTable)
            ntGeneFile = self.ntGeneFile + '.' + str(translationTable)
            gffFile = self.gffFile + '.' + str(translationTable)

            if totalBases < 100000:
                procedureStr = 'meta'  
            else:
                procedureStr = 'single'  

            if bNucORFs:
                cmd = ('prodigal -p %s -q -m -f gff -g %d -a %s -d %s -i %s > %s 2> /dev/null' % (procedureStr,
                                                                                                    translationTable,
                                                                                                    aaGeneFile,
                                                                                                    ntGeneFile,
                                                                                                    prodigal_input,
                                                                                                    gffFile))
            else:
                cmd = ('prodigal -p %s -q -m -f gff -g %d -a %s -i %s > %s 2> /dev/null' % (procedureStr,
                                                                                            translationTable,
                                                                                            aaGeneFile,
                                                                                            prodigal_input,
                                                                                            gffFile))

            os.system(cmd)

            if not self.__areORFsCalled(aaGeneFile) and procedureStr == 'single':
                cmd = cmd.replace('-p single', '-p meta')
                os.system(cmd)

            prodigalParser = ProdigalGeneFeatureParser(gffFile)

            codingBases = 0
            for seqId, seq in seqs.items():
                codingBases += prodigalParser.codingBases(seqId)

            if totalBases != 0:
                codingDensity = float(codingBases) / totalBases
            else:
                codingDensity = 0
            tableCodingDensity[translationTable] = codingDensity

        bestTranslationTable = 11
        if (tableCodingDensity[4] - tableCodingDensity[11] > 0.05) and tableCodingDensity[4] > 0.7:
            bestTranslationTable = 4

        shutil.copyfile(self.aaGeneFile + '.' + str(bestTranslationTable), self.aaGeneFile)
        shutil.copyfile(self.gffFile + '.' + str(bestTranslationTable), self.gffFile)
        if bNucORFs:
            shutil.copyfile(self.ntGeneFile + '.' + str(bestTranslationTable), self.ntGeneFile)

        for translationTable in [4, 11]:
            os.remove(self.aaGeneFile + '.' + str(translationTable))
            os.remove(self.gffFile + '.' + str(translationTable))
            if bNucORFs:
                os.remove(self.ntGeneFile + '.' + str(translationTable))

        if prodigal_input.endswith('.gz'):
            shutil.rmtree(tmp_dir)

        return bestTranslationTable

    def __areORFsCalled(self, aaGeneFile):
        return os.path.exists(aaGeneFile) and os.stat(aaGeneFile)[stat.ST_SIZE] != 0

    def areORFsCalled(self, bNucORFs):
        if bNucORFs:
            return os.path.exists(self.ntGeneFile) and os.stat(self.ntGeneFile)[stat.ST_SIZE] != 0

        return os.path.exists(self.aaGeneFile) and os.stat(self.aaGeneFile)[stat.ST_SIZE] != 0

    def checkForProdigal(self):
        try:
            subprocess.call(['prodigal', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("Make sure prodigal is on your system path.")
            sys.exit(1)


class ProdigalFastaParser():
    def __init__(self):
        pass

    def genePositions(self, filename):
        checkFileExists(filename)

        gp = {}
        for line in open(filename):
            if line[0] == '>':
                lineSplit = line[1:].split()

                geneId = lineSplit[0]
                startPos = int(lineSplit[2])
                endPos = int(lineSplit[4])

                gp[geneId] = [startPos, endPos]

        return gp


class ProdigalGeneFeatureParser():
    def __init__(self, filename):
        checkFileExists(filename)

        self.genes = {}
        self.lastCodingBase = {}

        self.__parseGFF(filename)

        self.codingBaseMasks = {}
        for seqId in self.genes:
            self.codingBaseMasks[seqId] = self.__buildCodingBaseMask(seqId)

    def __parseGFF(self, filename):
        self.translationTable = None
        for line in open(filename):
            if line.startswith('# Model Data') and not self.translationTable:
                lineSplit = line.split(';')
                for token in lineSplit:
                    if 'transl_table' in token:
                        self.translationTable = int(token[token.find('=') + 1:])

            if line[0] == '#' or line.strip() == '"':
                continue

            lineSplit = line.split('\t')
            seqId = lineSplit[0]
            if seqId not in self.genes:
                geneCounter = 0
                self.genes[seqId] = {}
                self.lastCodingBase[seqId] = 0

            geneId = seqId + '_' + str(geneCounter)
            geneCounter += 1

            start = int(lineSplit[3])
            end = int(lineSplit[4])

            self.genes[seqId][geneId] = [start, end]
            self.lastCodingBase[seqId] = max(self.lastCodingBase[seqId], end)

    def __buildCodingBaseMask(self, seqId):
        codingBaseMask = np.zeros(self.lastCodingBase[seqId])
        for pos in self.genes[seqId].values():
            codingBaseMask[pos[0]:pos[1] + 1] = 1

        return codingBaseMask

    def codingBases(self, seqId, start=0, end=None):
        if seqId not in self.genes:
            return 0

        if end == None:
            end = self.lastCodingBase[seqId]

        return np.sum(self.codingBaseMasks[seqId][start:end])
