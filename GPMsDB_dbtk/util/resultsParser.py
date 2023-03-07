#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'GPL3.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'


import sys
import os
import ast
from collections import defaultdict
import logging

from biolib.external.hmmer import HMMERParser

from GPMsDB_dbtk.defaultValues import DefaultValues
from GPMsDB_dbtk.common import checkFileExists
from GPMsDB_dbtk.util.pfam import PFAM


class ResultsParser():
    def __init__(self, binIdToModels):
        self.logger = logging.getLogger('GPMsDB_tk')
        self.results = {}
        self.models = binIdToModels
        self.genes = {}
        self.ribosomals = {}
        self.genesOthers = {}
        self.genesRibosomals = {}

    def analyseResults(self,
                       outDir,
                       hmmTableFile,
                       bIgnoreThresholds,
                       evalueThreshold,
                       lengthThreshold,
                       bSkipPseudoGeneCorrection,
                       bSkipAdjCorrection,
                       ):
        self.parseBinHits(outDir, hmmTableFile, bSkipAdjCorrection, bIgnoreThresholds, evalueThreshold, lengthThreshold, bSkipPseudoGeneCorrection)

    def cacheResults(self, outDir):
        markerGenesFile = self.__writeMarkerGeneStats(outDir)
        return markerGenesFile

    def __writeMarkerGeneStats(self, directory):
        markerGenesFile = os.path.join(directory, DefaultValues.MARKER_GENE_STATS)
        fout = open(markerGenesFile, 'w')
        header = "Genome Id\t# ribosomal peaks\t# other peaks\tribosomal list\tothers list\tname\ttaxonomy\n"
        fout.write(header)
        for binId in sorted(self.results.keys()):
            gene_list = ','.join(self.genesOthers[binId])
            ribo_list = ','.join(self.genesRibosomals[binId])
            fout.write(str(binId) + "\t" + str(len(self.genesRibosomals[binId])) + "\t" +
                          str(len(self.genesOthers[binId])) + "\t" + ribo_list + "\t" + gene_list + "\t\t\n")
        fout.close()

        return markerGenesFile

    def parseBinHits(self, outDir,
                     hmmTableFile,
                     bSkipAdjCorrection=False,
                     bIgnoreThresholds=False,
                     evalueThreshold=DefaultValues.E_VAL,
                     lengthThreshold=DefaultValues.LENGTH,
                     bSkipPseudoGeneCorrection=False):
        if not self.models:
            self.logger.error('Models must be parsed before identifying HMM hits.')
            sys.exit(1)

        self.logger.info('Parsing HMM hits to marker genes:')

        numBinsProcessed = 0
        for binId in self.models:
            self.genes[binId] = {}
            self.ribosomals[binId] = []
            self.genesOthers[binId] = []
            self.genesRibosomals[binId] = []
            geneTableFile = os.path.join(outDir, 'bins', binId, "mw_s.txt")
            checkFileExists(geneTableFile)
            for line in open(geneTableFile):
                if line.rstrip() == "":
                    break
                else:
                    element = line.split("\t")
                    try:
                        self.genes[binId][element[0].rstrip()] = element[1].rstrip()
                    except:
                        continue

            if self.logger.getEffectiveLevel() <= logging.INFO:
                numBinsProcessed += 1
                statusStr = '    Finished parsing hits for %d of %d (%.2f%%) bins.' % (numBinsProcessed, len(self.models), float(numBinsProcessed) * 100 / len(self.models))
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

            resultsManager = ResultsManager(binId, self.models[binId], bIgnoreThresholds, evalueThreshold, lengthThreshold, bSkipPseudoGeneCorrection)
            hmmerTableFile = os.path.join(outDir, 'bins', binId, hmmTableFile)
            self.parseHmmerResults(hmmerTableFile, resultsManager, bSkipAdjCorrection)
            self.results[binId] = resultsManager
            for marker, hitList in resultsManager.markerHits.items():
                for hit in hitList:
                    self.ribosomals[binId].append(hit.target_name)

            for n in self.genes[binId].keys():
                if n in self.ribosomals[binId]:
                    self.genesRibosomals[binId].append(self.genes[binId][n])
                else:
                    self.genesOthers[binId].append(self.genes[binId][n])

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

    def parseHmmerResults(self, fileName, resultsManager, bSkipAdjCorrection):
        try:
            with open(fileName, 'r') as hmmerHandle:
                try:
                    HP = HMMERParser(hmmerHandle)
                except:
                    print("Error opening HMM file: ", fileName)
                    raise

                while True:
                    hit = HP.next()
                    if hit is None:
                        break
                    resultsManager.addHit(hit)

                pfam = PFAM(DefaultValues.PFAM_CLAN_FILE)
                resultsManager.markerHits = pfam.filterHitsFromSameClan(resultsManager.markerHits)

        except IOError as detail:
            sys.stderr.write(str(detail) + "\n")

    def printSummary(self, anaFolder):
        header = "Genome Id\t# ribosomal peaks\t# other peaks"
        self.logger.info(header)

        seqsReported = 0
        for binId in sorted(self.results.keys()):
            gene_list = ','.join(self.genesOthers[binId])
            ribo_list = ','.join(self.genesRibosomals[binId])
            self.logger.info(str(binId) + "\t" + str(len(self.genesRibosomals[binId])) + "\t" +
                   str(len(self.genesOthers[binId])))

            

class ResultsManager():
    def __init__(self, binId, models,
                 bIgnoreThresholds=False,
                 evalueThreshold=DefaultValues.E_VAL,
                 lengthThreshold=DefaultValues.LENGTH,
                 bSkipPseudoGeneCorrection=False,
                 binStats=None):
        self.binId = binId
        self.markerHits = {}
        self.bIgnoreThresholds = bIgnoreThresholds
        self.evalueThreshold = evalueThreshold
        self.lengthThreshold = lengthThreshold
        self.bSkipPseudoGeneCorrection = bSkipPseudoGeneCorrection
        self.models = models

    def vetHit(self, hit):
        model = self.models[hit.query_accession]

        if not self.bSkipPseudoGeneCorrection:
            alignment_length = float(hit.ali_to - hit.ali_from)
            length_perc = alignment_length / float(hit.query_length)
            if length_perc < DefaultValues.PSEUDOGENE_LENGTH:
                return False

        if model.nc != None and not self.bIgnoreThresholds and 'TIGR' in model.acc:
            if model.nc[0] <= hit.full_score and model.nc[1] <= hit.dom_score:
                return True
        elif model.ga != None and not self.bIgnoreThresholds:
            if model.ga[0] <= hit.full_score and model.ga[1] <= hit.dom_score:
                return True
        elif model.tc != None and not self.bIgnoreThresholds:
            if model.tc[0] <= hit.full_score and model.tc[1] <= hit.dom_score:
                return True
        elif model.nc != None and not self.bIgnoreThresholds:
            if model.nc[0] <= hit.full_score and model.nc[1] <= hit.dom_score:
                return True
        else:
            if hit.full_e_value > self.evalueThreshold:
                return False

            alignment_length = float(hit.ali_to - hit.ali_from)
            length_perc = alignment_length / float(hit.query_length)
            if length_perc >= self.lengthThreshold:
                return True

        return False

    def addHit(self, hit):
        if self.vetHit(hit):
            if hit.query_accession in self.markerHits:
                previousHitToORF = None
                for h in self.markerHits[hit.query_accession]:
                    if h.target_name == hit.target_name:
                        previousHitToORF = h
                        break

                if not previousHitToORF:
                    self.markerHits[hit.query_accession].append(hit)
                else:
                    if previousHitToORF.dom_score < hit.dom_score:
                        self.markerHits[hit.query_accession].append(hit)
                        self.markerHits[hit.query_accession].remove(previousHitToORF)

            else:
                self.markerHits[hit.query_accession] = [hit]

    def countUniqueHits(self):
        uniqueHits = 0
        multiCopyHits = 0
        for hits in self.markerHits.values():
            if len(hits) == 1:
                uniqueHits += 1
            elif len(hits) > 1:
                multiCopyHits += 1

        return uniqueHits, multiCopyHits

    def hitsToMarkerGene(self, markerSet):
        ret = dict()
        for marker in markerSet.getMarkerGenes():
            try:
                ret[marker] = len(self.markerHits[marker])
            except KeyError:
                ret[marker] = 0

        return ret

    def printSummary(self, outputFormat, binMarkerSets, bIndividualMarkers, coverageBinProfiles=None, table=None, anaFolder=None):
        if outputFormat == 8:
            markerGenes = binMarkerSets.selectedMarkerSet().getMarkerGenes()

            genesWithMarkers = {}
            for marker, hit_list in self.markerHits.items():
                if marker not in markerGenes:
                    continue

                for hit in hit_list:
                    genesWithMarkers[hit.target_name] = genesWithMarkers.get(hit.target_name, []) + [hit]

            for geneId, hits in genesWithMarkers.items():
                rowStr = self.binId + '\t' + geneId
                for hit in hits:
                    rowStr += '\t' + hit.query_accession + ',' + str(hit.ali_from) + ',' + str(hit.ali_to)
                print(rowStr)

        return 0
