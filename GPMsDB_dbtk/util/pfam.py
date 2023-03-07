#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'GPL3.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'


from collections import defaultdict

from GPMsDB_dbtk.common import checkFileExists


class PFAM(object):
    def __init__(self, pfamClanFile):
        self.pfamClanFile = pfamClanFile
        self.idToAcc = {}  
        self.clan = {}  
        self.nested = {} 

    def __readClansAndNesting(self):
        checkFileExists(self.pfamClanFile)

        idNested = defaultdict(list)
        for line in open(self.pfamClanFile):
            if '#=GF ID' in line:
                ID = line.split()[2].strip()
            elif '#=GF AC' in line:
                pfamAcc = line.split()[2].strip()
                pfamAcc = pfamAcc[0:pfamAcc.rfind('.')]
                self.idToAcc[ID] = pfamAcc
            elif '#=GF CL' in line:
                clanId = line.split()[2].strip()
                self.clan[pfamAcc] = clanId
            elif '#=GF NE' in line:
                nestedId = line.split()[2].strip()
                idNested[nestedId].append(ID)
                idNested[ID].append(nestedId)

        for ID, nested in idNested.items():
            pfamAcc = self.idToAcc[ID]
            self.nested[pfamAcc] = set([self.idToAcc[x] for x in nested])

    def pfamIdToClanId(self):
        checkFileExists(self.pfamClanFile)

        d = {}
        for line in open(self.pfamClanFile):
            if '#=GF AC' in line:
                pfamAcc = line.split()[2].strip()
            elif '#=GF CL' in line:
                clanId = line.split()[2].strip()
                d[pfamAcc] = clanId

        return d

    def genesInClan(self):
        checkFileExists(self.pfamClanFile)

        d = defaultdict(set)
        for line in open(self.pfamClanFile):
            if '#=GF AC' in line:
                pfamAcc = line.split()[2].strip()
            elif '#=GF CL' in line:
                clanId = line.split()[2].strip()
                d[clanId].update([pfamAcc])

        return d

    def filterHitsFromSameClan(self, markerHits):
        if len(self.clan) == 0:
            self.__readClansAndNesting()

        filteredMarkers = defaultdict(list)
        hitsToORFs = defaultdict(list)
        for markerId, hits in markerHits.items():
            if markerId.startswith('PF'):
                for hit in hits:
                    hitsToORFs[hit.target_name].append(hit)
            else:
                filteredMarkers[markerId] = hits

        for target_name, hits in hitsToORFs.items():
            hits.sort(key=lambda x: (x.full_e_value, x.i_evalue))

            filtered = set()
            for i in range(0, len(hits)):
                if i in filtered:
                    continue

                pfamIdI = hits[i].query_accession
                pfamIdI = pfamIdI[0:pfamIdI.rfind('.')]
                clanI = self.clan.get(pfamIdI, None)
                startI = hits[i].ali_from
                endI = hits[i].ali_to

                for j in range(i + 1, len(hits)):
                    if j in filtered:
                        continue

                    pfamIdJ = hits[j].query_accession
                    pfamIdJ = pfamIdJ[0:pfamIdJ.rfind('.')]
                    clanJ = self.clan.get(pfamIdJ, None)
                    startJ = hits[j].ali_from
                    endJ = hits[j].ali_to

                    if pfamIdI != None and pfamIdJ != None and clanI == clanJ:
                        if (startI <= startJ and endI > startJ) or (startJ <= startI and endJ > startI):
                            if not (pfamIdI in self.nested and pfamIdJ in self.nested[pfamIdI]):
                                filtered.add(j)

            for i in range(0, len(hits)):
                if i in filtered:
                    continue

                filteredMarkers[hits[i].query_accession].append(hits[i])

        return filteredMarkers

    def genesInSameClan(self, genes):
        pfamIdToClanId = self.pfamIdToClanId()

        clans = set()
        for gene in genes:
            clanId = pfamIdToClanId.get(gene, None)
            if clanId != None:
                clans.add(clanId)

        genesInClan = self.genesInClan()

        allGenesInClans = set()
        for clan in clans:
            allGenesInClans.update(genesInClan[clan])

        return allGenesInClans - genes
