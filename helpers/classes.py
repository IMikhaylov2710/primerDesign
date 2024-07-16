import os
from Bio import SeqIO
import pandas as pd
import primer3
from colorama import Fore
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from DBlogic import getRsByCoord, getByRs

class Variant:

    def __init__(self, coordinate, chromosome, target, **kwargs):

        self.coordinate = coordinate
        self.target = target

        #needs sophisticated rewriting with data type assessment
        if not 'chr' in str(chromosome):
            self.chromosome = 'chr'+str(chromosome)
        else:
            self.chromosome = str(chromosome)
        self.ampliconMinLength = kwargs.get('minLen', 70)
        self.ampliconMaxLength = kwargs.get('maxLen', 200)
        self.delta = kwargs.get('delta', 180)
        self.ampliconStart = self.coordinate-self.delta
        self.ampliconEnd = self.coordinate+self.delta
        self.fwCandidates = []
        self.rvCandidates = []
        self.pairs = []
        self.reducedPairs = []

    def __str__(self):
        return f'{self.chromosome}|{self.coordinate}'

    def printPairs(self):
        if self.reducedPairs:
            print('Primer reduction was performed')
            for pair in self.reducedPairs:
                print(pair[0][0], pair[1][0])
        elif self.pairs:
            print('Primer reduction was not performed')
            for pair in self.pairs:
                print(pair[0][0], pair[1][0], pair[0][1], pair[1][1])
        else:
            print('no pairs selected yet, run Variant.findPrimers() and Variant.getPairs()')

    def getRegion(self, ref):
        for record in SeqIO.parse(f'/home/ivan/REFs/{ref}.fasta', 'fasta'):
            if record.id == self.chromosome:
                self.fwRegion = record.seq[self.ampliconStart:self.coordinate]
                self.rvRegion = record.seq[self.coordinate:self.ampliconEnd]
                return print('region set')
    
    def findPrimers(self, **kwargs):
        tmOptimum = kwargs.get('Tm', 58)
        tmDelta = kwargs.get('TmDelta', 1)
        lengths = kwargs.get('lengths', [18, 19, 20, 21, 22])

        for l in lengths:
            
            for i in range(0, len(self.fwRegion)-l):
                i = int(i)
                candidate = Primer(str(self.fwRegion[i:i+l]).upper(), self.ampliconStart+i, i, self.chromosome)
                candidate.checkAlleleDropout()
                if candidate.isUsable(tmLower = tmOptimum-tmDelta, tmHigher = tmOptimum+tmDelta) and not candidate.HasAlleleDropout:
                    self.fwCandidates.append([candidate, i])

            for i in range(0, len(self.rvRegion)-l):
                i = int(i)
                candidate = Primer(str(self.rvRegion[i:i+l]).upper(), self.coordinate+i, i, self.chromosome, orient='rv')
                candidate.checkAlleleDropout()
                if candidate.isUsable(tmLower = tmOptimum-tmDelta, tmHigher = tmOptimum+tmDelta) and not candidate.HasAlleleDropout:
                    self.rvCandidates.append([candidate, i+len(self.fwRegion)+len(candidate.seq)])


    def getPairs(self):
        self.pairs = []
        for forward in self.fwCandidates:
            for reverse in self.rvCandidates:
                pairTm = primer3.calc_heterodimer_tm(str(forward[0].seq), str(reverse[0].seq))
                if pairTm < 0 and self.ampliconMinLength < reverse[1]-forward[1] < self.ampliconMaxLength:
                    self.pairs.append([forward, reverse])
    
    def reducePairs(self):
        self.reducedPairs = []
        fwToRvs = {}
        for pair in self.pairs:
            try:
                fwToRvs[pair[0][0]].append([pair[0][0], pair[1][0], pair[0][1], pair[1][1], float(pair[1][1])-float(pair[0][1])])
            except:
                fwToRvs[pair[0][0]] = [[pair[0][0], pair[1][0], pair[0][1], pair[1][1], float(pair[1][1])-float(pair[0][1])]]
        for k in fwToRvs:
            print(k)
            optimalCandidate = sorted(fwToRvs[k], key=lambda l: l[-1]-180)[-1]
            print('\t REDUCED TO', optimalCandidate)
        self.reducedPairs.append([[optimalCandidate[0], optimalCandidate[2]], [optimalCandidate[1], optimalCandidate[3]]])
        for rP in self.reducedPairs:
            print(rP[0][0], rP[0][1], rP[1][0], rP[1][1])

    def writeToFasta(self, path):
        toWrite = []
        counter = 0
        self.lengths = {}
        self.targetPairs = {}
        for pair in self.pairs:
            toWrite.append(SeqRecord(str(pair[0][0].seq), id = f'{self.target}{counter}_{pair[0][0].orientation}'))
            self.targetPairs[f'{self.target}{counter}'] = [pair[0][0]]
            self.lengths[f'{self.target}{counter}_{pair[0][0].orientation}'] = len(pair[0][0].seq)
            toWrite.append(SeqRecord(str(pair[1][0].seq), id = f'{self.target}{counter}_{pair[1][0].orientation}'))
            self.targetPairs[f'{self.target}{counter}'].append(pair[1][0])
            self.lengths[f'{self.target}{counter}_{pair[1][0].orientation}'] = len(pair[1][0].seq)
            counter += 1
        SeqIO.write(toWrite, path, 'fasta')
        self.path = path

    def runBlast(self, tsvPath, dbpath='/home/ivan/REFs/hg38.fasta'):
        os.system(f'blastn -query {self.path} -outfmt 6 -word_size 8 -db {dbpath} -out {tsvPath}')
        self.blastOutPath = tsvPath

    def assessPairs(self):
        statistics = {}
        self.chosenPairs = []
        with open(self.blastOutPath, 'r') as handle:
            for i in handle:
                isplit = i.rstrip().split('\t')
                if str(isplit[3]) == str(self.lengths[isplit[0]]) and str(isplit[2]) == '100.000':
                    try:
                        statistics[isplit[0]].append(isplit[1]+'|'+isplit[-4])
                    except:
                        statistics[isplit[0]] = [isplit[1]+'|'+isplit[-4]]
        checkedPairs = {}
        for k in statistics:
            if len(statistics[k]) == 1:
                checkedPairs[k.split('_')[0]] = checkedPairs.get(k.split('_')[0], 0)+1
        for i in checkedPairs:
            if checkedPairs[i] == 2:
                self.chosenPairs.append(i)
        for chosenTarget in self.chosenPairs:
            print(self.targetPairs[chosenTarget][0])
            print(self.targetPairs[chosenTarget][1])
                            
class Primer:

    def __init__(self, seq, ampliconStart, iterator, chromosome, **kwargs):

        self.orientation = kwargs.get('orient', 'fw')
        if self.orientation == 'fw':
            self.seq = seq
        else:
            self.seq = Seq(seq.upper()).reverse_complement()

        self.Tm = primer3.calc_tm(seq)
        self.homodimer = primer3.calc_homodimer_tm(seq)
        self.hairpin = primer3.calc_hairpin_tm(seq)

        if float(self.homodimer) > 0:
            self.hasSecondaryStructures = True
        else:
            self.hasSecondaryStructures = False

        if float(self.hairpin) > self.Tm-25:
            self.hasHairpin = True
        else:
            self.hasHairpin = False
        
        self.primerStart = ampliconStart
        self.primerEnd = ampliconStart+iterator

        self.chromosome = chromosome

    def __str__(self):
        return f'{self.seq}|{self.Tm}|{self.orientation}'
    
    def calculateGC(self, input=None):
        if input == None:
            input = str(self.seq)
            
        gc = 0
        total = 0
        for i in input:
            total += 1
            if i.upper() == 'G' or i.upper() == 'C':
                gc += 1
        
        gcContent = gc/total

        return gcContent

    def isUsable(self, **kwargs):
        threeS = self.seq[-5:]
        if kwargs.get('tmLower', 57) < self.Tm < kwargs.get('tmHigher', 59) and not self.hasSecondaryStructures and self.calculateGC(threeS) < 0.8 and not self.hasHairpin:
            return True
        else:
            return False
        
    def checkAlleleDropout(self, getRsByCoord, conn):
        if getRsByCoord(self.chromosome, self.primerStart, self.primerEnd, conn):
            self.HasAlleleDropout = False
        else:
            self.HasAlleleDropout = True

class ProbeBatch:

    def __init__(self, seq, altSeq, coordinate, fwPrimer, rvPrimer):
        self.seq = seq
        self.altSeq = altSeq
        self.coordinate = coordinate
        self.fwPrimer = fwPrimer
        self.rvPrimer = rvPrimer
        self.probes = []
        self.goodProbes = []
        self.badProbes = []
        self.fwTm = primer3.calc_tm(fwPrimer)
        self.rvTm = primer3.calc_tm(rvPrimer)

    def generateBatches(self, **kwargs):
        self.probes = []
        lens = kwargs.get('lengths', [20, 21, 22, 23, 24, 25, 26])
        for i in range(max(0, int(self.coordinate)-19), min(len(self.seq), int(self.coordinate))):
            i = int(i)
            for l in lens:
                self.probes.append([
                    str(self.seq[i:i+l]).upper(), 
                    str(self.altSeq[i:i+l]).upper()
                ])
        return print(f'{len(self.probes)} probes candidates generated')

    
    def assessQuality(self):
        print('assessing probe quality')

        for element in self.probes:

            [ref, alt] = element

            #ref
            TmR = primer3.calc_tm(ref)
            homodimerR = primer3.calc_homodimer_tm(ref)
            heterodimerFwR = primer3.calc_heterodimer_tm(ref, self.fwPrimer)
            heterodimerRvR = primer3.calc_heterodimer_tm(ref, self.rvPrimer)

            #alt
            TmA = primer3.calc_tm(alt)
            homodimerA = primer3.calc_homodimer_tm(alt)
            heterodimerFwA = primer3.calc_heterodimer_tm(alt, self.fwPrimer)
            heterodimerRvA = primer3.calc_heterodimer_tm(alt, self.rvPrimer)

            if self.fwTm + 9.9 > TmR > self.fwTm + 4.9 \
                and self.rvTm + 9.9 > TmR > self.rvTm + 4.9 \
                and self.fwTm + 9.9 > TmA > self.rvTm + 4.9 \
                and self.rvTm + 9.9 > TmA > self.rvTm + 4.9 \
                and not str(ref.upper()).startswith('G'):
                self.goodProbes.append([ref, alt, self.fwTm, self.rvTm, 
                                        TmR, homodimerR, heterodimerFwR, heterodimerRvR, 
                                        TmA, homodimerA, heterodimerFwA, heterodimerRvA])
            else:
                self.badProbes.append([ref, alt, self.fwTm, self.rvTm, 
                                        TmR, homodimerR, heterodimerFwR, heterodimerRvR, 
                                        TmA, homodimerA, heterodimerFwA, heterodimerRvA])
    
    def selectProbe(self):

        if self.goodProbes:

            print('selecting from good probes')
            sortedFirstIteration = sorted(self.goodProbes, key=lambda lisOne: (abs(lisOne[-4]-67), abs(lisOne[-8]-67)))[:10]
            sortedProbes = sorted(sortedFirstIteration, key=lambda lis: (lis[-3], 
                                                                         lis[-7], 
                                                                         lis[-2], 
                                                                         lis[-1], 
                                                                         lis[-5], 
                                                                         lis[-6]), reverse=True )
            
        elif self.badProbes:

            print('selecting from bad probes')
            sortedFirstIteration = sorted(self.badProbes, key=lambda lisOne: (abs(lisOne[-4]-67), abs(lisOne[-8]-67)))[:10]
            sortedProbes = sorted(sortedFirstIteration, key=lambda lis: (lis[-3], 
                                                                         lis[-7], 
                                                                         lis[-2], 
                                                                         lis[-1], 
                                                                         lis[-5], 
                                                                         lis[-6]), reverse=True )

        else:
            
            return print('no probes selected')
        
        [
            self.resultSeqRef,
            self.resultSeqAlt,
            _, 
            _, 
            self.resultTmR, 
            self.resultHomodimerR, 
            self.resultHeterodimerFwR, 
            self.resultHeterodimerRvR,
            self.resultTmA, 
            self.resultHomodimerA, 
            self.resultHeterodimerFwA, 
            self.resultHeterodimerRvA
        ] = sortedProbes[-1]

        print(self.resultSeqRef, self.resultSeqAlt)

        return sortedProbes[-1]
        
    def printBad(self):
        for bp in self.badProbes:
            print(bp)

    def printGood(self):
        for gp in self.goodProbes:
            print(gp)