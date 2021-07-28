##################################
#                                #
# Last modified 2021/01/14       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import os
from sets import Set
import Levenshtein
import numpy as np

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','X':'X','a':'t','t':'a','g':'c','c':'g','n':'n','x':'x','R':'R','r':'r','M':'M','m':'m','Y':'Y','y':'y','S':'S','s':'s','K':'K','k':'k','W':'W','w':'w'}
    sequence=''
    for j in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-j-1]]
    return sequence

def run():

    if len(sys.argv) < 12:
        print 'usage: python %s BC1file fieldID pos1 lenBC1 BC2file fieldID2 pos2 lenBC2 BC3file fieldID3 pos3 lenBC3 [-BCedit N] [-revcompBC]' % sys.argv[0]
        print '\t stream the output of PEFastqToTabDelimited.py, then capture the output of this script with PEFastqToTabDelimited-reverse.py'
        print '\t if the barcode sequence are reverse complemented in the read IDs, use the [-revcompBC] option'
        print '\t the default [-BCedit] edit distance value is 1'
        print '\t it is expected that the i7 indexes look like this (dashed added for separation of the structures):'
        print '\t TCGGACGATCATGGG-TTGTTCGG-CAAGTATGCAGCGCGCTCAAGCACGTGGAT-TCTACGAC-AGTCGTACGCCGATGCGAAACATCGGCCACTCCACCAT+GTGTAGAT'
        print '\t specify pos1, pos2, pos2 to indicate the position of the indexes in the barcode array; all values are 0-based'
        sys.exit(1)

    BC1file = sys.argv[1]
    fieldID1 = int(sys.argv[2])
    pos1 = int(sys.argv[3])
    lenBC1 = int(sys.argv[4])
    BC2file = sys.argv[5]
    fieldID2 = int(sys.argv[6])
    pos2 = int(sys.argv[7])
    lenBC2 = int(sys.argv[4])
    BC3file = sys.argv[9]
    fieldID3 = int(sys.argv[10])
    pos3 = int(sys.argv[11])
    lenBC3 = int(sys.argv[12])

    BCedit = 1
    if '-BCedit' in sys.argv:
        BCedit = int(sys.argv[sys.argv.index('-BCedit') + 1])
#        print 'will used a barcoded edit distance of', BCedit

    doRevComp = False
    if '-revcompBC' in sys.argv:
        doRevComp = True
#        print 'will use reverse complemented barcodes'

    BCDict = {}
    BCDict[1] = {}
    BCDict[2] = {}
    BCDict[3] = {}

    lineslist = open(BC1file)
    for line in lineslist:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        BC = fields[fieldID1]
        if doRevComp:
            BC = getReverseComplement(BC)
        BCDict[1][BC] = 1

    lineslist = open(BC2file)
    for line in lineslist:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        BC = fields[fieldID1]
        if doRevComp:
            BC = getReverseComplement(BC)
        BCDict[2][BC] = 1

    lineslist = open(BC3file)
    for line in lineslist:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        BC = fields[fieldID1]
        if doRevComp:
            BC = getReverseComplement(BC)
        BCDict[3][BC] = 1

#    j = 0
    lineslist = sys.stdin
    for line in lineslist:
#        j+=1
#        if j % 1000000 == 0:
#            print j, 'lines processed'
        fields = line.strip().split('\t')
#        barcode7 = fields[0].strip().split(' ')[1].split(':')[-1].split('+')[0]
        barcode7 = fields[0].strip().split(':')[-1].split('+')[0]
        BC1seq = barcode7[pos1:pos1+lenBC1]
        BC2seq = barcode7[pos2:pos2+lenBC2]
        BC3seq = barcode7[pos3:pos3+lenBC3]

        if BCDict[1].has_key(BC1seq):
            BC1 = BC1seq
        else:
            EDist = lenBC1
            NearestRTIdx = []
            for BCindex in BCDict[1].keys():
                LDist = Levenshtein.distance(BC1seq,BCindex)
                if LDist <= BCedit: 
                    if LDist < EDist:
                        EDist = LDist
                        NearestRTIdx = [BCindex]
                    if LDist == EDist:
                        NearestRTIdx.append(BCindex)
            if len(NearestRTIdx) == 0:
                BC1 = 'nan'
            elif len(NearestRTIdx) == 1:
                BC1 = NearestRTIdx[0]
            else:
                BC1 = 'nan'

        if BCDict[1].has_key(BC2seq):
            BC2 = BC2seq
        else:
            EDist = lenBC2
            NearestRTIdx = []
            for BCindex in BCDict[2].keys():
                LDist = Levenshtein.distance(BC2seq,BCindex)
                if LDist <= BCedit: 
                    if LDist < EDist:
                        EDist = LDist
                        NearestRTIdx = [BCindex]
                    if LDist == EDist:
                        NearestRTIdx.append(BCindex)
            if len(NearestRTIdx) == 0:
                BC2 = 'nan'
            elif len(NearestRTIdx) == 1:
                BC2 = NearestRTIdx[0]
            else:
                BC2 = 'nan'

        if BCDict[3].has_key(BC3seq):
            BC3 = BC3seq
        else:
            EDist = lenBC3
            NearestRTIdx = []
            for BCindex in BCDict[3].keys():
                LDist = Levenshtein.distance(BC3seq,BCindex)
                if LDist <= BCedit: 
                    if LDist < EDist:
                        EDist = LDist
                        NearestRTIdx = [BCindex]
                    if LDist == EDist:
                        NearestRTIdx.append(BCindex)
            if len(NearestRTIdx) == 0:
                BC3 = 'nan'
            elif len(NearestRTIdx) == 1:
                BC3 = NearestRTIdx[0]
            else:
                BC3 = 'nan'

        newbarcode7 = '[' + BC1 + '+' + BC2 + '+' + BC3 + ']'
        newID = fields[0].split(' ')[0] + ':::' + newbarcode7
        print newID + '\t' + fields[1] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4]
            
run()
