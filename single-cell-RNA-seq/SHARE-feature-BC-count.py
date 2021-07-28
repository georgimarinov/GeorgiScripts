##################################
#                                #
# Last modified 2021/04/29       # 
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
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    sequence=''
    for j in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-j-1]]
    return sequence

def run():

    if len(sys.argv) < 9:
        print 'usage: python %s BClen BCpos read1|read2 for|rev UMIlen UMIpos read1|read2 for|rev outfilename [-UMIedit N] [-BCedit N] [-endogenousPromoter sequence]' % sys.argv[0]
        print '\t the script assumes that cell barcodes have already been annotated with SHARE-seq-barcode-annotate.py or SHARE-seq-barcode-annotate-UG.py'
        print '\t the default [-sgRNAedit] edit distance value is 1'
        print '\t the default [-UMIedit] edit distance value is 1'
        print '\t stream the output of PEFastqToTabDelimited.py, then capture the output of this script with PEFastqToTabDelimited-reverse.py'
        sys.exit(1)

    BCedit = 1
    if '-BCedit' in sys.argv:
        BCedit = int(sys.argv[sys.argv.index('-BCedit') + 1])
        print 'will used a BC edit distance of', BCedit

    UMIedit = 1
    if '-UMIedit' in sys.argv:
        UMIedit = int(sys.argv[sys.argv.index('-UMIedit') + 1])
        print 'will used a UMI edit distance of', UMIedit

    BClen = int(sys.argv[1])
    BCpos = int(sys.argv[2])
    BCread = sys.argv[3]
    BCstrand = sys.argv[4]
    UMIlen = int(sys.argv[5])
    UMIpos = int(sys.argv[6])
    UMIread = sys.argv[7]
    UMIstrand = sys.argv[8]
    outfilename = sys.argv[9]

    BCseqDict = {}

    if '-endogenousPromoter' in sys.argv:
        BCseq = sys.argv[sys.argv.index('-endogenousPromoter') + 1]
        BCseqDict[BCseq] = 1
        print 'will use the endogenous promoter sequence as a pre-annotated barcode', BCseq

    BCDict = {}

    LC = 0

    lineslist = sys.stdin 
    for line in lineslist:
        fields = line.strip().split('\t')
        LC += 1
        if LC % 1000000 == 0:
            print str(LC/1000000) + 'M reads processed'
        [BC1,BC2,BC3] = fields[0].split(':::[')[-1].split(']')[0].split('+')
        BC = (BC1,BC2,BC3)
        if BC1 == 'nan':
            continue
        if BC2 == 'nan':
            continue
        if BC3 == 'nan':
            continue

        read1 = fields[1]
        read2 = fields[3]

        if UMIread == 'read2':
           UMIseq = read2[UMIpos:UMIpos + UMIlen]
        if UMIread == 'read1':
           UMIseq = read1[UMIpos:UMIpos + UMIlen]
        if UMIstrand == 'rev':
           UMIseq = getReverseComplement(UMI)
        if BCread == 'read2':
           BCseq = read2[BCpos:BCpos + BClen]
        if BCread == 'read1':
           BCseq = read1[BCpos:BCpos + BClen]
        if BCstrand == 'rev':
           BCseq = getReverseComplement(BCseq)

        if UMIseq.count('N') > 1:
            continue

        if BCseq.count('N') > 1:
            continue

        if BCseqDict.has_key(BCseq):
            pass
        else:
            EDist = BCedit + 1
            Nearest = []
            for BCseq2 in BCseqDict.keys():
                LDist = Levenshtein.distance(BCseq2,BCseq)
                if LDist <= BCedit: 
                    Nearest.append(BCseq2)
                    break
            if len(Nearest) == 0:
                if BCseq.count('N') == 0:
                    BCseqDict[BCseq] = 1
            elif len(Nearest) == 1:
                BCseq = Nearest[0]
            else:
                continue

        if BCDict.has_key(BC):
            pass
        else:
            BCDict[BC] = {}
        if BCDict[BC].has_key(BCseq):
            pass
        else:
            BCDict[BC][BCseq] = {}
        
        if BCDict[BC][BCseq].has_key(UMIseq):
            BCDict[BC][BCseq][UMIseq] += 1
        else:
            EDist = UMIedit + 1
            Nearest = []
            for UMI in BCDict[BC][BCseq].keys():
                LDist = Levenshtein.distance(UMIseq,UMI)
                if LDist <= UMIedit: 
                    Nearest.append(UMI)
                    break
            if len(Nearest) == 0:
                if UMIseq.count('N') == 0:
                    BCDict[BC][BCseq][UMIseq] = 1
            else:
                continue

    outfile = open(outfilename, 'w')
    outline = '#barcode\tBCseq\tUMIs\treads'
    outfile.write(outline + '\n')

    barcodes = BCDict.keys()
    barcodes.sort()

    print len(barcodes)

    for BC in barcodes:
        (BC1,BC2,BC3) = BC
        for BCseq in BCDict[BC].keys():
            outline = BC1 + '+' + BC2 + '+' + BC3 + '\t' + BCseq + '\t' + str(len(BCDict[BC][BCseq].keys()))
            readcounts = 0
            for UMIseq in BCDict[BC][BCseq].keys():
                readcounts += BCDict[BC][BCseq][UMIseq]
            outline = outline + '\t' + str(readcounts)
            outfile.write(outline + '\n')
            
run()
