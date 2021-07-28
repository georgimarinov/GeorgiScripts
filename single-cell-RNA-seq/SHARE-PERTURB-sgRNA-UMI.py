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

    if len(sys.argv) < 12:
        print 'usage: python %s sgRNA_master_file sgRNA_fieldID sgRNA_label_fieldID sgRNAlen sgRNApos read1|read2 for|rev UMIlen UMIpos read1|read2 for|rev outfilename [-UMIedit N] [-sgRNAedit N]' % sys.argv[0]
        print '\t the script assumes that cell barcodes have already been annotated with SHARE-seq-barcode-annotate.py or SHARE-seq-barcode-annotate-UG.py'
        print '\t the default [-sgRNAedit] edit distance value is 1'
        print '\t the default [-UMIedit] edit distance value is 1'
        print '\t stream the output of PEFastqToTabDelimited.py, then capture the output of this script with PEFastqToTabDelimited-reverse.py'
        sys.exit(1)

    sgRNAedit = 1
    if '-sgRNAedit' in sys.argv:
        sgRNAedit = int(sys.argv[sys.argv.index('-sgRNAedit') + 1])
        print 'will used a sgRNA edit distance of', sgRNAedit

    UMIedit = 1
    if '-UMIedit' in sys.argv:
        UMIedit = int(sys.argv[sys.argv.index('-UMIedit') + 1])
        print 'will used a UMI edit distance of', UMIedit

    sgRNAs = sys.argv[1]
    sgRNAfieldID = int(sys.argv[2])
    sgRNAlabelfieldID = int(sys.argv[3])
    sgRNAlen = int(sys.argv[4])
    sgRNApos = int(sys.argv[5])
    sgRNAread = sys.argv[6]
    sgRNAstrand = sys.argv[7]
    UMIlen = int(sys.argv[8])
    UMIpos = int(sys.argv[9])
    UMIread = sys.argv[10]
    UMIstrand = sys.argv[11]
    outfilename = sys.argv[12]

    sgRNADict = {}

    lineslist = open(sgRNAs)
    for line in lineslist:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        sgRNA = fields[sgRNAfieldID]
        label = fields[sgRNAlabelfieldID]
        sgRNADict[sgRNA] = label

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
        if sgRNAread == 'read2':
           sgRNAseq = read2[sgRNApos:sgRNApos + sgRNAlen]
        if sgRNAread == 'read1':
           sgRNAseq = read1[sgRNApos:sgRNApos + sgRNAlen]
        if sgRNAstrand == 'rev':
           sgRNAseq = getReverseComplement(sgRNAseq)

        if UMIseq.count('N') > 1:
            continue

        if sgRNADict.has_key(sgRNAseq):
            sgRNA = sgRNAseq
        else:
            EDist = sgRNAedit + 1
            Nearest = []
            for sgRNA in sgRNADict.keys():
                LDist = Levenshtein.distance(sgRNAseq,sgRNA)
                if LDist <= sgRNAedit: 
                    if LDist < EDist:
                        EDist = LDist
                        Nearest = [sgRNA]
                    if LDist == EDist:
                        if sgRNA not in Nearest:
                            Nearest.append(sgRNA)
            if len(Nearest) == 0:
                sgRNA = 'nan'
            elif len(Nearest) == 1:
                sgRNA = Nearest[0]
            else:
                sgRNA = 'ambiguous'
        if sgRNA == 'nan':
            continue

        if BCDict.has_key(BC):
            pass
        else:
            BCDict[BC] = {}
        if BCDict[BC].has_key(sgRNA):
            pass
        else:
            BCDict[BC][sgRNA] = {}
        
        if BCDict[BC][sgRNA].has_key(UMIseq):
            BCDict[BC][sgRNA][UMIseq] += 1
        else:
            EDist = UMIedit + 1
            Nearest = []
            for UMI in BCDict[BC][sgRNA].keys():
                LDist = Levenshtein.distance(UMIseq,UMI)
                if LDist <= UMIedit: 
                    Nearest.append(UMI)
                    break
            if len(Nearest) == 0:
                if UMIseq.count('N') == 0:
                    BCDict[BC][sgRNA][UMIseq] = 1
            else:
                continue

    outfile = open(outfilename, 'w')
    outline = '#barcode\tsgRNA\tlabel\tUMIs\treads'
    outfile.write(outline + '\n')

    barcodes = BCDict.keys()
    barcodes.sort()

    print len(barcodes)

    for BC in barcodes:
        (BC1,BC2,BC3) = BC
        for sgRNA in BCDict[BC].keys():
            outline = BC1 + '+' + BC2 + '+' + BC3
            if sgRNA != 'ambiguous':
                outline = outline + '\t' + sgRNA + '\t' + sgRNADict[sgRNA] + '\t' + str(len(BCDict[BC][sgRNA].keys()))
                readcounts = 0
                for UMIseq in BCDict[BC][sgRNA].keys():
                    readcounts += BCDict[BC][sgRNA][UMIseq]
                outline = outline + '\t' + str(readcounts)
            else:
                outline = outline + '\t' + sgRNA + '\t' + 'ambiguous' + '\t' + str(len(BCDict[BC][sgRNA].keys())) + '\tnan'
            outfile.write(outline + '\n')
            
run()
