##################################
#                                #
# Last modified 2021/04/13       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import pysam
import os
from sets import Set
import Levenshtein
import numpy as np

# FLAG field meaning
# 0x0001 1 the read is paired in sequencing, no matter whether it is mapped in a pair
# 0x0002 2 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
# 0x0004 4 the query sequence itself is unmapped
# 0x0008 8 the mate is unmapped 1
# 0x0010 16 strand of the query (0 for forward; 1 for reverse strand)
# 0x0020 32 strand of the mate 1
# 0x0040 64 the read is the first read in a pair 1,2
# 0x0080 128 the read is the second read in a pair 1,2
# 0x0100 256 the alignment is not primary (a read having split hits may have multiple primary alignment records)
# 0x0200 512 the read fails platform/vendor quality checks
# 0x0400 1024 the read is either a PCR duplicate or an optical duplicate

def FLAG(FLAG):

    Numbers = [0,1,2,4,8,16,32,64,128,256,512,1024]

    FLAGList=[]

    MaxNumberList=[]
    for i in Numbers:
        if i <= FLAG:
            MaxNumberList.append(i)

    Residual=FLAG
    maxPos = len(MaxNumberList)-1

    while Residual > 0:
        if MaxNumberList[maxPos] <= Residual:
            Residual = Residual - MaxNumberList[maxPos]
            FLAGList.append(MaxNumberList[maxPos])
            maxPos-=1
        else:
            maxPos-=1
  
    return FLAGList

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s BC_dedup_BAMfilename chrom.sizes TSS_bed chrFieldID posFieldID radius window outfilename' % sys.argv[0]
        print '\tNote: '
        sys.exit(1)

    BAM = sys.argv[1]
    chrominfo = sys.argv[2]
    chromInfoList = []
    chromInfoDict = {}
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))
        chromInfoDict[chr] = end
    TSSbed = sys.argv[3]
    chrFieldID = int(sys.argv[4])
    posFieldID = int(sys.argv[5])
    tssR = int(sys.argv[6])
    tssW = int(sys.argv[7])
    outputfilename = sys.argv[8]

    print 'Parsing GTF'

    TSSDict = {}
    TSSFlankDict = {}
    if TSSbed.endswith('.gz'):
        listoflines = gzip.open(TSSbed)
    elif TSSbed.endswith('.bgz'):
        listoflines = gzip.open(TSSbed)
    else:
        listoflines = open(TSSbed)
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields = line.replace('\x00','').strip().split('\t')
        chr = fields[chrFieldID]
        TSS = int(fields[posFieldID])
        if TSSDict.has_key(chr):
            pass
        else:
            TSSDict[chr] = {}
            TSSFlankDict[chr] = {}
        for i in range(TSS - tssW/2,TSS + tssW/2):
            TSSDict[chr][i] = 0
        for i in range(TSS - tssR - tssW/2,TSS - tssR + tssW/2):
            TSSFlankDict[chr][i] = 0
        for i in range(TSS + tssR - tssW/2,TSS + tssR + tssW/2):
            TSSFlankDict[chr][i] = 0

    print 'Finished parsing GTF'

    BCDict = {}

    samfile = pysam.Samfile(BAM, "rb" )

    BCFragmentDict = {}
    RN = 0.0
    for (chr,start,end) in chromInfoList:
        print (chr,start,end)
        try:
            for alignedread in samfile.fetch(chr, start, end):
                fields = str(alignedread).split('\t')
                break
        except:
            print chr, start, end, 'not found in BAM file'
            continue
        for alignedread in samfile.fetch(chr, start, end):
            RN+=1
            if RN % 1000000 == 0:
                print str(RN/1000000) + 'M alignments processed'
            fields = str(alignedread).split('\t')
            [BC1,BC2,BC3] = fields[0].split(':::[')[-1].split(']')[0].split('+')
            BC = (BC1,BC2,BC3)
            pos = alignedread.pos
            if alignedread.is_reverse:
                pos = pos + len(alignedread.query)
            if BCFragmentDict.has_key(BC):
                pass
            else:
                BCFragmentDict[BC] = {}
                BCFragmentDict[BC]['fragments'] = 0
                BCFragmentDict[BC]['TSS'] = 0
                BCFragmentDict[BC]['flanks'] = 0
            BCFragmentDict[BC]['fragments'] += 1
            if TSSDict.has_key(chr):
                if TSSDict[chr].has_key(pos):
                    BCFragmentDict[BC]['TSS'] += 1
                if TSSFlankDict[chr].has_key(pos):
                    BCFragmentDict[BC]['flanks'] += 1

    outfile = open(outputfilename, 'w')
    outline = '#BC\tNumber_fragments\tTSS_reads\tFlank_reads\tTSS_ratio'
    outfile.write(outline + '\n')

    BCs = BCFragmentDict.keys()
    BCs.sort()

    for BC in BCs:
        (BC1,BC2,BC3) = BC
#        print BC1,BC2,BC3, BCFragmentDict[BC]['fragments'], BCFragmentDict[BC]['flanks']
        TSS = BCFragmentDict[BC]['TSS']
        if BCFragmentDict[BC]['flanks'] > 0:
            Sides = BCFragmentDict[BC]['flanks']/2.0
#            print BC1,BC2,BC3, BCFragmentDict[BC]['fragments'], BCFragmentDict[BC]['flanks'], TSS/Sides
            outline = BC1 + '+' + BC2 + '+' + BC3 + '\t' + str(BCFragmentDict[BC]['fragments']/2.0) + '\t' + str(BCFragmentDict[BC]['TSS']) + '\t' + str(BCFragmentDict[BC]['flanks']) + '\t' + str(TSS/Sides)
        else:
#            print BC1,BC2,BC3, BCFragmentDict[BC]['fragments'], BCFragmentDict[BC]['flanks'], 'nan'
            outline = BC1 + '+' + BC2 + '+' + BC3 + '\t' + str(BCFragmentDict[BC]['fragments']/2.0) + '\t' + str(BCFragmentDict[BC]['TSS']) + '\t' + str(BCFragmentDict[BC]['flanks']) + '\t' + 'nan'
        outfile.write(outline + '\n')

    outfile.close()
            
run()
