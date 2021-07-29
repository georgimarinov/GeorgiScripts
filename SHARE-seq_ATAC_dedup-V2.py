##################################
#                                #
# Last modified 2020/11/30       # 
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

    if len(sys.argv) < 2:
        print 'usage: python %s BAMfilename samtools [-addBC]' % sys.argv[0]
        print '\tThe script will print to stdout, capture it with samtools view to make the new BAM file'
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
    outputfilename = sys.argv[3]

    doAddBC = False
    if '-addBC' in sys.argv:
        doAddBC = True
        print 'will add barcode tags to alignments'

    samfile = pysam.Samfile(BAM, "rb" )
    outfile = pysam.Samfile(outputfilename, "wb", template=samfile)

    FilteredBC = 0
    FilteredDup = 0
    RN = 0.0
    for (chr,start,end) in chromInfoList:
        try:
            for alignedread in samfile.fetch(chr, start, end):
                fields = str(alignedread).split('\t')
                break
        except:
            print chr, start, end, 'not found in BAM file'
            continue
        BCFragmentDict = {}
        print (chr,start,end)
#        if RN > 0:
#            print 'filtered due to duplication:', FilteredDup, FilteredDup/RN
#            print 'filtered due to barcode NaN:', FilteredBC, FilteredBC/RN
#            print 'total alignments:', RN, RN/RN
        for alignedread in samfile.fetch(chr, start, end):
            RN+=1
            if RN % 1000000 == 0:
                print str(RN/1000000) + 'M alignments processed'
            fields = str(alignedread).split('\t')
            [BC1,BC2,BC3] = fields[0].split(':::[')[-1].split(']')[0].split('+')
            BC = (BC1,BC2,BC3)
            if BC1 == 'nan':
                FilteredBC += 1
                continue
            if BC2 == 'nan':
                FilteredBC += 1
                continue
            if BC3 == 'nan':
                FilteredBC += 1
                continue
            if alignedread.is_read1:
                E = 1
            else:
                E = 2
            pos = alignedread.pos
            matepos =  alignedread.pnext
            if BCFragmentDict.has_key(BC):
                pass
            else:
                BCFragmentDict[BC] = {}
            if BCFragmentDict[BC].has_key((E,pos,matepos)):
                FilteredDup += 1
                continue
            else:
                BCFragmentDict[BC][(E,pos,matepos)] = 0
                if doAddBC:
                    BCtag = BC1 + '+' + BC2 + '+' + BC3
                    alignedread.tags += [('BC', BCtag)]
                    outfile.write(alignedread)
                outfile.write(alignedread)

    if RN > 0:
        print 'filtered due to duplication:', FilteredDup, FilteredDup/RN
        print 'filtered due to barcode NaN:', FilteredBC, FilteredBC/RN
        print 'total alignments:', RN, RN/RN

    outfile.close()
            
run()
