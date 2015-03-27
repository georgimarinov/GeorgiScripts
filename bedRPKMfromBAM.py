##################################
#                                #
# Last modified 06/09/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import pysam
import math
from sets import Set
import os

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

    if len(sys.argv) < 5:
        print 'usage: python %s bedfilename chrField BAMfilename chrom.sizes outputfilename [-nomulti] [-RPM] [-stranded +|-] [-readLength min max] [-printSum] [-uniqueBAM] [-mappabilityNormalize mappability.wig readLength] [-noNH samtools]' % sys.argv[0]
        print 'Note: the script will divide multireads by their multiplicity'
        print '\t-printSum option only working together with the RPM option'
        print '\tuse the uniqueBAM option if the BAM file contains only unique alignments; this will save a lot of memory'
        print '\tuse the -mappabilityNormalize option to get mappability normalized RPKMs (it will not do anything to the RPMs; not that a mappability track that goes from 0 to the read length is assumed'
        print '\tuse the -noNH option and supply a path to samtools in order to have the file converted to one that has NH tags'
        print '\tthe stranded option will normalized against all reads, not just reads on the indicated strand'
        sys.exit(1)
    
    bed = sys.argv[1]
    fieldID = int(sys.argv[2])
    SAM = sys.argv[3]
    chromSize = sys.argv[4]
    outfilename = sys.argv[5]

    chromInfoList=[]
    linelist=open(chromSize)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))

    noMulti=False
    if '-nomulti' in sys.argv:
        noMulti=True
        print 'will discard multi-read alignments'

    doReadLength=False
    if '-readLength' in sys.argv:
        doReadLength=True
        minRL = int(sys.argv[sys.argv.index('-readLength')+1])
        maxRL = int(sys.argv[sys.argv.index('-readLength')+2])
        print 'will only consider reads between', minRL, 'and', maxRL, 'bp length'
        ORLL = 0

    doPrintSum=False

    doStranded=False
    if '-stranded' in sys.argv:
        doStranded=True
        thestrand = sys.argv[sys.argv.index('-stranded')+1]
        print 'will only consider', thestrand, 'strand reads'

    doRPM=False
    if '-RPM' in sys.argv:
        doRPM=True
        print 'will output RPMs'
        if '-printSum' in sys.argv:
            doPrintSum=True
            RPMSum=0

    doUniqueBAM = False
    if '-uniqueBAM' in sys.argv:
        print 'will treat all alignments as unique'
        doUniqueBAM = True
        TotalReads = 0
        pass

    samfile = pysam.Samfile(SAM, "rb" )
    try:
        print 'testing for NH tags presence'
        for alignedread in samfile.fetch():
            multiplicity = alignedread.opt('NH')
            print 'file has NH tags'
            break
    except:
        if '-noNH' in sys.argv:
            print 'no NH: tags in BAM file, will replace with a new BAM file with NH tags'
            samtools = sys.argv[sys.argv.index('-noNH')+1]
            BAMpreporcessingScript = sys.argv[0].rpartition('/')[0] + '/bamPreprocessing.py'
            cmd = 'python ' + BAMpreporcessingScript + ' ' + SAM + ' ' + SAM + '.NH'
            os.system(cmd)
            cmd = 'rm ' + SAM
            os.system(cmd)
            cmd = 'mv ' + SAM + '.NH' + ' ' + SAM
            os.system(cmd)
            cmd = samtools + ' index ' + SAM
            os.system(cmd)
        else:
            if doUniqueBAM:
                pass
            else:
                print 'no NH: tags in BAM file, exiting'
                sys.exit(1)

    doMappabilityCorrection = False
    if not doRPM and '-mappabilityNormalize' in sys.argv:
        doMappabilityCorrection = True
        print 'will correct for mappability'
        mappability = sys.argv[sys.argv.index('-mappabilityNormalize')+1]
        readLength = int(sys.argv[sys.argv.index('-mappabilityNormalize')+2])
        WantedDict = {}
        MappabilityRegionDict = {}
        lineslist = open(bed)
        i=0
        print 'inputting regions'
        for line in lineslist:
            if line[0]=='#':
                continue
            i+=1
            if i % 1000 == 0:
                print i, 'regions inputted'
            fields = line.strip().split('\t')
            if len(fields) < fieldID+2:
                continue
            chr = fields[fieldID]
            try:
                left = int(fields[fieldID+1])
                right = int(fields[fieldID+2])
            except:
                print 'problem with region, skipping:', line.strip()
            if left >= right:
                print 'problem with region, skipping:', chr, left, right
                continue
            if MappabilityRegionDict.has_key(chr):
                pass
            else:
                MappabilityRegionDict[chr]={}
                WantedDict[chr]={}
            MappabilityRegionDict[chr][(left,right)]=0
            for j in range(left,right):
                WantedDict[chr][j]=0
        lineslist = open(mappability)
        print 'inputting mappability'
        i=0
        for line in lineslist:
            if line.startswith('#'):
                continue
            i+=1
            if i % 1000000 == 0:
                print str(i/1000000) + 'M lines processed'
            fields = line.strip().split('\t')
            if len(fields) == 1:
                fields = line.strip().split(' ')
            chr = fields[0]
            left = int(fields[1])
            right = int(fields[2])
            score = float(fields[3])
            if WantedDict.has_key(chr):
                pass
            else:
                continue
            for j in range(left,right):
                if WantedDict[chr].has_key(j):
                    WantedDict[chr][j] = score
        print 'calculating mappable fractions'
        for chr in MappabilityRegionDict.keys():
            for (left,right) in MappabilityRegionDict[chr].keys():
                TotalScore = 0.0
                for j in range(left,right):
                    TotalScore += WantedDict[chr][j]
                Score = TotalScore / (right-left)
                MappabilityRegionDict[chr][(left,right)] = Score/readLength
        WantedDict = {}
     
    regionDict={}

    Unique=0
    UniqueSplices=0
    Multi=0
    MultiSplices=0

    if doUniqueBAM and not doReadLength:
        TotalReads = 0
        for chrStats in pysam.idxstats(SAM):
            fields = chrStats.strip().split('\t')
            chr = fields[0]
            reads = int(fields[2])
            if chr != '*':
                TotalReads += reads
        UniqueReads = TotalReads
    else:
        MultiplicityDict={}
        UniqueReads = 0
        i=0
        samfile = pysam.Samfile(SAM, "rb" )
        for (chr,start,end) in chromInfoList:
            try:
                for alignedread in samfile.fetch(chr, start, end):
                    i+=1
                    if i % 5000000 == 0:
                        print str(i/1000000) + 'M alignments processed', chr,start,end
                    fields=str(alignedread).split('\t')
                    if doReadLength:
                        if len(alignedread.seq) > maxRL or len(alignedread.seq) < minRL:
                            ORLL += 1
                            continue
                    if doUniqueBAM:
                        TotalReads+=1
                        continue
                    if alignedread.opt('NH') == 1:
                        UniqueReads += 1
                        continue
                    ID=fields[0]
                    if alignedread.is_read1:
                        ID = ID + '/1'
                    if alignedread.is_read2:
                        ID = ID + '/2'
                    if MultiplicityDict.has_key(ID):
                        MultiplicityDict[ID]+=1
                    else:
                        MultiplicityDict[ID]=1
            except:
                print 'problem with region:', chr, start, end, 'skipping'
        if doReadLength:
            print ORLL, 'alignments outside of read length limits'
        if doUniqueBAM:
            pass
        else:
            TotalReads = UniqueReads + len(MultiplicityDict.keys())

    print TotalReads, UniqueReads

    normalizeBy = TotalReads/1000000.

    outfile = open(outfilename, 'w')

    lineslist = open(bed)
    i=0
    for line in lineslist:
        i+=1
        if i % 10000 == 0:
            print i, 'regions processed'
        if line[0]=='#':
            continue
        fields = line.strip().split('\t')
        if len(fields) < fieldID+2:
            continue
        chr = fields[fieldID]
        try:
            left = int(fields[fieldID+1])
            right = int(fields[fieldID+2])
        except:
            print 'problem with region, skipping:', line.strip()
        if left >= right:
            print 'problem with region, skipping:', chr, left, right
            continue
        reads=0
        try:
            for alignedread in samfile.fetch(chr, left, right):
                fields2=str(alignedread).split('\t')
                if doReadLength:
                    if len(alignedread.seq) > maxRL or len(alignedread.seq) < minRL:
                        continue
                ID=fields2[0]
                if doStranded:
                    if alignedread.is_reverse:
                        s = '-'
                    else:
                        s = '+'
                    if s != thestrand:
                        continue
                if alignedread.is_read1:
                    ID = ID + '/1'
                if alignedread.is_read2:
                    ID = ID + '/2'
                if doUniqueBAM:
                    reads += 1
                else:
                    if noMulti and alignedread.opt('NH') > 1:
                        continue
                    reads += 1./alignedread.opt('NH')
#                    print 'NH, weight:', alignedread.opt('NH'), 1./alignedread.opt('NH')
        except:
            print 'problem with region:', chr, left, right, 'assigning 0 value'
            reads=0
        if doRPM:
            score = reads / normalizeBy
#            print chr, right - left, normalizeBy
        else:
            try:
                score = reads / (((right-left)/1000.)*normalizeBy)
            except:
                print 'region of size 0, skipping:', line.strip()
                continue
        if doPrintSum:
            RPMSum+=score
        outline = line.strip() +'\t' + str(score)
        if doMappabilityCorrection:
            outline = outline + '\t' + str(MappabilityRegionDict[chr][(left,right)])
            if MappabilityRegionDict[chr][(left,right)] == 0:
                outline = outline + '\t0'
            else:
                outline = outline + '\t' + str(score/MappabilityRegionDict[chr][(left,right)])
        outfile.write(outline + '\n')
         
    if doPrintSum:
        outfile.write('#Total RPM:' + str(RPMSum) + '\n')

    outfile.close()
   
run()
