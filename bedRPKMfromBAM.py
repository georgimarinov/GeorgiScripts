##################################
#                                #
# Last modified 2019/03/08       # 
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
        print 'usage: python %s bedfilename chrField BAMfilename chrom.sizes outputfilename [-endonly] [-end2only] [-end1only] [-nomulti] [-RPM] [-totalReadNumber number] [-stranded +|-] [-readLength min max] [-printSum] [-uniqueBAM] [-mappabilityNormalize mappability.wig readLength] [-noNH samtools] [-singleFieldRegion] [-excludeReadsMappingToOtherChromosomes]' % sys.argv[0]
        print 'Note: the script will divide multireads by their multiplicity'
        print '\tthe [-printSum] option only working together with the RPM option'
        print '\tuse the [-uniqueBAM] option if the BAM file contains only unique alignments; this will save a lot of memory'
        print '\tuse the [-mappabilityNormalize] option to get mappability normalized RPKMs (it will not do anything to the RPMs; not that a mappability track that goes from 0 to the read length is assumed'
        print '\tuse the [-noNH] option and supply a path to samtools in order to have the file converted to one that has NH tags'
        print '\tthe stranded option will normalized against all reads, not just reads on the indicated strand'
        print '\tthe [-endonly] option will only count reads whose leftmost coordinates are in a given region'
        print '\tuse the [-excludeReadsMappingToOtherChromosomes] option if you want to exclude multimappers that map to chromosomes other than what is included in the chrom.sizes file; note that it is incompatible with the [-noNHinfo] option'
        sys.exit(1)

    
    bed = sys.argv[1]
    fieldID = int(sys.argv[2])
    SAM = sys.argv[3]
    chromSize = sys.argv[4]
    outfilename = sys.argv[5]

    chromInfoList = []
    chromInfoDict = {}
    linelist = open(chromSize)
    for line in linelist:
        fields = line.strip().split('\t')
        chr = fields[0]
        start = 0
        end = int(fields[1])
        chromInfoList.append((chr,start,end))
        chromInfoDict[chr] = 1

    doERMTOC = False
    if '-excludeReadsMappingToOtherChromosomes' in sys.argv:
        print 'will exclude multimapping reads mapping to chromosomes other than those included in', chromSize
        doERMTOC = True
        ERMTOCDict = {}

    doEO = False
    if '-endonly' in sys.argv:
        doEO = True
        print '-endonly option enabled'

    doEnd1Only = False
    doEnd2Only = False
    if '-end1only' in sys.argv and '-end2only' in sys.argv:
        print 'both -end1only and -end2only option specified, a logical impossiblity, exiting'
        sys.exit(1)

    if '-end1only' in sys.argv:
        doEnd1Only = True
        print 'will only consider the first end of read pairs'

    if '-end2only' in sys.argv:
        doEnd2Only = True 
        print 'will only consider the second end of read pairs'
    
    noMulti=False
    if '-nomulti' in sys.argv:
        noMulti=True
        print 'will discard multi-read alignments'

    doDirectTRN = False
    if '-totalReadNumber' in sys.argv:
        doDirectTRN = True
        DTRN = int(sys.argv[sys.argv.index('-totalReadNumber')+1])

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

    doSFR = False
    if '-singleFieldRegion' in sys.argv:
        doSFR = True

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

    if doUniqueBAM and not doReadLength and not doEnd2Only and not doEnd1Only:
        TotalReads = 0
        try:
            for chrStats in pysam.idxstats(SAM):
                fields = chrStats.strip().split('\t')
                chr = fields[0]
                reads = int(fields[2])
                if chr != '*':
                    TotalReads += reads
        except:
            for chrStats in pysam.idxstats(SAM).strip().split('\n'):
#                print chrStats
                fields = chrStats.strip().split('\t')
                print fields
                chr = fields[0]
                reads = int(fields[2])
                if chr != '*':
                    TotalReads += reads
        UniqueReads = TotalReads
    elif doDirectTRN:
        TotalReads = DTRN
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
                    ID=fields[0]
                    if alignedread.is_read1:
                        if doEnd2Only:
                             continue
                        ID = ID + '/1'
                    if alignedread.is_read2:
                        if doEnd1Only:
                             continue
                        ID = ID + '/2'
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

    normalizeBy = TotalReads/1000000.

    print 'TotalReads', TotalReads
    print 'RPM normalization Factor =', normalizeBy

    if doERMTOC:
        i = 0
        samfile = pysam.Samfile(SAM, "rb" )
        for read in samfile.fetch(until_eof=True):
            i+=1
            if i % 5000000 == 0:
                print 'examining read cross-chromosome alignments, first pass', str(i/1000000) + 'M alignments processed processed'
            fields = str(read).split('\t')
            ID = read.qname
            if read.is_unmapped:
                continue
            if read.opt('NH') == 1:
                continue
            if ERMTOCDict.has_key(ID):
                pass
            else:
                ERMTOCDict[ID] = {}
            chr = samfile.getrname(read.tid)
            ERMTOCDict[ID][chr] = 1.0/read.opt('NH')
        i = 0
#        print 'found', len(ERMTOCDict.keys()), 'multimappers'
        Excluded = 0
        for ID in ERMTOCDict.keys():
            i+=1
            if i % 5000000 == 0:
                 print 'examining read cross-chromosome alignments, second pass', str(i/1000000) + 'M alignments processed processed'
            ToBeExcluded = False
            AlignsToWantedChromosomes = False
            for chr in ERMTOCDict[ID].keys():
                if chr not in chromInfoDict.keys():
                    ToBeExcluded = True
                if chr in chromInfoDict.keys():
                    AlignsToWantedChromosomes = True
            if ToBeExcluded:
                if AlignsToWantedChromosomes:
                    TotalReads = TotalReads - ERMTOCDict[ID][chr]
                    Excluded += 1
                del ERMTOCDict[ID]
        print 'excuded', Excluded, 'multimappers'
        print 'retained', len(ERMTOCDict.keys()), 'multimappers'

        normalizeBy = TotalReads/1000000.

        print 'TotalReads after excluding multimappers mapping to chromosomes outside provided list', TotalReads
        print 'RPM normalization Factor =', normalizeBy

    outfile = open(outfilename, 'w')

    i=0
    if bed.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + bed
    elif bed.endswith('.gz'):
        cmd = 'gunzip -c ' + bed
    else:
        cmd = 'cat ' + bed
    p = os.popen(cmd, "r")
    line = 'line'
    while line != '':
        line = p.readline()
        if line == '':
            break
#    lineslist = open(bed)
#    for line in lineslist:
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
#        reads2 =  0
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
                    if doEnd2Only:
                        continue
                    ID = ID + '/1'
                if alignedread.is_read2:
                    if doEnd1Only:
                        continue
                    ID = ID + '/2'
                if doUniqueBAM:
                    if doEO:
                        if alignedread.pos >= left and alignedread.pos < right:
                            reads += 1
                        else:
                            reads += 0
                    else:
                        reads += 1
                else:
                    if noMulti and alignedread.opt('NH') > 1:
                        continue
                    if doERMTOC and alignedread.opt('NH') > 1:
                        if ERMTOCDict.has_key(ID):
                            pass
                        else:
                            continue
                    if doEO:
                        if alignedread.pos >= left and alignedread.pos < right:
                            reads += 1./alignedread.opt('NH')
                        else:
                            reads += 0
                    else:
                        reads += 1./alignedread.opt('NH')
        except:
            print 'problem with region:', chr, left, right, 'assigning 0 value'
            reads=0
        if doRPM:
            score = reads/normalizeBy
        else:
            try:
                score = reads / (((right-left)/1000.)*normalizeBy)
            except:
                print 'region of size 0, skipping:', line.strip()
                continue
        if doPrintSum:
            RPMSum += score
        if doSFR:
            outline = chr + ':' + str(left) + '-' + str(right) +'\t' + str(score)
        else:
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
