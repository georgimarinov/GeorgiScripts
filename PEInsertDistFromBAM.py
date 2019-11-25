##################################
#                                #
# Last modified 2019/05/14       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import pysam

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
        print 'usage: python %s BAMfilename chrom.sizes outputfilename [-noNHinfo] [-nomulti] [-uniqueBAM] [-firstN number_pairs] [-chr chr1,...,chrN] [-regions file chrFiledID leftFieldID rightFieldID] [-normalize]' % sys.argv[0]
        print '\Note: the -regions option and the -chr option will be integrated if both run, i.e. only the regions within the wanted chromosomes will be used'
        sys.exit(1)

    doUniqueBAM = False
    if '-uniqueBAM' in sys.argv:
        print 'will treat all alignments as unique'
        doUniqueBAM = True
        pass

    doChr = False
    if '-chr' in sys.argv:
        doChr = True
        chromosomes = sys.argv[sys.argv.index('-chr') + 1].split(',')
        WantedDict = {}
        for chr in chromosomes:
            WantedDict[chr] = ''

    doNoNHinfo = False
    if '-noNHinfo' in sys.argv:
        doNoNHinfo = True
        MultiplicityDict = {}

    doNorm = False
    if '-normalize' in sys.argv:
        doNorm = True
        print 'will normalize counts'

    doRegions = False
    if '-regions' in sys.argv:
        doRegions = True
        regionsFile = sys.argv[sys.argv.index('-regions') + 1]
        regionsChr = int(sys.argv[sys.argv.index('-regions') + 2])
        regionsLeft = int(sys.argv[sys.argv.index('-regions') + 3])
        regionsRight = int(sys.argv[sys.argv.index('-regions') + 4])
        linelist = open(regionsFile)
        Regions = []
        for line in linelist:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chr = fields[regionsChr]
            left = int(fields[regionsLeft])
            right = int(fields[regionsRight])
            if doChr:
                if WantedDict.has_key(chr):
                    Regions.append((chr,left,right))
            else:
                Regions.append((chr,left,right))
    
    BAM = sys.argv[1]
    chrominfo=sys.argv[2]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        if doChr:
            if WantedDict.has_key(chr):
                chromInfoList.append((chr,start,end))
        else:
            chromInfoList.append((chr,start,end))
    outfilename = sys.argv[3]

    noMulti = False
    if '-nomulti' in sys.argv:
        noMulti = True

    doFirstN = False
    if '-firstN' in sys.argv:
        doFirstN = True
        FN = int(sys.argv[sys.argv.index('-firstN') + 1])

    InsertLengthDistribution = {}
    InsertLengthDistribution['singleton'] = 0

    samfile = pysam.Samfile(BAM, "rb" )

    if doNoNHinfo:
        MultiplicityDict = {}
        i=0
        for (chr,start,end) in chromInfoList:
            try:
                jj=0
                for alignedread in samfile.fetch(chr, start, end):
                    jj+=1
                    if jj==1:
                        break
            except:
                print 'problem with region:', chr, start, end, 'skipping'
                continue
            for alignedread in samfile.fetch(chr, start, end):
                i+=1
                if i % 5000000 == 0:
                    print str(i/1000000) + 'M alignments processed in multiplicity assessment', chr,start,alignedread.pos,end
                fields=str(alignedread).split('\t')
                ID=fields[0]
                if alignedread.is_read1:
                    ID = ID + '/1'
                if alignedread.is_read2:
                    ID = ID + '/2'
                if MultiplicityDict.has_key(ID):
                    pass
                else:
                    MultiplicityDict[ID] = 0
                MultiplicityDict[ID] += 1

    RN=0
    PN=0
    if doRegions:
        for (chr,start,end) in Regions:
            try:
                for alignedread in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                print 'region', chr,start,end, 'not found in bam file, skipping'
                continue
            currentPos=0
            if doFirstN and PN >= FN:
                break
            for alignedread in samfile.fetch(chr, start, end):
                RN+=1
                if RN % 5000000 == 0:
                    print str(RN/1000000) + 'M alignments processed', chr, currentPos, end
                if doFirstN and PN >= FN:
                    break
                if doUniqueBAM:
                    multiplicity = 1
                elif doNoNHinfo:
                    fields=str(alignedread).split('\t')
                    ID=fields[0]
                    if alignedread.is_read1:
                        ID = ID + '/1'
                    if alignedread.is_read2:
                        ID = ID + '/2'
                    multiplicity = MultiplicityDict[ID]
                else:
                    try:
                        multiplicity = alignedread.opt('NH')
                    except:
                        print 'no NH: tags in BAM file, exiting'
                        sys.exit(1)
                if noMulti and multiplicity > 1:
                    continue
                fields=str(alignedread).split('\t')
                FLAGfields = FLAG(int(fields[1]))
                pos = alignedread.pos
                if 8 in FLAGfields:
                    InsertLengthDistribution['singleton'] += 1
                    PN+=1
                    continue
                matepos = alignedread.pnext
                if matepos > pos:
                    continue
                IL = pos - matepos + len(alignedread.query)
                if InsertLengthDistribution.has_key(IL):
                    pass
                else:
                    InsertLengthDistribution[IL] = 0
                InsertLengthDistribution[IL] += 1
                PN+=1
    else:
        for (chr,start,end) in chromInfoList:
            try:
                for alignedread in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                print 'region', chr,start,end, 'not found in bam file, skipping'
                continue
            currentPos=0
            if doFirstN and PN >= FN:
                break
            for alignedread in samfile.fetch(chr, start, end):
                RN+=1
                if RN % 5000000 == 0:
                    print str(RN/1000000) + 'M alignments processed', chr, currentPos, end
                if doFirstN and PN >= FN:
                    break
                if doUniqueBAM:
                    multiplicity = 1
                elif doNoNHinfo:
                    fields=str(alignedread).split('\t')
                    ID=fields[0]
                    if alignedread.is_read1:
                        ID = ID + '/1'
                    if alignedread.is_read2:
                        ID = ID + '/2'
                    multiplicity = MultiplicityDict[ID]
                else:
                    try:
                        multiplicity = alignedread.opt('NH')
                    except:
                        print 'no NH: tags in BAM file, exiting'
                        sys.exit(1)
                if noMulti and multiplicity > 1:
                    continue
                fields=str(alignedread).split('\t')
                FLAGfields = FLAG(int(fields[1]))
                pos = alignedread.pos
                if 8 in FLAGfields:
                    InsertLengthDistribution['singleton'] += 1
                    PN+=1
                    continue
                matepos = alignedread.pnext
                if matepos > pos:
                    continue
                IL = pos - matepos + len(alignedread.query)
                if InsertLengthDistribution.has_key(IL):
                    pass
                else:
                    InsertLengthDistribution[IL] = 0
                InsertLengthDistribution[IL] += 1
                PN+=1

    outfile = open(outfilename, 'w')

    outline = '#Length\tNumberPairs'
    outfile.write(outline + '\n')

    keys = InsertLengthDistribution.keys()
    keys.sort()

    if doNorm:
        Total = 0.0
        for IL in keys:
            Total += InsertLengthDistribution[IL]

    for IL in keys:
        if doNorm:
            outline = str(IL) + '\t' + str(InsertLengthDistribution[IL]/Total)
        else:
            outline = str(IL) + '\t' + str(InsertLengthDistribution[IL])
        outfile.write(outline + '\n')

    outfile.close()
            
run()
