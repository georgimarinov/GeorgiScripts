##################################
#                                #
# Last modified 2018/04/11       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import pysam
import os

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

    if len(sys.argv) < 7:
        print 'usage: python %s BAMfilename positions_file chrFieldID positionFieldID radius maxFragmentLlength outputfilename [-stranded fieldID] [-nomulti] [-narrowPeak] [-chr chr1,...,chrN] ' % sys.argv[0]
        print '\tNote: the regions file can be zipped'
        print '\tNote: positionFieldID can be also middle'
        sys.exit(1)

    BAM = sys.argv[1]
    regions = sys.argv[2]
    chrFieldID = int(sys.argv[3])
    if sys.argv[4] == 'middle':
        posFieldID = 'middle'
    else:
        posFieldID = int(sys.argv[4])
    radius = int(sys.argv[5])
    mFL = int(sys.argv[6])
    outfilename = sys.argv[7]

    doStranded = False
    if '-stranded' in sys.argv:
        doStranded = True
        strandFieldID = int(sys.argv[sys.argv.index('-stranded') + 1])
        print 'will treat regions as stranded, strand field:', strandFieldID

    noMulti = False
    if '-nomulti' in sys.argv:
        noMulti = True

    doNP = False
    if '-narrowPeak' in sys.argv:
        doNP = True

    doChr = False
    if '-chr' in sys.argv:
        doChr = True
        chromosomes = sys.argv[sys.argv.index('-chr') + 1].split(',')
        WantedDict = {}
        for chr in chromosomes:
            WantedDict[chr] = ''

    InsertLengthMatrix = {}
    for i in range(-radius,radius+1):
        InsertLengthMatrix[i] = {}
        for j in range(mFL):
            InsertLengthMatrix[i][j] = 0

    samfile = pysam.Samfile(BAM, "rb" )

    if regions.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + regions
    elif regions.endswith('.gz'):
        cmd = 'gunzip -c ' + regions
    elif regions.endswith('.zip'):
        cmd = 'unzip -p ' + regions
    else:
        cmd = 'cat ' + regions
    p = os.popen(cmd, "r")
    line = 'line'
    RP = 0
    while line != '':
        line = p.readline().strip()
        fields = line.split('\t')
        RP += 1
        if RP % 1000 == 0:
            print RP, 'regions processed'
        if line == '':
            break
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chr = fields[chrFieldID]
        if doChr:
            if WantedDict.has_key(chr):
                pass
            else:
                continue
        if doNP:
            peak = int(fields[1]) + int(fields[9])
        elif posFieldID == 'middle':
            peak = int((int(fields[chrFieldID + 1]) + int(fields[chrFieldID + 2]))/2.)
        else:
            peak = int(fields[posFieldID])
        strand = '+'
        if doStranded:
            strand = fields[strandFieldID]
        TEST = 0
        try:
            for alignedread in samfile.fetch(chr, max(0,peak-radius-mFL), peak+radius+mFL):
                TEST += 1
                if TEST >= 1:
                    break
        except:
            print 'problem with region:'
            print fields
            print chr, max(0,peak-radius-mFL), peak+radius+mFL
            continue
        for alignedread in samfile.fetch(chr, max(0,peak-radius-mFL), peak+radius+mFL):
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
                continue
            matepos = alignedread.pnext
            if matepos > pos:
                continue
            IL = pos - matepos + len(alignedread.query)
            FP = int(matepos + IL/2.)
            if FP < peak-radius or FP >= peak+radius or IL >= mFL:
                continue
            if strand == '+':
                relativepos = FP - peak
            if strand == '-':
                relativepos = peak - FP
            InsertLengthMatrix[relativepos][IL] += 1

    TotalFrags = 0.0
    for i in InsertLengthMatrix.keys():
        for j in InsertLengthMatrix[i].keys():
            TotalFrags += InsertLengthMatrix[i][j]

    NormFactor = TotalFrags/1000000

    outfile = open(outfilename, 'w')

    outline = '#'
    for i in range(-radius,radius):
        outline = outline + '\t' + str(i)
    outfile.write(outline + '\n')

    for j in range(mFL):
        IL = mFL - j - 1
        outline = str(IL)
        for i in range(-radius,radius):
            outline = outline + '\t' + str(InsertLengthMatrix[i][IL]/NormFactor)
        outfile.write(outline + '\n')

    outfile.close()
            
run()
