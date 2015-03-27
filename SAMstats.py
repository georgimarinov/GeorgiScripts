##################################
#                                #
# Last modified 03/21/2015       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import gc
import pysam
import string
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

    if len(sys.argv) < 2:
        print 'usage: python %s SAMfilename outputfilename [-bam chrom.sizes samtools] [-paired]' % sys.argv[0]
        print '       BAM file has to be indexed'
        print '       Complexity will be calculated only for BAM files'
        sys.exit(1)

    SAM = sys.argv[1]
    outputfilename = sys.argv[2]

    doBAM=False
    if '-bam' in sys.argv:
        doBAM=True
        chrominfo=sys.argv[sys.argv.index('-bam')+1]
        samtools=sys.argv[sys.argv.index('-bam')+2]
        chromInfoList=[]
        linelist=open(chrominfo)
        for line in linelist:
            fields=line.strip().split('\t')
            chr=fields[0]
            start=0
            end=int(fields[1])
            chromInfoList.append((chr,start,end))
        samfile = pysam.Samfile(SAM, "rb" )
        try:
            print 'testing for NH tags presence'
            for alignedread in samfile.fetch():
                multiplicity = alignedread.opt('NH')
                print 'file has NH tags'
                break
        except:
            print 'no NH: tags in BAM file, will replace with a new BAM file with NH tags'
            BAMpreporcessingScript = sys.argv[0].rpartition('/')[0] + '/bamPreprocessing.py'
            cmd = 'python ' + BAMpreporcessingScript + ' ' + SAM + ' ' + SAM + '.NH'
            os.system(cmd)
            cmd = 'rm ' + SAM
            os.system(cmd)
            cmd = 'mv ' + SAM + '.NH' + ' ' + SAM
            os.system(cmd)
            cmd = samtools + ' index ' + SAM
            os.system(cmd)

    ReadDict={}

    Unique=0
    UniqueSplices=0
    Multi=0
    MultiSplices=0

    SeenDict={}
    SeenTwiceDict={}
    doPaired=False
    if '-paired' in sys.argv:
        doPaired=True
        print 'will treat reads as paired'
        SeenDictPaired={}
        SeenDictPairedSpliced={}

    print 'examining read multiplicty'

    ReadLengthDict={}

    if doBAM:
        i=0
        samfile = pysam.Samfile(SAM, "rb" )
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
                    print str(i/1000000) + 'M alignments processed', chr,start,alignedread.pos,end
                fields=str(alignedread).split('\t')
                ID=fields[0]
                length=len(alignedread.seq)
                if ReadLengthDict.has_key(length):
                    ReadLengthDict[length]+=1
                else:
                    ReadLengthDict[length]=1
                if alignedread.is_read1:
                    ID = ID + '/1'
                if alignedread.is_read2:
                    ID = ID + '/2'
                multiplicity = alignedread.opt('NH')
                if multiplicity > 1:
                    if SeenTwiceDict.has_key(ID):
                        continue
                    SeenTwiceDict[ID]=''
                    if len(alignedread.cigar) > 1:
                        Splice=False
                        for (m,bp) in alignedread.cigar:
                            if m == 3:
                               MultiSplices+=1
                               Splice=True
                               break
                        if not Splice:
                            Multi+=1
                    else:
                        Multi+=1
                if multiplicity == 1:
                    if alignedread.cigar == None:
                        Unique+=1
                        continue
                    if len(alignedread.cigar) > 1:
                        Splice=False
                        for (m,bp) in alignedread.cigar:
                            if m == 3:
                                UniqueSplices+=1
                                Splice=True
                                break
                        if not Splice:
                            Unique+=1
                    else:
                        Unique+=1
    else:
        i=0
        lineslist = open(SAM)
        for line in lineslist:
            if line.startswith('@'):
                continue
            i+=1
            if i % 5000000 == 0:
                print str(i/1000000) + 'M alignments processed'
            fields = line.strip().split('\t')
            length=len(fields[6])
            if ReadLengthDict.has_key(length):
                ReadLengthDict[length]+=1
            else:
                ReadLengthDict[length]=1
            if fields[2] == '*':
                continue
            ID=fields[0]
            FLAGfields = FLAG(int(fields[1]))
            if 64 in FLAGfields:
                ID = ID + '/1'
            elif 128 in FLAGfields:
                ID = ID + '/2'
            else:
                 print 'paired information incorrectly specified, exiting'
                 sys.exit(1)
            if SeenDict.has_key(ID):
                if SeenTwiceDict.has_key(ID):
                    continue
                SeenTwiceDict[ID]=''
                if 'N' in fields[5]:
                    UniqueSplices-=1
                    MultiSplices+=1
                else:
                    Unique-=1
                    Multi+=1
            else:
                SeenDict[ID]=''
                if 'N' in fields[5]:
                    UniqueSplices+=1
                else:
                    Unique+=1

    TotalUniqueReads = 0.0
    DistinctUniqueReads = 0

    Complexity='N\A'
    if doBAM:
        TotalUniqueReads = 0
        DistinctUniqueReads = 0
        i=0
        for (chr,start,end) in chromInfoList:
            try:
                has_reads = False
                for alignedread in samfile.fetch(chr, start, end):
                    has_reads = True
                    if has_reads:
                        break
            except:
                continue
            currentPos=0
            CurrentPosDictPlus=[]
            CurrentPosDictMinus=[]
            for alignedread in samfile.fetch(chr, start, end):
                i+=1
                if i % 5000000 == 0:
                    print str(i/1000000) + 'M alignments processed in complexity calculation', chr,start,alignedread.pos,end
                    gc.collect()
                fields=str(alignedread).split('\t')
                if 16 in FLAG(alignedread.flag):
                    strand = '-'
                else:
                    strand = '+'
                ID=fields[0]
                if alignedread.is_read1:
                    ID = ID + '/1'
                if alignedread.is_read2:
                    ID = ID + '/2'
                multiplicity = alignedread.opt('NH')
                if multiplicity > 1:
                    continue
                start = alignedread.pos
                if start != currentPos:
                    TotalUniqueReads += len(CurrentPosDictPlus)
                    DistinctUniqueReads += len(Set(CurrentPosDictPlus))
                    TotalUniqueReads += len(CurrentPosDictMinus)
                    DistinctUniqueReads += len(list(Set(CurrentPosDictMinus)))
                    CurrentPosDictPlus=[]
                    CurrentPosDictMinus=[]
                    currentPos=start
                if strand == '+':
                    CurrentPosDictPlus.append(str(alignedread.cigar))
                    if doPaired:
                         CurrentPosDictPlus.append((str(alignedread.cigar),alignedread.mpos))
                if strand == '-':
                    CurrentPosDictMinus.append(str(alignedread.cigar))
                    if doPaired:
                         CurrentPosDictMinus.append((str(alignedread.cigar),alignedread.mpos))
        Complexity = DistinctUniqueReads/float(TotalUniqueReads) 

    SeenDict=''
    SeenTwiceDict=''

    outfile=open(outputfilename, 'w')

    outline='Unique:\t'+str(Unique)
    print outline
    outfile.write(outline+'\n')

    outline='Unique Splices:\t'+str(UniqueSplices)
    print outline
    outfile.write(outline+'\n')

    outline='Multi:\t'+str(Multi)
    print outline
    outfile.write(outline+'\n')

    outline='Multi Splices:\t'+str(MultiSplices)
    print outline
    outfile.write(outline+'\n')

    outline='Complexity:\t'+str(Complexity)[0:4]
    print outline
    outfile.write(outline+'\n')

    outline='Read Length, Minimum:\t'+str(min(ReadLengthDict.keys()))
    print outline
    outfile.write(outline+'\n')

    outline='Read Length, Maximum:\t'+str(max(ReadLengthDict.keys()))
    print outline
    outfile.write(outline+'\n')

    TotalReads=0.0
    TotalLength=0.0
    for length in ReadLengthDict.keys():
        TotalReads += ReadLengthDict[length]
        TotalLength += length*ReadLengthDict[length]

    outline='Read Length, Average:\t'+str(TotalLength/TotalReads)
    print outline
    outfile.write(outline+'\n')
             
    outfile.close()

run()
