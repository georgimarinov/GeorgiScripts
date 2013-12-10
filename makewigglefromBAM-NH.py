##################################
#                                #
# Last modified 09/06/2013       # 
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

    if len(sys.argv) < 4:
        print 'usage: python %s title BAMfilename chrom.sizes outputfilename [-stranded + | -] [-fragments first-read-strand | second-read-strand] [-shift bp] [-nomulti] [-RPM] [-notitle] [-singlebasepair] [-mismatchesMD M] [-mismatches M] [-end2only] [-end1only] [-readLength min max] [-chr chrN1(,chrN2....)] [-uniqueBAM] [-noNH samtools]' % sys.argv[0]
        print '\tUse the -mismatches option to specify the maximum number of mismatches allowed for an alignment to be considered; use the -mimatchesMD option is mismatches are specified with the MD special tag'
        print '\tuse the -noNH option and supply a path to samtools in order to have the file converted to one that has NH tags'
        sys.exit(1)

    title = sys.argv[1]
    BAM = sys.argv[2]
    chrominfo=sys.argv[3]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))
    outfilename = sys.argv[4]

    doEnd1Only = False
    doEnd2Only = False
    if '-end1only' in sys.argv and '-end2only' in sys.argv:
        print 'both -end1only and -end2only option specified, a logical impossiblity, exiting'
        sys,exit(1)

    if '-end1only' in sys.argv:
        doEnd1Only = True
        print 'will only consider the first end of read pairs'

    if '-end2only' in sys.argv:
        doEnd2Only = True 
        print 'will only consider the second end of read pairs'
    
    doSingleBP=False
    if '-singlebasepair' in sys.argv:
        doSingleBP=True

    doTitle=True
    if '-notitle' in sys.argv:
        doTitle=False

    doUniqueBAM = False
    if '-uniqueBAM' in sys.argv:
        print 'will treat all alignments as unique'
        doUniqueBAM = True
        TotalReads = 0
        pass

    samfile = pysam.Samfile(BAM, "rb" )
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
            cmd = 'python ' + BAMpreporcessingScript + ' ' + BAM + ' ' + BAM + '.NH'
            os.system(cmd)
            cmd = 'rm ' + BAM
            os.system(cmd)
            cmd = 'mv ' + BAM + '.NH' + ' ' + BAM
            os.system(cmd)
            cmd = samtools + ' index ' + BAM
            os.system(cmd)
        elif doUniqueBAM:
            pass
        else:
            print 'no NH: tags in BAM file, exiting'
            sys.exit(1)

    doReadLength=False
    if '-readLength' in sys.argv:
        doReadLength=True
        minRL = int(sys.argv[sys.argv.index('-readLength')+1])
        maxRL = int(sys.argv[sys.argv.index('-readLength')+2])
        print 'will only consider reads between', minRL, 'and', maxRL, 'bp length'

    doMaxMMMD=False
    if '-mismatchesMD' in sys.argv:
        doMaxMMMD=True
        maxMM = int(sys.argv[sys.argv.index('-mismatchesMD')+1])
        print 'Will only consider alignments with', maxMM, 'or less mismatches'

    doMaxMM=False
    if '-mismatches' in sys.argv:
        doMaxMM=True
        maxMM = int(sys.argv[sys.argv.index('-mismatches')+1])
        print 'Will only consider alignments with', maxMM, 'or less mismatches'

    shift = 0
    if '-shift' in sys.argv:
        shift = int(sys.argv[sys.argv.index('-shift')+1])
        print 'Will shfit reads by ', shift, 'bp'

    doChrSubset=False
    if '-chr' in sys.argv:
        doChrSubset=True
        WantedChrDict={}
        for chr in sys.argv[sys.argv.index('-chr')+1].split(','):
            WantedChrDict[chr]=''

    noMulti=False
    if '-nomulti' in sys.argv:
        print 'will only consider unique alignments'
        noMulti=True

    doRPM=False
    if '-RPM' in sys.argv:
        doRPM=True

    doStranded=False
    if '-stranded' in sys.argv:
        doStranded=True
        strand=sys.argv[sys.argv.index('-stranded')+1]
        print 'will only consider', strand, 'strand reads'

    doFragments=False
    if '-fragments' in sys.argv:
        doFragments=True
        FragmentStrandAssignment = sys.argv[sys.argv.index('-fragments')+1]
        print 'will assign strand for both reads based on', FragmentStrandAssignment

    TotalNumberRead=0

    if doUniqueBAM and not doReadLength and not doMaxMMMD and not doMaxMM:
        TotalNumberRead = 0
        for chrStats in pysam.idxstats(BAM):
            fields = chrStats.strip().split('\t')
            chr = fields[0]
            reads = int(fields[2])
            if chr != '*':
                TotalNumberRead += reads
    else:
        samfile = pysam.Samfile(BAM, "rb" )
        RN=0
        for (chr,start,end) in chromInfoList:
            try:
                for alignedread in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                print 'region', chr,start,end, 'not found in bam file, skipping'
                continue
            currentPos=0
            for alignedread in samfile.fetch(chr, start, end):
                RN+=1
                if RN % 5000000 == 0:
                    print 'counting total number of reads', str(RN/1000000) + 'M alignments processed', chr, currentPos, end
                fields=str(alignedread).split('\t')
                FLAGfields = FLAG(int(fields[1]))
                if doEnd1Only:
                    if 128 in FLAGfields:
                        continue
                if doEnd2Only:
                    if 64 in FLAGfields:
                        continue
                if doReadLength:
                    if len(alignedread.seq) > maxRL or len(alignedread.seq) < minRL:
                        continue
                if doMaxMM:
                    mismatches = 0
                    for (m,bp) in alignedread.cigar:
                        if m == 8:
                            mismatches+=1
                    if mismatches > maxMM:
                        continue
                if doMaxMMMD:
                    MM = alignedread.opt('MD')
                    mismatches = 0
                    if MM.isdigit():
                        pass
                    else:
                        for s in range(len(MM)):
                            if MM[s].isalpha():
                                mismatches+=1
                    if mismatches > maxMM:
                        continue
                if doUniqueBAM:
                    TotalNumberRead+=1
                    continue
                multiplicity = alignedread.opt('NH')
                if noMulti and multiplicity > 1:
                    continue
                TotalNumberRead += 1.0/multiplicity

        TotalNumberRead = round(TotalNumberRead)

    print 'found', TotalNumberRead, 'reads'
    normFactor = TotalNumberRead/1000000.
    print 'RPM normalization Factor =', normFactor

    outfile = open(outfilename, 'w')
    if doTitle:
        outline='track type=bedGraph name="' + title + '"'
        outfile.write(outline+'\n')

    RN=0

    for (chr,start,end) in chromInfoList:
        coverageDict={}
        if doChrSubset:
            if WantedChrDict.has_key(chr):
                pass
            else:
                continue
        try:
            for alignedread in samfile.fetch(chr, 0, 100):
                a='b'
        except:
            print 'region', chr,start,end, 'not found in bam file, skipping'
            continue
        currentPos=0
        for alignedread in samfile.fetch(chr, start, end):
            RN+=1
            if RN % 5000000 == 0:
                print str(RN/1000000) + 'M alignments processed', chr, currentPos, end
            if doReadLength:
                if len(alignedread.seq) > maxRL or len(alignedread.seq) < minRL:
                    continue
            fields=str(alignedread).split('\t')
            FLAGfields = FLAG(int(fields[1]))
            if doEnd1Only:
                if 128 in FLAGfields:
                    continue
            if doEnd2Only:
                if 64 in FLAGfields:
                    continue
            ID = fields[0]
            if doMaxMM:
                mismatches = 0
                for (m,bp) in alignedread.cigar:
                    if m == 8:
                        mismatches+=1
                if mismatches > maxMM:
                    continue
            if doMaxMMMD:
                MM = alignedread.opt('MD')
                mismatches = 0
                if MM.isdigit():
                    pass
                else:
                    for s in range(len(MM)):
                        if MM[s].isalpha():
                            mismatches+=1
                if mismatches > maxMM:
                    continue
            if doUniqueBAM:
                multiplicity = 1
            else:
                multiplicity = alignedread.opt('NH')
            if noMulti and multiplicity > 1:
                continue
            scaleby=1.0/multiplicity
            FLAGfields = FLAG(int(fields[1]))
            if alignedread.is_reverse:
                s = '-'
                shift = (0 - shift)
            else:
                s = '+'
            if doFragments:
                if FragmentStrandAssignment == 'first-read-strand':
                    if alignedread.is_read2:
                        if s == '+':
                            s = '-'
                        elif s == '-':
                            s = '+'
                if FragmentStrandAssignment == 'second-read-strand':
                    if alignedread.is_read1:
                        if s == '+':
                            s = '-'
                        elif s == '-':
                            s = '+'
            if doStranded:
                if s != strand:
                    continue
            currentPos=alignedread.pos
            for (m,bp) in alignedread.cigar:
                if m == 0:
                    for j in range(currentPos,currentPos+bp):
                        if coverageDict.has_key(j+1 + shift):
                            coverageDict[j+1 + shift]+=scaleby
                        else:
                            coverageDict[j+1 + shift]=scaleby
                elif m == 2:
                    pass
                elif m == 3:
                    pass
                else:
                    continue
                currentPos=currentPos+bp
        posKeys=coverageDict.keys()
        posKeys.sort()
        if len(posKeys) == 0:
            continue
        initial=(posKeys[0],coverageDict[posKeys[0]])
        previous=(posKeys[0],coverageDict[posKeys[0]])
        written=['']
        if doSingleBP:
            for i in range(1,max(posKeys)+1):
                if coverageDict.has_key(i):
                    if doStranded and strand == '-':
                        if doRPM:
                            outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t-' + str(coverageDict[i]/normFactor)
                        else:
                            outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t-' + str(coverageDict[i])
                    else:
                        if doRPM:
                            outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t' + str(coverageDict[i]/normFactor)
                        else:
                            outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t' + str(coverageDict[i])
                    outfile.write(outline+'\n')
                else:
                    outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t' + str(0)
                    outfile.write(outline+'\n')
        else:
            for i in posKeys[1:len(posKeys)]:
                if previous[0]+1 == i and previous[1]==coverageDict[i]:
                     previous=(i,coverageDict[i])
                else:
                     if written[0]==initial[0]:
                         print '####', written, initial, previous
                     try:
                         if doStranded and strand == '-':
                             if doRPM:
                                 outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t-'+str('{0:.10f}'.format(initial[1]/normFactor)).split('.')[0] + '.' + str('{0:.10f}'.format(initial[1]/normFactor)).split('.')[1][0:4]
                             else:
                                 outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t-'+str('{0:.10f}'.format(initial[1])).split('.')[0] + '.' + str('{0:.10f}'.format(initial[1])).split('.')[1][0:4]
                         else:
                             if doRPM:
                                 outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t'+str('{0:.10f}'.format(initial[1]/normFactor)).split('.')[0] + '.' + str('{0:.10f}'.format(initial[1]/normFactor)).split('.')[1][0:4]
                             else:
                                 outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t'+str('{0:.10f}'.format(initial[1])).split('.')[0] + '.' + str('{0:.10f}'.format(initial[1])).split('.')[1][0:4]
                     except:
                         print initial[0]-1
                         print previous[0]+1-1
                         print initial[1]
                         print str(initial[1]).split('.')[0]
                         print str(initial[1]).split('.')[1]
                         print str(initial[1]).split('.')[1][0:4]
                         sys.exit(1)
#                     if doStranded and strand == '-':
#                         if doRPM:
#                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t-'+str(initial[1]/normFactor).split('.')[0] + '.' + str(initial[1]/normFactor).split('.')[1][0:4]
#                         else:
#                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t-'+str(initial[1]).split('.')[0] + '.' + str(initial[1]).split('.')[1][0:4]
#                     else:
#                         if doRPM:
#                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t'+str(initial[1]/normFactor).split('.')[0] + '.' + str(initial[1]/normFactor).split('.')[1][0:4]
#                         else:
#                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t'+str(initial[1]).split('.')[0] + '.' + str(initial[1]).split('.')[1][0:4]
                     written=(initial[0],previous[0]+1)
                     outfile.write(outline+'\n')
                     initial=(i,coverageDict[i])
                     previous=(i,coverageDict[i])

    outfile.close()
            
run()
