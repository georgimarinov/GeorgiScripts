##################################
#                                #
# Last modified 2019/07/15       # 
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
# 0x0800 2048 supplementary alignment

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','X':'X','a':'t','t':'a','g':'c','c':'g','n':'n','x':'x','R':'R','r':'r','M':'M','m':'m','Y':'Y','y':'y','S':'S','s':'s','K':'K','k':'k','W':'W','w':'w'}
    sequence=''
    for j in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-j-1]]
    return sequence

def FLAG(FLAG):

    Numbers = [0,1,2,4,8,16,32,64,128,256,512,1024,2048]

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
        print 'usage: python %s title BAMfilename chrom.sizes outputfilename [-stranded + | -] [-end2only] [-end1only] [-nomulti] [-RPM] [-notitle] [-shift bp] [-shift_minus bp] [-shift_plus bp] [-mismatchesMD M] [-singlebasepair] [-mismatches M] [-readLength min max] [-chr chrN1(,chrN2....)] [-absValue] [-uniqueBAM] [-noNHinfo] [-biasCorrect genome.fa model.bias k-mer_size] [-readIDstring string present|absent]' % sys.argv[0]
        print '\tUse the [-mismatches] option to specify the maximum number of mismatches allowed for an alignment to be considered; use the -mimatchesMD option is mismatches are specified with the MD special tag'
        print '\tuse the [-readIDstring] option to indluce/exclude reads which contain this string; normalization will still be carried out on all reads'
        sys.exit(1)
    
    doSingleBP=False
    if '-singlebasepair' in sys.argv:
        doSingleBP=True

    doTitle=True
    if '-notitle' in sys.argv:
        doTitle=False

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
    
    doReadLength=False
    if '-readLength' in sys.argv:
        doReadLength=True
        minRL = int(sys.argv[sys.argv.index('-readLength')+1])
        maxRL = int(sys.argv[sys.argv.index('-readLength')+2])
        print 'will only consider reads between', minRL, 'and', maxRL, 'bp length'

    doRID = False
    if '-readIDstring' in sys.argv:
        doRID = True
        RID = sys.argv[sys.argv.index('-readIDstring') + 1]
        RIDpresence = sys.argv[sys.argv.index('-readIDstring') + 2]
        print 'will only output alignments with read IDs where the string ', RID, ' is ', RIDpresence

    doAbs = False
    if '-absValue' in sys.argv:
        doAbs = True
        print 'will output absolute values'

    doBCorr = False
    if '-biasCorrect' in sys.argv:
        doBCorr = True
        fasta = sys.argv[sys.argv.index('-biasCorrect') + 1]
        biasmodelfile = sys.argv[sys.argv.index('-biasCorrect') + 2]
        KS = int(sys.argv[sys.argv.index('-biasCorrect') + 3])

        KmerDict = {}
        linelist = open(biasmodelfile)
        for line in linelist:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            kmer = fields[0]
            KmerDict[kmer] = float(fields[3])

        GenomeDict={}
        sequence=''
        inputdatafile = open(fasta)
        for line in inputdatafile:
            if line[0]=='>':
                if sequence != '':
                    GenomeDict[chr] = ''.join(sequence).upper()
                chr = line.strip().split('>')[1]
                print chr
                sequence=[]
                Keep=False
                continue
            else:
                sequence.append(line.strip())
        GenomeDict[chr] = ''.join(sequence).upper()

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

    doChrSubset=False
    if '-chr' in sys.argv:
        doChrSubset=True
        WantedChrDict={}
        for chr in sys.argv[sys.argv.index('-chr')+1].split(','):
            WantedChrDict[chr]=''

    shiftplus = 0
    shiftminus = 0
    if '-shift' in sys.argv:
        shiftplus = int(sys.argv[sys.argv.index('-shift')+1])
        shiftminus = int(sys.argv[sys.argv.index('-shift')+1])
        print 'Will shfit reads by ', shiftplus, 'bp'

    if '-shift_plus' in sys.argv:
        shiftplus = int(sys.argv[sys.argv.index('-shift_plus')+1])
        print 'Will shfit + strand reads by ', shiftplus, 'bp'

    if '-shift_minus' in sys.argv:
        shiftminus = int(sys.argv[sys.argv.index('-shift_minus')+1])
        print 'Will shfit - strand reads by ', shiftminus, 'bp'

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

    doUniqueBAM = False
    if '-uniqueBAM' in sys.argv:
        print 'will treat all alignments as unique'
        doUniqueBAM = True
        TotalReads = 0
        pass

    TotalNumberRead=0

    samfile = pysam.Samfile(BAM, "rb" )

    doNoNHinfo = False
    if '-noNHinfo' in sys.argv:
        doNoNHinfo = True
        print 'will directly evaluate read mulitplicity'
        MultiplicityDict = {}
    elif doUniqueBAM:
        pass
    else:
        try:
            print 'testing for NH tags presence'
            for alignedread in samfile.fetch():
                multiplicity = alignedread.opt('NH')
                print 'file has NH tags'
                break
        except:
            print 'no NH: tags in BAM file, exiting'
            sys.exit(1)

    RN=0
    if doNoNHinfo:
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
                fields = str(alignedread).split('\t')
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
        TotalNumberRead = len(MultiplicityDict.keys())
    else:
        for (chr,start,end) in chromInfoList:
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
                    print 'counting total number of reads', str(RN/1000000) + 'M alignments processed', chr, currentPos, end
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
                fields=str(alignedread).split('\t')
                FLAGfields = FLAG(int(fields[1]))
                if 4 in FLAGfields:
                    continue
                if doEnd1Only:
                    if 128 in FLAGfields:
                        continue
                if doEnd2Only:
                    if 64 in FLAGfields:
                        continue
                if doUniqueBAM:
                    TotalNumberRead+=1
                    continue
                try:
                    multiplicity = alignedread.opt('NH')
                except:
                    print 'no NH: tags in BAM file, exiting'
                    sys.exit(1)
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
            if 4 in FLAGfields:
                continue
            if doEnd1Only:
                if 128 in FLAGfields:
                    continue
            if doEnd2Only:
                if 64 in FLAGfields:
                    continue
            ID = fields[0]
            if doRID:
                if RIDpresence == 'present':
                    if RID in ID:
                        pass
                    else:
                        continue
                if RIDpresence == 'absent':
                    if RID in ID:
                        continue
                    else:
                        pass
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
            elif doNoNHinfo:
                if alignedread.is_read1:
                    ID = ID + '/1'
                if alignedread.is_read2:
                    ID = ID + '/2'
                multiplicity = MultiplicityDict[ID]
            else:
                multiplicity = alignedread.opt('NH')
            if noMulti and multiplicity > 1:
                continue
            scaleby = 1.0/multiplicity
            FLAGfields = FLAG(int(fields[1]))
            if 16 in FLAGfields:
                s = '-'
            else:
                s = '+'
            if doStranded:
                if s != strand:
                    continue
            currentPos = alignedread.pos
            if doBCorr:
                kmer = GenomeDict[chr][currentPos - (KS/2):currentPos + (KS/2)]
                if s == '-':
                    kmer = getReverseComplement(kmer)
                if KmerDict.has_key(kmer):
                    scaleby = 1/KmerDict[kmer]
            if (doStranded and strand == '+') or s == '+':
                currentPos = currentPos + shiftplus
                if currentPos + 1 > 0 and currentPos + 1 <= end:
                    if coverageDict.has_key(currentPos+1):
                        coverageDict[currentPos+1] += scaleby
                    else:
                        coverageDict[currentPos+1] = scaleby
            if (doStranded and strand == '-') or s == '-':
                endPos = currentPos
                for (m,bp) in alignedread.cigar:
                    if m == 0:
                        endPos+=bp
                    elif m == 2:
                        endPos+=bp
                    elif m == 3:
                        endPos+=bp
                    else:
                        continue
                endPos = endPos - shiftminus
                if endPos + 1 > 0 and endPos + 1 <= end:
                    if coverageDict.has_key(endPos+1):
                        coverageDict[endPos+1]+=scaleby
                    else:
                        coverageDict[endPos+1]=scaleby

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
                    if doStranded and strand == '-' and not doAbs:
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
#                else:
#                    outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t' + str(0)
#                    outfile.write(outline+'\n')
        else:
            for i in posKeys[1:len(posKeys)]:
                if previous[0]+1 == i and previous[1]==coverageDict[i]:
                     previous=(i,coverageDict[i])
                else:
                     if written[0]==initial[0]:
                         print '####', written, initial, previous
                     if doStranded and strand == '-' and not doAbs:
                         if doRPM:
                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t-'+str(initial[1]/normFactor).split('.')[0] + '.' + str(initial[1]/normFactor).split('.')[1][0:4]
                         else:
                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t-'+str(initial[1]).split('.')[0] + '.' + str(initial[1]).split('.')[1][0:4]
                     else:
                         if doRPM:
                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t'+str(initial[1]/normFactor).split('.')[0] + '.' + str(initial[1]/normFactor).split('.')[1][0:4]
                         else:
                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t'+str(initial[1]).split('.')[0] + '.' + str(initial[1]).split('.')[1][0:4]
                     written=(initial[0],previous[0]+1)
                     outfile.write(outline+'\n')
                     initial=(i,coverageDict[i])
                     previous=(i,coverageDict[i])

    outfile.close()
            
run()
