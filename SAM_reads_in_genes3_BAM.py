##################################
#                                #
# Last modified 05/19/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
import pysam

def run():

    if len(sys.argv) < 4:
        print 'usage: python %s gtf BAM chrom.sizes outputfilename [-nomulti] [-nounique] [-noNH] ' % sys.argv[0]
        print '  by default, the script is not designed to deal with multi reads unless the NH: field is present; use the [-noNH] option if it is not; multi reads will be weighed by their multiplicity' 
        sys.exit(1)
    
    gtf = sys.argv[1]
    SAM = sys.argv[2]
    chrominfo = sys.argv[3]
    outfilename = sys.argv[4]

    noMulti=False
    if '-nomulti' in sys.argv:
        noMulti=True
        print 'will discard multi reads'
    noUnique=False
    if '-nounique' in sys.argv:
        noUnique=True
        print 'will discard unique reads'

    PosCountsDict={}

    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))

    samfile = pysam.Samfile(SAM, "rb" )

    noNH=False
    if '-noNH' in sys.argv:
        noNH=True
        i=0
        ReadMultiplcityDict={}
        for (chr,start,end) in chromInfoList:
            try:
                jj=0
                for alignedread in samfile.fetch(chr, start, end):
                    jj+=1
                    if jj==1:
                        break
            except:
                print chr, start, end, 'not found in BAM file, skipping'
                continue
            for alignedread in samfile.fetch(chr, start, end):
                i+=1
                if i % 5000000 == 0:
                        print str(i/1000000) + 'M alignments processed in multiplicity examination', chr,start,alignedread.pos,end
                fields=str(alignedread).split('\t')
                ID=fields[0]
                if alignedread.is_read1:
                    ID = ID + '/1'
                if alignedread.is_read2:
                    ID = ID + '/2'
                if ReadMultiplcityDict.has_key(ID):
                    pass
                else:
                    ReadMultiplcityDict[ID]=0
                ReadMultiplcityDict[ID]+=1

    i=0
    for (chr,start,end) in chromInfoList:
        try:
            jj=0
            for alignedread in samfile.fetch(chr, start, end):
                jj+=1
                if jj==1:
                    break
        except:
            print chr, start, end, 'not found in BAM file, skipping'
            continue
        if PosCountsDict.has_key(chr):
            pass
        else:
            PosCountsDict[chr]={}
        for alignedread in samfile.fetch(chr, start, end):
            i+=1
            if i % 5000000 == 0:
                print str(i/1000000) + 'M alignments processed'
            fields=str(alignedread).split('\t')
            ID=fields[0]
            if alignedread.is_read1:
                ID = ID + '/1'
            if alignedread.is_read2:
                ID = ID + '/2'
            if noNH:
                scaleby = ReadMultiplcityDict[ID]
            else:
                try:
                    scaleby = alignedread.opt('NH')
                except:
                    print 'multireads not specified with the NH tag, exiting'
                    sys.exit(1)
            if noMulti and scaleby > 1:
                continue
            weight = 1.0/scaleby
            pos=alignedread.pos
            if PosCountsDict[chr].has_key(pos):
                PosCountsDict[chr][pos]+=weight
            else:
                PosCountsDict[chr][pos]=weight
            if chr == 'PDF1':
                print fields

    print '....................................'

    TotalReads = 0
    for chr in PosCountsDict.keys():
        print chr
        for pos in PosCountsDict[chr].keys():
            TotalReads += PosCountsDict[chr][pos]

    ExonicReads=0
    IntronicReads=0

    ExonicPosDict={}

    GeneDict={}
    lineslist = open(gtf)
    for line in lineslist:
        if line[0]=='#':
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        chr=fields[0]
        if GeneDict.has_key(chr):
            pass
        else:
            ExonicPosDict[chr]={}
            GeneDict[chr]={}
        start=int(fields[3])
        stop=int(fields[4])
        geneID=fields[8].split('gene_id "')[1].split('";')[0]
        if GeneDict[chr].has_key(geneID):
            pass
        else:
            GeneDict[chr][geneID]=[]
        GeneDict[chr][geneID].append((start,stop))
        for i in range(start,stop):
            ExonicPosDict[chr][i]=''

    print 'finished inputting annotation'

    keys=GeneDict.keys()
    keys.sort()

    for chr in keys:
        print chr
        for geneID in GeneDict[chr].keys():
            coordinates=[]
            for (start,stop) in GeneDict[chr][geneID]:
                coordinates.append(start)
                coordinates.append(stop)
                for i in range(start,stop):
                    ExonicPosDict[chr][i]=''
            for i in range(min(coordinates),max(coordinates)):
                if ExonicPosDict[chr].has_key(i):
                    if PosCountsDict.has_key(chr) and PosCountsDict[chr].has_key(i):
                        ExonicReads+=PosCountsDict[chr][i]
                        del PosCountsDict[chr][i]
                else:
                    if PosCountsDict.has_key(chr) and PosCountsDict[chr].has_key(i):
                        IntronicReads+=PosCountsDict[chr][i]
                        del PosCountsDict[chr][i]

    IntergenicReads = TotalReads - ExonicReads - IntronicReads

    outfile=open(outfilename,'w')
        
    outline='#Class\tFraction'
    outfile.write(outline+'\n')
    outline='Exonic:' +'\t'+str(ExonicReads/TotalReads)
    outfile.write(outline+'\n')
    outline='Intronic:' +'\t'+str(IntronicReads/TotalReads)
    outfile.write(outline+'\n')
    outline='Intergenic:' +'\t'+str(IntergenicReads/TotalReads)
    outfile.write(outline+'\n')
    outfile.close()
   
run()
