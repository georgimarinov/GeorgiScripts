##################################
#                                #
# Last modified 04/05/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import numpy
from sets import Set

try:
	import psyco
	psyco.full()
except:
	pass

def run():

    if len(sys.argv) < 4:
        print 'usage: python %s gtf wig minGeneLength outfile [-field1 biotype] [-normalize] [-maxGeneLength bp] [-singlemodelgenes] [-printlist]' % sys.argv[0]
        sys.exit(1)

    gtf = sys.argv[1]
    wig = sys.argv[2]
    minGeneLength = int(sys.argv[3])
    outputfilename = sys.argv[4]

    doPrintList = False
    if '-printlist' in sys.argv:
        doPrintList = True

    doBioType=False
    BioType=''
    if '-field1' in sys.argv:
        doBioType=True
        BioType=sys.argv[sys.argv.index('-field1')+1]
        print 'will only consider', BioType, 'genes'

    doMaxGeneLength=False
    if '-maxGeneLength' in sys.argv:
        doMaxGeneLength=True
        maxGeneLength=int(sys.argv[sys.argv.index('-maxGeneLength')+1])
        print 'will only consider genes longer than', minGeneLength, 'and shorter than', maxGeneLength

    doNormalize=False
    if '-normalize' in sys.argv:
        doNormalize=True
        print 'will normalize scores'

    doSingleModel=False
    if '-singlemodelgenes' in sys.argv:
        doSingleModel=True
        print 'will only use genes with one isoform'

    listoflines = open(gtf)
    GeneDict={}
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        if doBioType:
            if fields[1] != BioType:
                continue
        chr=fields[0]
        strand=fields[6]
        left=int(fields[3])
        right=int(fields[4])
        geneID=fields[8].split('gene_id "')[1].split('";')[0]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if GeneDict.has_key(geneID):
            pass
        else:
            GeneDict[geneID]={}
        if GeneDict[geneID].has_key(transcriptID):
            if chr != GeneDict[geneID][transcriptID][0][0]:
                continue
        else:
            GeneDict[geneID][transcriptID]=[]
        GeneDict[geneID][transcriptID].append((chr,left,right,strand))

    print 'finished inputting annotation'

    CoverageDict={}

    i=0
    for geneID in GeneDict.keys():
        if doSingleModel and len(GeneDict[geneID].keys()) > 1:
            del GeneDict[geneID]
            continue
        i+=1
        for transcriptID in GeneDict[geneID].keys():
            for (chr,left,right,strand) in GeneDict[geneID][transcriptID]:
                if CoverageDict.has_key(chr):
                    pass
                else:
                    CoverageDict[chr]={}
                for j in range(left,right):
                    CoverageDict[chr][j]=0

    listoflines = open(wig)
    for line in listoflines:
        if line.startswith('track'):
            continue
        if line.startswith('#'):
            continue
        fields=line.replace(' ','\t').strip().split('\t')
        chr=fields[0]
        left=int(fields[1])
        right=int(fields[2])
        score=float(fields[3])
        if CoverageDict.has_key(chr):
            pass
        else:
            continue
        for j in range(left,right):
            if CoverageDict[chr].has_key(j):
                CoverageDict[chr][j]=score

    print 'finished inputting wiggle'

    output_Array={}
    for i in range(100):
        output_Array[i]=0

    print len(GeneDict.keys())

    if doPrintList:
        outfile=open(outputfilename + '.geneList','w')

    geneNumber=0.0
    for geneID in GeneDict.keys():
        NucleotideList=[]
        for transcriptID in GeneDict[geneID].keys():
            for (chr,left,right,strand) in GeneDict[geneID][transcriptID]:
                for i in range(left,right):
                    NucleotideList.append(i)
        NucleotideList=list(Set(NucleotideList))
        NucleotideList.sort()
        if strand=='-' or strand=='R':
            NucleotideList.reverse()
        geneLength=len(NucleotideList)
        if geneLength < minGeneLength:
            continue
        if doMaxGeneLength:
            if geneLength > maxGeneLength:
                continue
        stepsize = geneLength/100.0
        k=0
        final_vector=[]
        while k < geneLength-stepsize:
            b=k+stepsize
            counts=[]
            for i in range(int(k),int(b)):
                counts.append(CoverageDict[chr][NucleotideList[i]])
            final_vector.append(numpy.mean(counts))
            k=b
        i=0
        for v in final_vector:
            output_Array[i]+=v
            i+=1
        geneNumber+=1
        if doPrintList:
            outfile.write(geneID + '\n')

    if doPrintList:
        outfile.close()

    print geneNumber, 'genes considered'

    outfile=open(outputfilename,'w')

    for i in range(100):
        if doNormalize:
            outfile.write(str(i) + '\t' + str(output_Array[i]/geneNumber)+'\n')
        else:
            outfile.write(str(i) + '\t' + str(output_Array[i])+'\n')

    outfile.close()

run()
