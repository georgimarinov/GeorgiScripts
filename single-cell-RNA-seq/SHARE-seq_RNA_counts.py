##################################
#                                #
# Last modified 2021/04/26       # 
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

    if len(sys.argv) < 4:
        print 'usage: python %s read1_BAMfilename GTF chrom.sizes outprefix [-UMIedit N] [-UMICutoff N_UMIs_per_cell] [-exonicOnly] [-uniqueBAM] [-noNH samtools] [-noNHinfo] [-medianExpCellsCutoff fraction] [-3UTRextend bp] [-refFlat] [-splitAmbiguousReads]' % sys.argv[0]
        print '\t the script will only use uniquely mapped reads'
        print '\t the [-medianExpCellsCutoff] option will remove all cells with UMIs fewer than the indicated fraction of the median for the top N cells, where N is the expected number of cells (default: 0.20)'
        print '\t the [-3UTRextend] option will extend the 3UTRs of genes by the indicated number of bp. Default: 100'
        print '\t #### [-exonicOnly] option not implemented yet'
        print '\t The [-splitAmbiguousReads] will give each overlapping gene that a uniquely mapped read could be assigned a fractional UMI count'
        print '\t The default [-UMICutoff] is 0'
        print '\t UMIs with more than one N will be excluded'
        sys.exit(1)

    BAM = sys.argv[1]
    GTF = sys.argv[2]
    chrominfo=sys.argv[3]
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
    outprefix = sys.argv[4]

    noMulti=True
#    if '-nomulti' in sys.argv:
#        print 'will only consider unique alignments'
#        noMulti=True

    doExonicOnly = False
    if '-exonicOnly' in sys.argv:
        doExonicOnly = True

    UMIedit = 1
    if '-UMIedit' in sys.argv:
        UMIedit = int(sys.argv[sys.argv.index('-UMIedit') + 1])

    UMICutoff = 0
    if '-UMICutoff' in sys.argv:
        UMICutoff = float(sys.argv[sys.argv.index('-UMICutoff') + 1])

    ReadIDDict = {}

    samfile = pysam.Samfile(BAM, "rb" )

    UTR3ext = 100
    if '-3UTRextend' in sys.argv:
        UTR3ext = int(sys.argv[sys.argv.index('-3UTRextend') + 1])

    doUniqueBAM = False
    if '-uniqueBAM' in sys.argv:
        print 'will treat all alignments as unique'
        doUniqueBAM = True
        TotalReads = 0
        pass

    doRF = False
    if '-refFlat' in sys.argv:
        print 'will treat GTF file as if it is a refFlat one'
        doRF = True

    doSAR = False
    if '-splitAmbiguousReads' in sys.argv:
        print 'will assign fractional weights to genes with amniguous alignments'
        doSAR = True

    doNoNHinfo = False
    if '-noNHinfo' in sys.argv:
        doNoNHinfo = True
        print 'will directly evaluate read mulitplicity'
        MultiplicityDict = {}
    else:
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
#                if alignedread.is_read1:
#                    ID = ID + '/1'
#                if alignedread.is_read2:
#                    ID = ID + '/2'
                if MultiplicityDict.has_key(ID):
                    pass
                else:
                    MultiplicityDict[ID] = 0
                MultiplicityDict[ID] += 1

    j=0
    lineslist = open(GTF)
    TranscriptDict = {}
    GeneDict = {}
    GeneIDDict = {}
    GeneNameDict = {}
    G = 0
    for line in lineslist:
        j+=1
        if j % 1000000 == 0:
            print j, 'lines processed'
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        chr = fields[0]
        if 'gene_name "' in fields[8]:
            geneName = fields[8].split('gene_name "')[1].split('";')[0]
        else:
            geneName = fields[8].split('gene_id "')[1].split('";')[0]
        geneID = fields[8].split('gene_id "')[1].split('";')[0]
        if doRF:
            geneID = geneID + '-' + chr
        if 'transcript_name "' in fields[8]:
            transcriptName = fields[8].split('transcript_name "')[1].split('";')[0]
        else:
            transcriptName = fields[8].split('transcript_id "')[1].split('";')[0]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        transcript = (geneID, geneName, transcriptName, transcriptID)
        if GeneDict.has_key(geneID):
            pass
        else:
            GeneDict[geneID] = {}
            GeneDict[geneID]['transcripts'] = {}
            GeneDict[geneID]['transcripts'][transcript] = 1
            GeneDict[geneID]['coordinates'] = []
            G+=1
            GeneDict[geneID]['G'] = G
            GeneIDDict[G] = geneID
            GeneNameDict[G] = geneName
        if TranscriptDict.has_key(transcript):
            pass
        else:
            TranscriptDict[transcript]=[]
        left = int(fields[3])
        right = int(fields[4])
        strand = fields[6]
        TranscriptDict[transcript].append((chr,left,right,strand))
        GeneDict[geneID]['chr'] = chr
        GeneDict[geneID]['strand'] = strand
        GeneDict[geneID]['coordinates'].append(left)
        GeneDict[geneID]['coordinates'].append(right)

    for transcript in TranscriptDict.keys():
        TranscriptDict[transcript].sort()

    print 'finished parsing GTF file'

    print 'found', len(GeneDict.keys()), 'genes'

    ReadDict = {}

    RN=0
    GN = 0
    for geneID in GeneDict.keys():
        GN += 1
        chr = GeneDict[geneID]['chr'] 
        start = min(GeneDict[geneID]['coordinates'])
        end = max(GeneDict[geneID]['coordinates'])
        strand = GeneDict[geneID]['strand']
        if strand == '+':
            end = min(end + UTR3ext, chromInfoDict[chr])
        if strand == '-':
            start = max(0,start - UTR3ext)
        for alignedread in samfile.fetch(chr, start, end):
            RN+=1
            if RN % 1000000 == 0:
                print str(RN/1000000) + 'M alignments processed;', GN, 'genes processed'
            fields=str(alignedread).split('\t')
            FLAGfields = FLAG(int(fields[1]))
            if alignedread.is_reverse:
                s = '-'
            else:
                s = '+'
            if s != strand:
                continue
            ID = fields[0].split('_2:N:0:')[0]
            if doUniqueBAM:
                multiplicity = 1
            elif doNoNHinfo:
                multiplicity = MultiplicityDict[ID]
            else:
                multiplicity = alignedread.opt('NH')
            if multiplicity > 1:
                continue
            cigar = alignedread.cigar
            pos = alignedread.pos
            FivePEnd = alignedread.pos
            if strand == '-':
                currentPos = alignedread.pos
                for (m,bp) in cigar:
                    currentPos = currentPos + bp
                FivePEnd = currentPos
            if FivePEnd < start or FivePEnd > end:
                continue
            if ReadDict.has_key(ID):
                pass
            else:
                ReadDict[ID] = {}
                ReadDict[ID]['al'] = []
                ReadDict[ID]['gs'] = []
            ReadDict[ID]['al'] = (chr,pos,s,str(cigar))
            ReadDict[ID]['gs'].append(GeneDict[geneID]['G'])

    print 'found alignments for', len(ReadDict.keys()), 'reads'
    GeneExpressionMatrix = {}

    SeenGenes = {}

    i = 0
    ambiguous = 0
    for ID in ReadDict.keys():
        i += 1
        if i % 1000000 == 0:
            print i
        [BC1,BC2,BC3,UMI] = ID.split(':::[')[-1].split(']')[0].split('+')
        if UMI.count('N') > 1:
            continue
        if BC1 == 'nan':
            continue
        if BC2 == 'nan':
            continue
        if BC3 == 'nan':
            continue
        BC = (BC1,BC2,BC3)
        AL = ReadDict[ID]['al']
        genes = ReadDict[ID]['gs']
        if len(genes) == 1:
            TheG = [genes[0]]
        else:
            (chr, pos, strand, cigar) = AL
            cigarfields = cigar.split('), (')
            if len(cigarfields) == 1:
                ambiguous += 1
                if doSAR:
                    TheG = genes
                else:
                    continue
            else:
                splicesites = {}
                KCF = 0
                for CF in cigarfields:
                    KCF+=1
                    CF = CF.replace('[','').replace(']','').replace('(','').replace(')','')
                    m = int(CF.split(', ')[0])
                    bp = int(CF.split(', ')[1])
                    if m == 0 and KCF > 1:
                        splicesites[currentPos] = 1
                    elif m == 3:
                        if strand == '+':
                            splicesites[currentPos+1] = 1
                        if strand == '-':
                            splicesites[currentPos+1] = 1
                    else:
                        pass
                    currentPos = currentPos+bp
                matchedSS = {}
                for G in genes:
                    matchedSS[G] = {}
                    geneID = GeneIDDict[G]
                    for transcript in GeneDict[geneID]['transcripts']:
                        KTP = 0
                        for (chr,LL,RR,strand) in TranscriptDict[transcript]:
                            KTP += 1
                            if KTP == 1:
                                if splicesites.has_key(RR):
                                    matchedSS[G][RR] = 1
                            elif KTP == len(TranscriptDict[transcript]):
                                if splicesites.has_key(LL):
                                    matchedSS[G][LL] = 1
                            else:
                                if splicesites.has_key(RR):
                                    matchedSS[G][RR] = 1
                                if splicesites.has_key(LL):
                                    matchedSS[G][LL] = 1
                matchedSSgenes = []
                for G in genes:
                    matchedSSgenes.append((len(matchedSS[G].keys()),G))
                matchedSSgenes.sort()
                matchedSSgenes.reverse()
                if matchedSSgenes[0][0] > matchedSSgenes[1][0]:
                    TheG = [matchedSSgenes[0][0]]
                else:
                    ambiguous += 1
                    if doSAR:
                        TheG = []
                        for G in matchedSSgenes:
                            if G[0] == matchedSSgenes[0]:
                                TheG.append(G[1])
                    else:
                        continue
        if GeneExpressionMatrix.has_key(BC):
            pass
        else:
            GeneExpressionMatrix[BC] = {}
            GeneExpressionMatrix[BC]['UMIs'] = {}
            GeneExpressionMatrix[BC]['gene_counts'] = {}
        if GeneExpressionMatrix[BC]['UMIs'].has_key(UMI):
            GeneExpressionMatrix[BC]['UMIs'][UMI].append((len(TheG),tuple(TheG)))
            continue
        else:
            EDist = len(UMI)
            NearestUMI = []
            for U in GeneExpressionMatrix[BC]['UMIs'].keys():
                LDist = Levenshtein.distance(U,UMI)
                if LDist <= UMIedit:
                    if LDist < EDist:
                        EDist = LDist
                        NearestUMI = [U]
                    if LDist == EDist:
                        NearestUMI.append(U)
            NearestUMI = list(Set(NearestUMI))
            if len(NearestUMI) == 0:
#                if UMI.count('N') == 0:
#                    GeneExpressionMatrix[BC]['UMIs'][UMI] = []
                GeneExpressionMatrix[BC]['UMIs'][UMI] = []
                GeneExpressionMatrix[BC]['UMIs'][UMI].append((len(TheG),tuple(TheG)))
#            else len(NearestUMI) == 1:
            else:
                GeneExpressionMatrix[BC]['UMIs'][NearestUMI[0]].append((len(TheG),tuple(TheG)))
#            else:
#                print 'multiple matching UMIs detected for the following read:'
#                print AL
#                print NearestUMI
#               print BC

    print 'finished collapsing UMIs'
    if doSAR:
        print 'found', ambiguous, 'reads with ambiguous gene assignment'
    else:
        print 'discarded', ambiguous, 'reads with ambiguous gene assignment'  
    print 'found', len(GeneExpressionMatrix.keys()), 'cell barcodes in total'

    print 'starting dismbiguating UMIs'

    for BC in GeneExpressionMatrix.keys():
        for UMI in GeneExpressionMatrix[BC]['UMIs']:
            TheG = []
            GeneExpressionMatrix[BC]['UMIs'][UMI].sort()
            genes = list(Set(GeneExpressionMatrix[BC]['UMIs'][UMI]))
            if len(genes) == 1:
                G = GeneExpressionMatrix[BC]['UMIs'][UMI][0][1]
                if GeneExpressionMatrix[BC]['UMIs'][UMI][0][0] > 1:
                    G = list(G)
                    for GG in G:
                        TheG.append(GG)
                else:
                    TheG.append(G[0])
            else:
                i = 0
                while i < len(GeneExpressionMatrix[BC]['UMIs'][UMI]) - 2:
#                    print i, BC, UMI, GeneExpressionMatrix[BC]['UMIs'][UMI][i]
#                    print GeneNameDict[GeneExpressionMatrix[BC]['UMIs'][UMI][i][1][0]]
                    for G in GeneExpressionMatrix[BC]['UMIs'][UMI][i][1]:
                        if GeneExpressionMatrix[BC]['UMIs'][UMI][i][0] > 1:
                            G = list(G)
                            for GG in G:
                                TheG.append(GG)
                        else:
                            TheG.append(G)
                    if GeneExpressionMatrix[BC]['UMIs'][UMI][i][0] < GeneExpressionMatrix[BC]['UMIs'][UMI][i+1][0]:
                        break
                    i+=1
            TheG = list(Set(TheG))
            for G in TheG:
                SeenGenes[G] = 1
                if GeneExpressionMatrix[BC]['gene_counts'].has_key(G):
                    if len(TheG) == 1:
                        GeneExpressionMatrix[BC]['gene_counts'][G] +=1
                    else:
                        GeneExpressionMatrix[BC]['gene_counts'][G] += 1./len(TheG)
                else:
                    if len(TheG) == 1:
                        GeneExpressionMatrix[BC]['gene_counts'][G] = 1
                    else:
                        GeneExpressionMatrix[BC]['gene_counts'][G] = 1./len(TheG)

    print 'finished dismbiguating UMIs'

    outfile = open(outprefix + '.UMIs_per_cell', 'w')
    outline = '#BC1+BC2+BC3\trank\tUMIs\tAligned_Positions\tgenes'
    outfile.write(outline + '\n')

    outlines = []

    for BC in GeneExpressionMatrix.keys():
        UMIs = len(GeneExpressionMatrix[BC]['UMIs'].keys())
        ALs = 0
        for UMI in GeneExpressionMatrix[BC]['UMIs'].keys():
            ALs += len(GeneExpressionMatrix[BC]['UMIs'][UMI])
        outlines.append((UMIs,ALs,BC))

    outlines.sort()
    outlines.reverse()

    R=0
    for (UMIs,ALs,BC) in outlines:
        R+=1
        outline = BC[0] + '+' + BC[1] + '+' + BC[2] + '\t' + str(R) + '\t' + str(UMIs) + '\t' + str(ALs) + '\t' + str(len(GeneExpressionMatrix[BC]['gene_counts'].keys()))
        outfile.write(outline + '\n')

    outfile.close()

    print 'filtering poor quality cells'

    print 'cells before filtering poor quality cells:', len(GeneExpressionMatrix.keys())

    for BC in GeneExpressionMatrix.keys():
        if len(GeneExpressionMatrix[BC]['UMIs'].keys()) < UMICutoff:
            del GeneExpressionMatrix[BC]

    print 'cells remaining after filtering poor quality cells:', len(GeneExpressionMatrix.keys())

    outfile = open(outprefix + '.table', 'w')

    outline = '#geneID\tgeneName'
    BCs = GeneExpressionMatrix.keys()
    BCs.sort()
    for BC in BCs:
        (BC1,BC2,BC3) = BC
        outline = outline + '\t' + BC1 + '+' + BC2 + '+' + BC3
    outfile.write(outline+'\n')

    Gs = SeenGenes.keys()
    for G in Gs:
        geneID = GeneIDDict[G]
        geneName = GeneNameDict[G]
        outline = geneID + '\t' + geneName
        for BC in BCs:
            if GeneExpressionMatrix[BC]['gene_counts'].has_key(G):
                umicounts = GeneExpressionMatrix[BC]['gene_counts'][G]
            else:
                umicounts = 0
            outline = outline + '\t' + str(umicounts)
        outfile.write(outline+'\n')

    print 'finished outputting counts'

    outfile.close()
            
run()
