##################################
#                                #
# Last modified 2019/09/14       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
from sets import Set

def run():

    if len(sys.argv) < 4:
        print 'usage: python %s gtf upstream_bp downstream_bp outputfilename [-bioType biotype1,biotype2...biotype3] [-CDS] [-5UTR]' % sys.argv[0]
        print '\tUse the [-CDS] option if the gtf file only contains CDS entires; note that then we cannot really speak of TSSs, of course'
        sys.exit(1)
    
    GTF = sys.argv[1]
    upstreamBP = int(sys.argv[2])
    downstreamBP = int(sys.argv[3])
    outfilename = sys.argv[4]

    doCDS = False
    if '-CDS' in sys.argv:
        doCDS = True

    doBioType=False
    if '-bioType' in sys.argv:
        doBioType=True
        bioTypeDict={}
        bioTypes = sys.argv[sys.argv.index('-bioType')+1].split(',')
        print 'will only consider', bioTypes
        for b in bioTypes:
            bioTypeDict[b]=''

    do5UTR = False
    if '-5UTR' in sys.argv:
        do5UTR = True

    TranscriptDict={}
    
    linelist = open(GTF)
    i=0
    for line in linelist:
        i+=1
        if i % 100000 == 0:
            print i, 'lines processed'
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if do5UTR:
            if fields[2]!='5UTR':
                continue
        else:
            if doCDS:
                if fields[2] != 'CDS':
                    continue
            else:
                if fields[2] != 'exon':
                    continue
        if doBioType:
            if 'transcript_type "' in fields[8]:
                if bioTypeDict.has_key(fields[8].split('transcript_type "')[1].split('";')[0]):
                    pass
                else:
                    continue
            elif 'transcript_biotype "' in fields[8]:
                if bioTypeDict.has_key(fields[8].split('transcript_biotype "')[1].split('";')[0]):
                    pass
                else:
                    continue
            elif 'gene_biotype "' in fields[8]:
                if bioTypeDict.has_key(fields[8].split('gene_biotype "')[1].split('";')[0]):
                    pass
                else:
                    continue
            elif 'gene_type "' in fields[8]:
                if bioTypeDict.has_key(fields[8].split('gene_type "')[1].split('";')[0]):
                    pass
                else:
                    continue
            else:
                print 'no biotype specified, exiting'
                print fields
                sys.exit(1)
        chr=fields[0]
        start=int(fields[3])
        stop=int(fields[4])
        strand=fields[6]
        geneID=fields[8].split('gene_id "')[1].split('";')[0]
        if 'gene_name "' in fields[8]:
            geneName=fields[8].split('gene_name "')[1].split('";')[0]
        else:
            geneName=geneID
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if TranscriptDict.has_key((geneName,geneID,transcriptID)):
            pass
        else:
            TranscriptDict[(geneName,geneID,transcriptID)]=[]
        TranscriptDict[(geneName,geneID,transcriptID)].append((chr,start,stop,strand))

    TSSDict={}

    outfile = open(outfilename,'w')

    outfile.write('#chr\tleft\tright\tstrand\tgeneName(s)\tgeneID(s)\ttranscript(s)\n')

    for (geneName,geneID,transcriptID) in TranscriptDict:
        TranscriptDict[(geneName,geneID,transcriptID)].sort()
        chr = TranscriptDict[(geneName,geneID,transcriptID)][0][0]
        strand = TranscriptDict[(geneName,geneID,transcriptID)][0][3]
        if strand == '+':
            TSS=(chr,TranscriptDict[(geneName,geneID,transcriptID)][0][1],strand)
        if strand == '-':
            TSS=(chr,TranscriptDict[(geneName,geneID,transcriptID)][-1][2],strand)
        if TSSDict.has_key(TSS):
            pass
        else:
            TSSDict[TSS]={}
            TSSDict[TSS]['transcripts']=[]
            TSSDict[TSS]['genes']=[]
        TSSDict[TSS]['transcripts'].append(transcriptID)
        TSSDict[TSS]['genes'].append((geneName,geneID))
        
    keys=TSSDict.keys()
    keys.sort()
    for (chr,TSS,strand) in keys:
        TSSDict[(chr,TSS,strand)]['genes'] = list(Set(TSSDict[(chr,TSS,strand)]['genes']))
        TSSDict[(chr,TSS,strand)]['transcripts'] = list(Set(TSSDict[(chr,TSS,strand)]['transcripts']))
        if strand == '+':
            outline=chr + '\t' + str(max(0,TSS - upstreamBP)) + '\t' + str(TSS + downstreamBP) + '\t'+strand + '\t'
        if strand == '-':
            outline=chr + '\t' + str(max(0,TSS - downstreamBP)) + '\t' + str(TSS + upstreamBP) + '\t'+strand + '\t'
        for (geneName,geneID) in TSSDict[(chr,TSS,strand)]['genes']:
            outline = outline + geneName + ','
        outline=outline[0:-1] + '\t'
        for (geneName,geneID) in TSSDict[(chr,TSS,strand)]['genes']:
            outline = outline + geneID + ','
        outline=outline[0:-1] + '\t'
        for transcriptID in TSSDict[(chr,TSS,strand)]['transcripts']:
            outline = outline + transcriptID+','
        outline=outline[0:-1]
        outfile.write(outline+'\n')
   
    outfile.close()
   
run()
