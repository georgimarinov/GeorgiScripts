##################################
#                                #
# Last modified 2019/03/22       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
import os
import pyBigWig
from sets import Set

def run():

    if len(sys.argv) < 6:
        print 'usage: python %s inputfilename chrFieldID posField [strandField | -noStrand] upstream downstream bigWig outputfilename [-strand +|-] [-average bp] [-window bp] [-sortby fieldID] [-fullRegionInfo] [-narrowPeak]' % sys.argv[0]
        print '\tInput format: <fields .. tabs> chr <tab> position' 
        print '\tthe wig file can be in .bz2 or .gz format' 
        sys.exit(1)
    
    regionfilename = sys.argv[1]
    chrFieldID = int(sys.argv[2])
    posFieldID = int(sys.argv[3])
    noStrand=False
    if sys.argv[4]=='-noStrand':
        noStrand=True
    else:
        strandFieldID = int(sys.argv[4])
    upstream = int(sys.argv[5])
    downstream = int(sys.argv[6])
    bigWig = sys.argv[7]
    outfilename = sys.argv[8]

    doNP = False
    if '-narrowPeak' in sys.argv:
        doNP = True

    doStrand = False
    if '-strand' in sys.argv:
        doStrand = True
        WantedStrand = sys.argv[sys.argv.index('-strand') + 1]

    sortFieldID = posFieldID
    sortList = []
    if '-sortby' in sys.argv:
        sortFieldID = int(sys.argv[sys.argv.index('-sortby')+1])

    doFullRegionInfo = False
    if '-fullRegionInfo' in sys.argv:
        doFullRegionInfo = True

    window=1
    averageRadius = 0
    doAverage=False
    if '-average' in sys.argv:
        doAverage=True
        averageRadius=int(int(sys.argv[sys.argv.index('-average')+1])/2.0)
        print 'will average signal over', 2*averageRadius, 'bp'

    doWindow=False
    if '-window' in sys.argv:
        doWindow=True
        window=int(sys.argv[sys.argv.index('-window')+1])
        print 'will split into windows of size', window, 'bp'

    bw = pyBigWig.open(bigWig)

    RegionDict={}
    ScoreDict={}
    header=''
    if regionfilename.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + regionfilename
    elif regionfilename.endswith('.gz'):
        cmd = 'zcat ' + regionfilename
    else:
        cmd = 'cat ' + regionfilename
    p = os.popen(cmd, "r")
    LC = 0
    line = 'line'
    while line != '':
        line = p.readline()
        if line == '':
            break
        LC += 1
        if LC % 10000 == 0:
            print LC, 'lines processed'
        if line.startswith('#'):
            header=line.strip()
            continue
        fields = line.strip().split('\t')
        chr = fields[chrFieldID]
        if doNP:
            pos = int(fields[1]) + int(fields[9])
        else:
            pos = int(fields[posFieldID])
        if noStrand:
            strand = '+'
        else:
            strand = fields[strandFieldID]
        if doStrand and strand != WantedStrand:
            continue
        RegionDict[(chr,pos,strand)]={}
        RegionDict[(chr,pos,strand)]['line'] = line.strip()
        if ScoreDict.has_key(chr):
            pass
        else:
            ScoreDict[chr]={}
        if strand=='+' or strand=='F':
            try:
                values = bw.values(chr,pos-upstream-averageRadius,pos+downstream+averageRadius)
            except:
                print 'skipping region:'
                print fields, pos, upstream, downstream
                print chr, pos, pos-upstream-averageRadius, pos+downstream+averageRadius
                continue
        if strand=='-' or strand=='R':
            try:
                values = bw.values(chr,pos-downstream-averageRadius,pos+upstream+averageRadius)
            except:
                print 'skipping region:'
                print fields
                print chr, pos, pos-upstream-averageRadius, pos+downstream+averageRadius
                continue
            values.reverse()
        for i in range(len(values)):
            if math.isnan(values[i]):
                values[i] = 0
        sortList.append((float(fields[sortFieldID]),chr,pos,strand,values))

    outfile=open(outfilename,'w')
    if doFullRegionInfo:
        outline='#chr\tleft\tright\tstrand'
    else:
        outline='#'
    for i in range(0-upstream,0+downstream+1,window):
        outline=outline+'\t'+str(i)
    outfile.write(outline+'\n')

    sortList.sort()
    sortList.reverse()
    for (score,chr,pos,strand,values) in sortList:
        FinalDict = {}
        if doWindow:
            for i in range(0 - upstream, 0 + downstream + 1, window):
                FinalDict[i]=0
        elif doAverage:
            for i in range(0 - upstream - averageRadius, 0 + downstream + averageRadius + 1):
                FinalDict[i]=0
        else:
            for i in range(0 - upstream, 0 + downstream + 1):
                FinalDict[i]=0
        ks = FinalDict.keys()
        ks.sort()
#        print ks
        if doAverage:
            for i in range(averageRadius,len(values) - averageRadius):
#                print i, i-upstream, averageRadius, len(values)
                FinalDict[i - upstream] += sum(values[i - averageRadius:i + averageRadius])/(2.0*averageRadius)
        elif doWindow:
            for i in range(0,len(values)+1,window):
                FinalDict[i - upstream] += sum(values[i:i + window])/(window + 0.0)
        else:
            for i in range(len(values)):
                FinalDict[i - upstream] += values[i]
        if doFullRegionInfo:
            outline = chr + '\t' + str(pos-upstream) + '\t' + str(pos+downstream) + '\t'+strand
        else:
            outline = chr + '|' + str(pos) + '|' + strand
        keys = FinalDict.keys()
        keys.sort()
        for i in keys:
            outline = outline + '\t ' +str(FinalDict[i])
        outfile.write(outline + '\n')

    outfile.close()
   
run()
