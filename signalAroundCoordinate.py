##################################
#                                #
# Last modified 2018/11/11       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
import os
import gzip
import pyBigWig
from sets import Set

def run():

    if len(sys.argv) < 7:
        print 'usage: python %s inputfilename chrFieldID posField strandField radius bigWig outputfilename [-normalize] [-bed] [-unstranded] [-ERANGE_hts] [-narrowPeak] [-first number]' % sys.argv[0]
        print '\tInput format: <fields .. tabs> chr <tab> position <tab> strandField' 
        print '\tThis script outputs the average signal over all regions within the given radius' 
        print '\tif the the -bed option is used, the middle point of a bed region will be used; specifiy the posField as the left coordinate of the region' 
        print '\tif the the -narrowPeak option is used, the posFielf will be ignored and strand will be assumed to be +' 
        sys.exit(1)
    
    regionfilename = sys.argv[1]
    chrFieldID = int(sys.argv[2])
    posFieldID = int(sys.argv[3])
    strandFieldID = int(sys.argv[4])
    radius = int(sys.argv[5])
    bigWig = sys.argv[6]
    outfilename = sys.argv[7]

    doFirst=False
    if '-first' in sys.argv:
        firstN=int(sys.argv[sys.argv.index('-first')+1])
        doFirst=True
        print 'will only look at the first', firstN, 'locations'

    doNarrowPeak=False
    if '-narrowPeak' in sys.argv:
        doNarrowPeak=True
        print 'will treat regions as being in narrowPeak format'

    doBed=False
    if '-bed' in sys.argv:
        print 'will treat input as bed file and center around the midpoint of reigons'
        doBed=True

    noStrand = False
    if '-unstranded' in sys.argv:
        print 'will treat all regions as + strand'
        noStrand = True

    doERANGE=False
    if '-ERANGE_hts' in sys.argv:
        doERANGE=True

    doNormalize=False
    if '-normalize' in sys.argv:
        doNormalize=True

    bw = pyBigWig.open(bigWig)

    RegionDict={}
    ScoreDict={}
    if regionfilename.endswith('.gz'):
        listoflines = gzip.open(regionfilename)
    else:
        listoflines = open(regionfilename)
    k=0
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields=line.replace('\x00','').strip().split('\t')
        if doNarrowPeak:
            pass
        else:
            if len(fields) < max(chrFieldID, posFieldID, strandFieldID, 3):
                continue
        k+=1
        if doFirst and k > firstN:
            continue
        if len(fields)<3:
           continue
        if doBed:
            chr=fields[chrFieldID]
            left=int(fields[posFieldID])
            right=int(fields[posFieldID+1])
            pos=int((right+left)/2.0)
            if noStrand:
                strand='+'
            else:
                strand=fields[strandFieldID]
        elif doERANGE:
            chr=fields[1]
            pos=int(fields[9])
            strand='+'
        elif doNarrowPeak:
            chr=fields[0]
            pos=int(fields[1]) + int(fields[9])
            strand='+'
        else:
            chr=fields[chrFieldID]
            pos=int(fields[posFieldID])
            if noStrand:
                strand='+'
            else:
                strand=fields[strandFieldID]
        RegionDict[(chr,pos,strand)]=[]
        if ScoreDict.has_key(chr):
            pass
        else:
            ScoreDict[chr]={}
        try:
            RegionDict[(chr,pos,strand)] = bw.values(chr,pos-radius,pos+radius)
        except:
            print 'problem with region', chr, pos, strand, 'skipping'
            del RegionDict[(chr,pos,strand)]
            continue
        if strand == 'R' or strand == '-':
            RegionDict[(chr,pos,strand)].reverse()
        if k % 10000 == 0:
            print k

    keys=RegionDict.keys()
    keys.sort()

    FinalDict={}
    for i in range(0-radius,0+radius):
        FinalDict[i]=0.0
    for (chr,pos,strand) in keys:
        for i in range(-radius,+radius):
            if math.isnan(RegionDict[(chr,pos,strand)][i+radius]):
                pass
            else:
                FinalDict[i] += RegionDict[(chr,pos,strand)][i+radius]

    outfile=open(outfilename,'w')
    
    keys=FinalDict.keys()
    keys.sort()
    for i in keys:
        outline = str(i) + '\t' + str(FinalDict[i])
        if doNormalize:
            outline = str(i) + '\t' + str(FinalDict[i]/len(RegionDict))
        outfile.write(outline + '\n')

    outfile.close()
   
run()
