##################################
#                                #
# Last modified 09/19/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
from sets import Set

try:
    import psyco
    psyco.full()
except:
    print 'psyco not running'

def run():

    if len(sys.argv) < 5:
        print 'usage: python %s file1 chromField1 file2 chromField2 outfilenameprefix [-minOverlap fraction-of-2nd-file-region] [-minOverlapBP bp] [-combinedOutput]' % sys.argv[0]
        print "       Note: have file1 be the smaller file; the script will store its coordinates in the memory in a dictionary and check against the other file as it goes over its entries"
        print "       Note: works with 0bp-sized regions"
        print "       Note: this script does not work with overlapping regions in the first file!!!"
        print "       Note: if you use the minimal overlap option, only the regions from the second files will be intersected under that criteria, i.e. what is in the intersection2 files; regions in intersection1 will be intersected using the 1bp overlap crtieria"
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[3]
    outfileprefix = sys.argv[5]
    minOverlap=0
    if '-minOverlap' in sys.argv:
        minOverlap=float(sys.argv[sys.argv.index('-minOverlap') + 1])
        print 'will require overlap of', minOverlap, 'fraction of entries in the second file'
    chromField1 = int(sys.argv[2])
    chromField2 = int(sys.argv[4])

    doMinBasePairOverlap = False
    if '-minOverlapBP' in sys.argv:
        doMinBasePairOverlap = True
        minOverlapBP = int(sys.argv[sys.argv.index('-minOverlapBP') + 1])
        print 'will require overlap of', minOverlapBP, 'bp of entries in the second file'

    File1DictCoverage={}
    File1Dict={}

    p=0
    listoflines = open(file1)
    i=0
    for line in listoflines:
        i+=1
        if line.startswith('#') or line.startswith('track ') or line.strip() == '':
            continue
        if len(line.strip())==0:
            continue
        fields=line.strip().split('\t')
        chr=fields[chromField1]
        if File1DictCoverage.has_key(chr):
            pass
        else:
            File1DictCoverage[chr]={}
        start = int(fields[chromField1+1])
        end = int(fields[chromField1+2])
        if end < start:
            print 'region of negative length, skipping'
            print line.strip()
            continue
        for j in range(start,end):
            File1DictCoverage[chr][j]=i
            p+=1
        if end == start:
            File1DictCoverage[chr][start]=i
        File1Dict[i]=line

    print 'total bases covered:', p

    OverlappedListDict={}

    outfilename_intersection1=open(outfileprefix+'-intersection1','w')
    outfilename_intersection2=open(outfileprefix+'-intersection2','w')
    outfilename_outersection1=open(outfileprefix+'-outersection1','w')
    outfilename_outersection2=open(outfileprefix+'-outersection2','w')

    listoflines = open(file2)
    k=0
    for line in listoflines:
        k+=1
        if k % 100000 == 0:
            print k, 'lines in file2 processed'
        if line.startswith('#') or line.startswith('track ') or line.strip() == '':
            continue
        fields=line.strip().split('\t')
        chr=fields[chromField2].split(':')[0]
        start = int(fields[chromField2+1])
        end = int(fields[chromField2+2])
        if File1DictCoverage.has_key(chr):
            pass
        else:
            outfilename_outersection2.write(line)
            continue
        overlapBP=0
        for i in range(start,end):
            if File1DictCoverage[chr].has_key(i):
                OverlappedListDict[File1DictCoverage[chr][i]]=0
                overlapBP+=1
        if doMinBasePairOverlap:
            if overlapBP >= minOverlapBP:
                outfilename_intersection2.write(line)
            else:
                outfilename_outersection2.write(line)
        else:
            if overlapBP != 0 and (overlapBP/(end - start + 0.0) >= minOverlap):
                outfilename_intersection2.write(line)
            else:
                outfilename_outersection2.write(line)

    for i in File1Dict.keys():
        if OverlappedListDict.has_key(i):
            outfilename_intersection1.write(File1Dict[i])
        else:
            outfilename_outersection1.write(File1Dict[i])

    outfilename_intersection1.close()
    outfilename_intersection2.close()    
    outfilename_outersection1.close()
    outfilename_outersection2.close()

run()
