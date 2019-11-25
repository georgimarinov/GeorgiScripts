##################################
#                                #
# Last modified 2017/07/15       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import os
import random

try:
	import psyco
	psyco.full()
except:
	pass

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s end1 end2 [-renameReads prefix] [-shuffle] [-trim end1length end2legnth] [-trim5 bptoberemoved1 bptoberemoved2] [-minReadLen bp]' % sys.argv[0]
        print '\tThe script can read compressed files as long as they have the correct suffix - .zip, .bz2 or .gz'
        print '\tYou can have multiple files in each end1 or end2 argument, they have to be comma-separated'
        print '\tit will print to standard output'
        sys.exit(1)

    fastq1files = sys.argv[1].split(',')
    fastq2files = sys.argv[2].split(',')

    doShuffle = False
    if '-shuffle' in sys.argv:
        doShuffle = True
        DataList = []

    doRename = False
    if '-renameReads' in sys.argv:
        doRename = True
        prefix = sys.argv[sys.argv.index('-renameReads') + 1]

    if len(fastq1files) != len(fastq1files):
        print 'incorrect number of input files'
        sys.exit(1)

    minRL = 0
    if '-minReadLen' in sys.argv:
        minRL = int(sys.argv[sys.argv.index('-minReadLen') + 1])
    
    doTrim = False
    if '-trim' in sys.argv:
        doTrim = True
        end1Trim = int(sys.argv[sys.argv.index('-trim') + 1])
        end2Trim = int(sys.argv[sys.argv.index('-trim') + 2])

    doTrim5 = False
    if '-trim5' in sys.argv:
        doTrim5 = True
        end1Trim5 = int(sys.argv[sys.argv.index('-trim5') + 1])
        end2Trim5 = int(sys.argv[sys.argv.index('-trim5') + 2])

    i=0
    for f in range(len(fastq1files)):
        fastq1 = fastq1files[f]
        fastq2 = fastq2files[f]
        if fastq1.endswith('.bz2'):
            cmd1 = 'bzip2 -cd ' + fastq1
        elif fastq1.endswith('.gz'):
            cmd1 = 'zcat ' + fastq1
        elif fastq1.endswith('.zip'):
            cmd1 = 'unzip -p ' + fastq1
        else:
            cmd1 = 'cat ' + fastq1
        p1 = os.popen(cmd1, "r")

        if fastq2.endswith('.bz2'):
            cmd2 = 'bzip2 -cd ' + fastq2
        elif fastq2.endswith('.gz'):
            cmd2 = 'gunzip -c ' + fastq2
        elif fastq2.endswith('.zip'):
            cmd1 = 'unzip -p ' + fastq2
        else:
            cmd2 = 'cat ' + fastq2
        p2 = os.popen(cmd2, "r")

        line = 'line'
        lines = []
        while line != '':
            line = p1.readline().strip()
            lines.append(line)
            if line == '':
                break
            i+=1
            if i % 4 == 0:
               ID1 = lines[0]
               seq1 = lines[1]
               q1 = lines[3]
               ID2 = p2.readline().strip()
               seq2 = p2.readline().strip()
               ID22 = p2.readline().strip()
               q2 = p2.readline().strip()
               if doTrim5:
                   seq1 = seq1[end1Trim5:]
                   q1 = q1[end1Trim5:]
                   seq2 = seq2[end1Trim5:]
                   q2 = q2[end1Trim5:]
               if doTrim:
                   seq1 = seq1[0:end1Trim]
                   q1 = q1[0:end1Trim]
                   seq2 = seq2[0:end2Trim]
                   q2 = q2[0:end2Trim]
               if doRename:
                   outline = prefix + str(i/4) + '\t' + seq1 + '\t' + q1 + '\t' + seq2 + '\t' + q2
               else:
#                   outline = ID1[1:len(ID1)].split('/')[0].split(' ')[0] + '\t' + seq1 + '\t' + q1 + '\t' + seq2 + '\t' + q2
                   outline = ID1[1:len(ID1)].split('/')[0] + '\t' + seq1 + '\t' + q1 + '\t' + seq2 + '\t' + q2
               if len(seq1) >= minRL and len(seq2) >= minRL:
                   if doShuffle:
                       DataList.append(outline)
                   else:
                       print outline
               lines = []

    if doShuffle:
        random.shuffle(DataList)
        for outline in DataList:
            print outline

run()