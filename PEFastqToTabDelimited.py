##################################
#                                #
# Last modified 05/26/2014       # 
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
        print 'usage: python %s end1 end2 [-renameReads prefix] [-shuffle]' % sys.argv[0]
        print '\tThe script can read compressed files as long as they have the correct suffix - .bz2 or .gz'
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

    i=0
    for f in range(len(fastq1files)):
        fastq1 = fastq1files[f]
        fastq2 = fastq2files[f]
        if fastq1.endswith('.bz2'):
            cmd1 = 'bzip2 -cd ' + fastq1
        elif fastq1.endswith('.gz'):
            cmd1 = 'gunzip -c ' + fastq1
        else:
            cmd1 = 'cat ' + fastq1
        p1 = os.popen(cmd1, "r")

        if fastq2.endswith('.bz2'):
            cmd2 = 'bzip2 -cd ' + fastq2
        elif fastq2.endswith('.gz'):
            cmd2 = 'gunzip -c ' + fastq2
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
               line1 = p2.readline().strip()
               line2 = p2.readline().strip()
               line3 = p2.readline().strip()
               line4 = p2.readline().strip()
               if doRename:
                   outline = prefix + str(i/4) + '\t' + lines[1] + '\t' + lines[3] + '\t' + line2 + '\t' + line4
               else:
                   outline = lines[0][1:len(lines[0])].split('/')[0].split(' ')[0] + '\t' + lines[1] + '\t' + lines[3] + '\t' + line2 + '\t' + line4
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