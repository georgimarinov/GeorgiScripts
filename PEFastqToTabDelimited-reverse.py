##################################
#                                #
# Last modified 2018/12/30       # 
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
        print 'usage: python %s input prefix [-slashend]' % sys.argv[0]
        print '\tThe script can read compressed files as long as they have the correct suffix - .bz2 or .gz'
        print '\tThe script can also read stdin (-)'
        sys.exit(1)

    input = sys.argv[1]
    outprefix = sys.argv[2]

    doSlashEnd = False
    if '-slashend' in sys.argv:
        doSlashEnd = True

    outfile1 = open(outprefix + '.end1.fastq', 'w')
    outfile2 = open(outprefix + '.end2.fastq', 'w')

    doStdIn = False
    if input == '-':
        doStdIn = True
    elif input.endswith('.bz2'):
        cmd1 = 'bzip2 -cd ' + input
        p1 = os.popen(cmd1, "r")
    elif input.endswith('.gz'):
        cmd1 = 'gunzip -c ' + input
        p1 = os.popen(cmd1, "r")
    else:
        cmd1 = 'cat ' + input
        p1 = os.popen(cmd1, "r")
    line = 'line'
    lines = []
    while line != '':
        if doStdIn:
            line = sys.stdin.readline()
        else:
            line = p1.readline()
        if line == '':
            break
        fields = line.strip().split('\t')
        ID = fields[0]
        seq1 = fields[1]
        q1 = fields[2]
        seq2 = fields[3]
        q2 = fields[4]
        if doSlashEnd:
            ID1 = ID + '/1'
            ID2 = ID + '/2'
        else:
            if '_1:N:0:' in ID:
                ID1 = ID.replace('_1:N:0:',' 1:N:0:')
                ID2 = ID.replace('_1:N:0:',' 2:N:0:')
            elif ' 1:N:0:' in ID:
                ID1 = ID
                ID2 = ID.replace(' 1:N:0:',' 2:N:0:')
            else:
                ID1 = ID + ' 1:N:0:'
                ID2 = ID + ' 2:N:0:'
        outfile1.write('@' + ID1 + '\n')
        outfile2.write('@' + ID2 + '\n')
        outfile1.write(seq1 + '\n')
        outfile2.write(seq2 + '\n')
        outfile1.write('+' + '\n')
        outfile2.write('+' + '\n')
        outfile1.write(q1 + '\n')
        outfile2.write(q2 + '\n')

    outfile1.close()
    outfile2.close()

run()