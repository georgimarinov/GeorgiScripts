##################################
#                                #
# Last modified 2018/05/10       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import os
import random
import string

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s BAM samtools genome.fa' % sys.argv[0]
        print 'Only run this script for files containing uniquely aligned reads; the script will not check for alignment multiplicity!'
        sys.exit(1)

    BAM = sys.argv[1]
    samtools = sys.argv[2]
    genome = sys.argv[3]

    outfile1 = open(BAM.split('.bam')[0] + '.pseudoRep1.SAM.temp', 'w')
    outfile2 = open(BAM.split('.bam')[0] + '.pseudoRep2.SAM.temp', 'w')

    cmd = samtools + ' view -H ' + BAM
    p = os.popen(cmd, "r")
    line = 'line'
    while line != '':
        line = p.readline().strip()
        if line == '':
            break
        outfile1.write(line.strip() + '\n')
        outfile2.write(line.strip() + '\n')

    cmd = samtools + ' view ' + BAM
    p = os.popen(cmd, "r")
    line = 'line'
    while line != '':
        line = p.readline().strip()
        if line == '':
            break
        fields = line.strip().split('\t')
        if fields[2] == '*':
            continue
        if random.random() >= 0.5:
            outfile1.write(line.strip() + '\n')
        else:
            outfile2.write(line.strip() + '\n')

    outfile1.close()
    outfile2.close()

    cmd = samtools + ' view -bT ' + genome + ' ' + BAM.split('.bam')[0] + '.pseudoRep1.SAM.temp' + ' -o ' + BAM.split('.bam')[0] + '.pseudoRep1.bam'
    os.system(cmd)
    cmd = samtools + ' view -bT ' + genome + ' ' + BAM.split('.bam')[0] + '.pseudoRep2.SAM.temp' + ' -o ' + BAM.split('.bam')[0] + '.pseudoRep2.bam'
    os.system(cmd)

    cmd = 'rm ' + BAM.split('.bam')[0] + '.pseudoRep2.SAM.temp'
    os.system(cmd)
    cmd = 'rm ' + BAM.split('.bam')[0] + '.pseudoRep1.SAM.temp'
    os.system(cmd)

run()

