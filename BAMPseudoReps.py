##################################
#                                #
# Last modified 03/09/2015       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import os
import random
import pysam
import string

def run():

    if len(sys.argv) < 1:
        print 'usage: python %s BAM' % sys.argv[0]
        print 'Only run this script for files containing uniquely aligned reads; the script will not check for alignment multiplicity!'
        sys.exit(1)

    BAM = sys.argv[1]

    samfile = pysam.Samfile(BAM, "rb" )

    outfile1 = pysam.Samfile(BAM.split('.bam')[0] + '.pseudoRep1.bam', "wb", template=samfile)
    outfile2 = pysam.Samfile(BAM.split('.bam')[0] + '.pseudoRep2.bam', "wb", template=samfile)

    for alignedread in samfile.fetch():
         chr = samfile.getrname(alignedread.tid)
         if chr == '*':
             continue
         if random.random() >= 0.5:
             outfile1.write(alignedread)
         else:
             outfile2.write(alignedread)

    outfile1.close()
    outfile2.close()

run()

