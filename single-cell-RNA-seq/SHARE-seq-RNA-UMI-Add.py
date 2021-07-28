##################################
#                                #
# Last modified 2020/11/24       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import os
from sets import Set
import Levenshtein
import numpy as np

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s UMIlen read1|read2' % sys.argv[0]
        print '\t stream the output of PEFastqToTabDelimited.py, then capture the output of this script with PEFastqToTabDelimited-reverse.py'
        sys.exit(1)

    lenUMI = int(sys.argv[1])
    readUMI = sys.argv[2]

    lineslist = sys.stdin
    for line in lineslist:
        fields = line.strip().split('\t')
        if readUMI == 'read1':
            UMI = fields[1][0:lenUMI]
        if readUMI == 'read2':
            UMI = fields[3][0:lenUMI]
        print fields[0].split(' ')[0][0:-1] + '+' + UMI + ']' + '\t' + fields[1] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4]
            
run()
