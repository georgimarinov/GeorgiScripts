##################################
#                                #
# Last modified 2021/04/08       # 
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
        print 'usage: python %s list_of_UMIs_per_cell_files outfile' % sys.argv[0]
        sys.exit(1)

    BCDict = {}

    lineslist = open(sys.argv[1])
    for L in lineslist:
        file = L.strip().split('\t')[0]
        lines = open(file)
        print file
        for line in lines:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            BC = fields[0]
            UMI = int(fields[2])
            AP = int(fields[3])
            G = int(fields[4])
            if BCDict.has_key(BC):
                pass
            else:
                BCDict[BC] = {}
                BCDict[BC]['UMIs'] = 0
                BCDict[BC]['AP'] = 0
                BCDict[BC]['genes'] = 0
            BCDict[BC]['UMIs'] += UMI
            BCDict[BC]['AP'] += AP
            BCDict[BC]['genes'] += G

    UMIList = []
    for BC in BCDict.keys():
        UMIList.append((BCDict[BC]['UMIs'],BC))

    UMIList.sort()
    UMIList.reverse()

    outfile = open(sys.argv[2],'w')

    outline = '#BC1+BC2+BC3\trank3\tUMIs3\tAligned_Positions3\tgenes'
    outfile.write(outline + '\n')

    R = 0
    for (UMIs,BC) in UMIList:
        R += 1
        outline = BC + '\t' + str(R) + '\t' + str(BCDict[BC]['UMIs']) + '\t' + str(BCDict[BC]['AP']) + '\t' + str(BCDict[BC]['genes'])
        outfile.write(outline + '\n')
            
run()
