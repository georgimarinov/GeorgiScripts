##################################
#                                #
# Last modified 11/04/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string

def run():

    if len(sys.argv) < 4:
        print 'usage: python %s wanted fieldID fasta outputfilename' % sys.argv[0]
        sys.exit(1)

    wanted = sys.argv[1]
    fieldID = int(sys.argv[2])
    fasta = sys.argv[3]
    outfilename = sys.argv[4]

    WantedDict = {}
    linelist = open(wanted)
    for line in linelist:
        fields = line.strip().split('\t')
        WantedDict[fields[fieldID]] = 1

    outfile = open(outfilename, 'w')
    
    inputdatafile = open(fasta)
    Keep = False
    for line in inputdatafile:
        if line[0]=='>':
            ID = line.strip().split('>')[1]
            if WantedDict.has_key(ID):
                Keep = True
            else:
                Keep = False
        else:
            pass
        if Keep:
            outfile.write(line)   

    outfile.close()
   
run()
