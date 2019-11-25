##################################
#                                #
# Last modified 2018/11/19       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import gzip

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s inputfilename outfilename' % sys.argv[0]
        sys.exit(1)

    inputfilename = sys.argv[1]
    outputfilename = sys.argv[2]

    outfile = open(outputfilename, 'w')

    if inputfilename.endswith('.gz'):
        input_stream = gzip.open(inputfilename)
    else:
        input_stream = open(inputfilename)
    LenDict={}
    seq=0    
    chr=''
    for line in input_stream:
        if line[0]=='>':
            if chr != '' and seq != 0:
                outline = chr+'\t'+str(seq) 
                print outline
                outfile.write(outline+'\n')
                seq=0
            chr=line[1:len(line)].strip()
            if seq=='':
                continue
        else:
            seq=seq+len(line.strip())
    outline = chr+'\t'+str(seq) 
    print outline
    outfile.write(outline+'\n')

    outfile.close()

run()

