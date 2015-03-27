##################################
#                                #
# Last modified 10/30/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys

def run():

    if len(sys.argv) < 4:
        print 'usage: python %s genome_fasta readlength outifle' % sys.argv[0]
        print '\tuse - instead of a filename to indicated printing to stdout'
        sys.exit(1)

    genome = sys.argv[1]
    readlength = int(sys.argv[2])
    outputfilename = sys.argv[3]

    doStdOut = False
    if outputfilename == '-':
        doStdOut = True

    SeqDict={}
    lineslist = open(genome)
    currentChr=''
    i=0
    for line in lineslist:
        i+=1
        if i % 1000000 == 0:
            if doStdOut:
                pass
            else:
                print i, 'lines processed'
        if line.startswith('>'):
            if currentChr != '':
                SeqDict[currentChr]=sequence
            currentChr=line.strip().split('>')[1]
            sequence=''
        else:
            sequence+=line.strip()

    SeqDict[currentChr]=sequence

    keys=SeqDict.keys()
    keys.sort()

    if doStdOut:
        pass
    else:
        outfile = open(outputfilename, 'w')

    for chr in keys:
        if doStdOut:
            pass
        else:
            print chr
        for i in range(len(SeqDict[chr])-readlength+1):
            read=SeqDict[chr][i:i+readlength]
            if doStdOut:
                print '>' + chr + ':' + str(i) + '-' + str(i+readlength)
                print read
            else:
                outfile.write('>'+chr+':'+str(i)+'-'+str(i+readlength)+'\n')
                outfile.write(read+'\n')

    if doStdOut:
        pass
    else:
        outfile.close()

run()

