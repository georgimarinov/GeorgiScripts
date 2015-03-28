'''
User specs:
1. Checks if there are NH tags presents
2. If there aren't, it goes through the BAM file, calculates the read multiplicity
3. Then makes a temp BAM file into which all alignments (and unaligned reads) are written but this time with the NH tags
4. Then deletes the original file and replaces with the temp file / renames the temp file to have the name of the original file
5. Also, if we can write the number of reads into the header and check for its presence and if present, use it, if not replace
the file with one that has that information in the same manner, it will save the need to count all the reads every time.
It will have to be several numbers - total reads, unique reads, multi reads, reads mapping on the plus and minus strand.
'''
import sys
import pysam


def main(argv=None):


    if len(sys.argv) < 3:
        print 'usage: python %s BAMfilename outputfilename' % sys.argv[0]
        print '       BAM file has to be indexed'
        sys.exit(1)

    BAM = sys.argv[1]
    outputfilename = sys.argv[2]

    samfile = pysam.Samfile(BAM, "rb" )
    readMultiplicityDict = getReadMultiplicity(samfile)

    counts = getReadCounts(samfile, readMultiplicityDict)
    outheader = buildHeader(samfile.header, counts)
    outfile = pysam.Samfile(outputfilename, "wb", header=outheader)
    for alignedread in samfile.fetch(until_eof=True):
        ID = getReadID(alignedread)
        try:
            scaleby = alignedread.opt('NH')
            outfile.write(alignedread)
        except KeyError:
            scaleby = readMultiplicityDict[ID]
            writeBAMEntry(outfile, alignedread, scaleby)


def getReadMultiplicity(samfile):
    readMultiplicityDict = {}
    processedReads = 0
    for alignedread in samfile.fetch(until_eof=True):
        processedReads += 1
        if processedReads % 5000000 == 0:
                print str(processedReads/1000000) + 'M alignments processed in multiplicity examination'

        ID = getReadID(alignedread)
        try:
            readMultiplicityDict[ID] = alignedread.opt('NH')
        except KeyError:
            try:
                readMultiplicityDict[ID] += 1
            except KeyError:
                readMultiplicityDict[ID] = 1

    return readMultiplicityDict


def getReadCounts(samfile, readMultiplicityDict):
    uniques = 0
    multis = 0
    minusUnique = 0
    plusUnique = 0
    uniqueSplice = 0
    multiSplice = 0
    plusMulti = 0
    minusMulti = 0
    for alignedread in samfile.fetch(until_eof=True):
        ID = getReadID(alignedread)
        if readMultiplicityDict[ID] == 1:
            if alignedread.is_reverse:
                minusUnique += 1
            else:
                plusUnique += 1

            if readIsASplice(alignedread.cigar):
                uniqueSplice += 1
            else:
                uniques += 1
        else:
            if alignedread.is_reverse:
                minusMulti += 1
            else:
                plusMulti += 1

            if readIsASplice(alignedread.cigar):
                multiSplice += 1
            else:
                multis += 1

    totalReads = uniques + multis + uniqueSplice + multiSplice
    counts = ["Total\t%d" % totalReads,
              "Unique\t%d" % uniques,
              "Multis\t%d" % multis,
              "UniqueSplices\t%d" % uniqueSplice,
              "Multisplices\t%d" % multiSplice,
              "PlusUnique\t%d" % plusUnique,
              "PlusMulti\t%d" % plusMulti,
              "MinusUnique\t%d" % minusUnique,
              "MinusMulti\t%d" % minusMulti
    ]

    return counts


def buildHeader(templateheader, counts):
    header = templateheader
    header["CO"] = counts

    return header


def getReadID(alignedread):
    ID = alignedread.qname
    if alignedread.is_read1:
        ID = ID + '/1'

    if alignedread.is_read2:
        ID = ID + '/2'

    return ID


def writeBAMEntry(outfile, originalRead, multiplicity):
    newAlignedRead = originalRead
    newAlignedRead.tags = newAlignedRead.tags + [("NH",multiplicity)]
    outfile.write(newAlignedRead)


def readIsASplice(cigar):
    isSplice = False
    if cigar is not None:
        for (m,bp) in cigar:
            if m == 3:
                isSplice = True
                break

    return isSplice


if __name__ == "__main__":
    main(sys.argv)