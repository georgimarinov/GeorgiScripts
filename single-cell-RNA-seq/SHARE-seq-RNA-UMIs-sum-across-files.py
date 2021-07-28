##################################
#                                #
# Last modified 2021/04/08       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string

def run():

    if len(sys.argv) < 4:
        print 'usage: python %s list_of_UMIs_tables(output of SHARE-seq_RNA_counts.py) wantedBCs fieldID outfile [-sparse]' % sys.argv[0]
        sys.exit(1)

    filelist = sys.argv[1]
    wantedBC = sys.argv[2]
    wantedBCID = int(sys.argv[3])
    outfilename = sys.argv[4]

    GeneDict = {}
    TotalUMIs = 0

    doSparse = False
    if '-sparse' in sys.argv:
        doSparse = True
        print 'will output a sparse matrix'

    BCDict = {}
    lineslist = open(wantedBC)
    for line in lineslist:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        BC = fields[wantedBCID]
        BCDict[BC] = 1

    lineslist = open(filelist)
    for L in lineslist:
        file = L.strip().split('\t')[0]
        print file
        lines = open(file)
        for line in lines:
            fields = line.strip().split('\t')
            if line.startswith('#'):
                BCPosDict = {}
                for BCID in range(2,len(fields)):
                    BCPosDict[BCID] = fields[BCID]
                continue
            gene = (fields[0],fields[1])
            for i in range(2,len(fields)):
                UMI = fields[i]
                BC = BCPosDict[i]
                if BCDict.has_key(BC):
                    pass
                else:
                    continue
                if GeneDict.has_key(gene):
                    pass
                else:
                    GeneDict[gene] = {}
                if float(UMI) > 0:
                    GeneDict[gene][BC] = UMI
                    TotalUMIs += float(UMI)

    BCs = BCDict.keys()
    BCs.sort()
    BCs.reverse()

    genes = GeneDict.keys()
    genes.sort()
    genes.reverse()

    if doSparse:

        outfile = open(outfilename + '.barcodes.tsv','w')

        for BC in BCs:
            outfile.write(BC + '\n')

        outfile.close()

        outfile = open(outfilename + '.genes.tsv','w')

        for gene in genes:
            (geneID,geneName) = gene
            outline = geneID + '\t' + geneName
            outfile.write(outline + '\n')

        outfile.close()

        outfile = open(outfilename + '.matrix.mtx','w')
        outfile.write('%%MatrixMarket matrix coordinate real general\n')
        outfile.write('%\n')
        outline = str(len(genes)) + ' ' + str(len(BCs)) + ' ' + str(int(TotalUMIs))
        outfile.write(outline + '\n')
        for gene in genes:
            for BC in BCs:
               if GeneDict[gene].has_key(BC):
                   outline = str(genes.index(gene) + 1) + ' ' + str(BCs.index(BC) + 1) + ' ' + GeneDict[gene][BC]
                   outfile.write(outline + '\n')
        outfile.close()
    else:
        outfile = open(outfilename,'w')

        outline = '#geneID\tgeneName'
        for BC in BCs:
            outline = outline + '\t' + BC
        outfile.write(outline + '\n')

        for gene in genes:
            (geneID,geneName) = gene
            outline = geneID + '\t' + geneName
            for BC in BCs:
               if GeneDict[gene].has_key(BC):
                   outline = outline + '\t' + GeneDict[gene][BC]
               else:
                   outline = outline + '\t' + '0'
            outfile.write(outline + '\n')

        outfile.close()
            
run()
