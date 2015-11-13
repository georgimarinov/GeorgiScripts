##################################
#                                #
# Last modified 11/11/2015       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import os

def run():

    if len(sys.argv) < 5:
        print 'usage: python %s config SPP_location samtools BAMPseudoReps.py N_cpus' % sys.argv[0]
        print '\tconfig format:'
        print '\tlabel\tChIP-Rep1.bam\tInput-Rep1.bam\tChIP-Rep2.bam\tInput-Rep2.bam'
        print '\tThe script will print to stdout'
        sys.exit(1)

    config = sys.argv[1]
    SPP = sys.argv[2]
    samtools = sys.argv[3]
    BAMpseudoreps = sys.argv[4]
    P = sys.argv[5]

    CommandDict = {'mkdirRep':[],
                   'mkdirPooled':[],
                   'mkdirPooledPseudoRep':[],
                   'mkdirIndidividualPseudoRep':[],
                   'PeakCallRep12':[],
                   'PeakCallPooled':[],
                   'PeakCallPooledPseudoRep':[],
                   'PeakCallIndividualPseudoRep':[],
                   'IDRRep12':[],
                   'IDRPooledPseudoRep':[],
                   'IDRIndividualPseudoRep':[],
                   'merge':[],
                   'sort':[],
                   'index':[],
                   'BAMPseudoRep':[],
                   'PlotsRep12':[],
                   'PlotsPooledPseudoRep':[],
                   'PlotsIndividualPseudoRep':[]
                    }
    DataList = []

    linelist = open(config)
    for line in linelist:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        label = fields[0]
        ChIPRep1 = fields[1]
        InputRep1 = fields[2]
        ChIPRep2 = fields[3]
        InputRep2 = fields[4]
        DataList.append((label,ChIPRep1,InputRep1,ChIPRep2,InputRep2))

    for (label,ChIPRep1,InputRep1,ChIPRep2,InputRep2) in DataList:
        Rep1Dir = 'SPP-300K-' + label + '-Rep1'
        Rep2Dir = 'SPP-300K-' + label + '-Rep2'
        CommandDict['mkdirRep'].append('mkdir ' + Rep1Dir)
        CommandDict['mkdirRep'].append('mkdir ' + Rep2Dir)

        PooledDir = 'SPP-300K-' + label + '-Pooled'
        CommandDict['mkdirPooled'].append('mkdir ' + PooledDir)

        PooledPseudoRep1Dir = 'SPP-300K-' + label + '-PooledPseudoRep1'
        PooledPseudoRep2Dir = 'SPP-300K-' + label + '-PooledPseudoRep2'
        CommandDict['mkdirPooledPseudoRep'].append('mkdir ' + PooledPseudoRep1Dir)
        CommandDict['mkdirPooledPseudoRep'].append('mkdir ' + PooledPseudoRep2Dir)

        Rep1IndividualPseudoRep1Dir = 'SPP-300K-' + label + '-Rep1PseudoRep1'
        Rep1IndividualPseudoRep2Dir = 'SPP-300K-' + label + '-Rep1PseudoRep2'
        Rep2IndividualPseudoRep1Dir = 'SPP-300K-' + label + '-Rep2PseudoRep1'
        Rep2IndividualPseudoRep2Dir = 'SPP-300K-' + label + '-Rep2PseudoRep2'
        CommandDict['mkdirIndidividualPseudoRep'].append('mkdir ' + Rep1IndividualPseudoRep1Dir)
        CommandDict['mkdirIndidividualPseudoRep'].append('mkdir ' + Rep1IndividualPseudoRep2Dir)
        CommandDict['mkdirIndidividualPseudoRep'].append('mkdir ' + Rep2IndividualPseudoRep1Dir)
        CommandDict['mkdirIndidividualPseudoRep'].append('mkdir ' + Rep2IndividualPseudoRep2Dir)

        CommandDict['PeakCallRep12'].append('Rscript ' + SPP + ' -c=' + ChIPRep1 + ' -i=' + InputRep1 + ' -p=' + P + ' -npeak=300000 -savr -savp -rf -odir=' + Rep1Dir)
        CommandDict['PeakCallRep12'].append('Rscript ' + SPP + ' -c=' + ChIPRep2 + ' -i=' + InputRep2 + ' -p=' + P + ' -npeak=300000 -savr -savp -rf -odir=' + Rep2Dir)

        PooledChIPBAMmerged = label + '.ChIP.pooled.bam'
        PooledInputBAMmerged = label + '.Control.pooled.bam'
        CommandDict['merge'].append(samtools + ' merge ' + PooledChIPBAMmerged + ' ' + ChIPRep1 + ' ' + ChIPRep2)
        CommandDict['merge'].append(samtools + ' merge ' + PooledInputBAMmerged + ' ' + InputRep1 + ' ' + InputRep2)

        PooledChIPBAMsorted = label + '.ChIP.pooled.sorted.bam'
        PooledInputBAMsorted = label + '.Control.pooled.sorted.bam'
        CommandDict['sort'].append(samtools + ' sort ' + PooledChIPBAMmerged + ' ' + PooledChIPBAMsorted[0:-4])
        CommandDict['sort'].append(samtools + ' sort ' + PooledInputBAMmerged + ' ' + PooledInputBAMsorted[0:-4])

        CommandDict['index'].append(samtools + ' index ' + PooledChIPBAMsorted)
        CommandDict['index'].append(samtools + ' index ' + PooledInputBAMsorted)

        CommandDict['BAMPseudoRep'].append('python ' + BAMpseudoreps + ' ' + PooledChIPBAMsorted)
        CommandDict['BAMPseudoRep'].append('python ' + BAMpseudoreps + ' ' + PooledInputBAMsorted)
        CommandDict['BAMPseudoRep'].append('python ' + BAMpseudoreps + ' ' + ChIPRep1)
        CommandDict['BAMPseudoRep'].append('python ' + BAMpseudoreps + ' ' + ChIPRep2)
        CommandDict['BAMPseudoRep'].append('python ' + BAMpseudoreps + ' ' + InputRep1)
        CommandDict['BAMPseudoRep'].append('python ' + BAMpseudoreps + ' ' + InputRep2)

        CommandDict['PeakCallPooled'].append('Rscript ' + SPP + ' -c=' + PooledChIPBAMsorted + ' -i=' + PooledInputBAMsorted + ' -p=' + P + ' -npeak=300000 -savr -savp -rf -odir=' + PooledDir)

        PooledChIPBAMPseudoRep1 = PooledChIPBAMsorted[0:-4] + '.pseudoRep1.bam'
        PooledChIPBAMPseudoRep2 = PooledChIPBAMsorted[0:-4] + '.pseudoRep2.bam'
        PooledInputBAMPseudoRep1 = PooledInputBAMsorted[0:-4] + '.pseudoRep1.bam'
        PooledInputBAMPseudoRep2 = PooledInputBAMsorted[0:-4] + '.pseudoRep2.bam'

        CommandDict['PeakCallPooledPseudoRep'].append('Rscript ' + SPP + ' -c=' + PooledChIPBAMPseudoRep1 + ' -i=' + PooledInputBAMPseudoRep1 + ' -p=' + P + ' -npeak=300000 -savr -savp -rf -odir=' + PooledPseudoRep1Dir)
        CommandDict['PeakCallPooledPseudoRep'].append('Rscript ' + SPP + ' -c=' + PooledChIPBAMPseudoRep2 + ' -i=' + PooledInputBAMPseudoRep2 + ' -p=' + P + ' -npeak=300000 -savr -savp -rf -odir=' + PooledPseudoRep2Dir)

        Rep1ChIPIndividualPseudoRep1 = ChIPRep1[0:-4] + '.pseudoRep1.bam'
        Rep2ChIPIndividualPseudoRep1 = ChIPRep2[0:-4] + '.pseudoRep1.bam'
        Rep1InputIndividualPseudoRep1 = InputRep1[0:-4] + '.pseudoRep1.bam'
        Rep2InputIndividualPseudoRep1 = InputRep2[0:-4] + '.pseudoRep1.bam'
        Rep1ChIPIndividualPseudoRep2 = ChIPRep1[0:-4] + '.pseudoRep2.bam'
        Rep2ChIPIndividualPseudoRep2 = ChIPRep2[0:-4] + '.pseudoRep2.bam'
        Rep1InputIndividualPseudoRep2 = InputRep1[0:-4] + '.pseudoRep2.bam'
        Rep2InputIndividualPseudoRep2 = InputRep2[0:-4] + '.pseudoRep2.bam'
        CommandDict['PeakCallIndividualPseudoRep'].append('Rscript ' + SPP + ' -c=' + Rep1ChIPIndividualPseudoRep1 + ' -i=' + Rep1InputIndividualPseudoRep1 + ' -p=' + P + ' -npeak=300000 -savr -savp -rf -odir=' + Rep1IndividualPseudoRep1Dir)
        CommandDict['PeakCallIndividualPseudoRep'].append('Rscript ' + SPP + ' -c=' + Rep1ChIPIndividualPseudoRep2 + ' -i=' + Rep1InputIndividualPseudoRep2 + ' -p=' + P + ' -npeak=300000 -savr -savp -rf -odir=' + Rep1IndividualPseudoRep2Dir)
        CommandDict['PeakCallIndividualPseudoRep'].append('Rscript ' + SPP + ' -c=' + Rep2ChIPIndividualPseudoRep1 + ' -i=' + Rep2InputIndividualPseudoRep1 + ' -p=' + P + ' -npeak=300000 -savr -savp -rf -odir=' + Rep2IndividualPseudoRep1Dir)
        CommandDict['PeakCallIndividualPseudoRep'].append('Rscript ' + SPP + ' -c=' + Rep2ChIPIndividualPseudoRep2 + ' -i=' + Rep2InputIndividualPseudoRep2 + ' -p=' + P + ' -npeak=300000 -savr -savp -rf -odir=' + Rep2IndividualPseudoRep2Dir)

        CommandDict['IDRRep12'].append('Rscript batch-consistency-analysis.r' + ' ' + Rep1Dir + '/' + ChIPRep1.split('/')[-1][0:-4] + '_VS_' + InputRep1.split('/')[-1][0:-4] + '.regionPeak' + 
                                                                                ' ' + Rep2Dir + '/' + ChIPRep2.split('/')[-1][0:-4] + '_VS_' + InputRep2.split('/')[-1][0:-4] + '.regionPeak' + 
                                                                            ' -1 IDR-SPP-'  + label + ' 0 F signal.value')

        CommandDict['IDRPooledPseudoRep'].append('Rscript batch-consistency-analysis.r' + ' ' + PooledPseudoRep1Dir + '/' + PooledChIPBAMPseudoRep1[0:-4] + '_VS_' + PooledInputBAMPseudoRep1[0:-4] + '.regionPeak' + 
                                                                                          ' ' + PooledPseudoRep2Dir + '/' + PooledChIPBAMPseudoRep2[0:-4] + '_VS_' + PooledInputBAMPseudoRep2[0:-4] + '.regionPeak' + 
                                                                            ' -1 IDR-SPP-'  + label + '-PooledPseudoReps 0 F signal.value')

        CommandDict['IDRIndividualPseudoRep'].append('Rscript batch-consistency-analysis.r' + ' ' + Rep1IndividualPseudoRep1Dir + '/' + Rep1ChIPIndividualPseudoRep1[0:-4] + '_VS_' + Rep1InputIndividualPseudoRep1[0:-4] + '.regionPeak' + 
                                                                                              ' ' + Rep1IndividualPseudoRep2Dir + '/' + Rep1ChIPIndividualPseudoRep2[0:-4] + '_VS_' + Rep1InputIndividualPseudoRep2[0:-4] + '.regionPeak' + 
                                                                            ' -1 IDR-SPP-'  + label + '-Rep1PseudoReps 0 F signal.value')
        CommandDict['IDRIndividualPseudoRep'].append('Rscript batch-consistency-analysis.r' + ' ' + Rep2IndividualPseudoRep1Dir + '/' + Rep2ChIPIndividualPseudoRep1[0:-4] + '_VS_' + Rep2InputIndividualPseudoRep1[0:-4] + '.regionPeak' + 
                                                                                              ' ' + Rep2IndividualPseudoRep2Dir + '/' + Rep2ChIPIndividualPseudoRep2[0:-4] + '_VS_' + Rep2InputIndividualPseudoRep2[0:-4] + '.regionPeak' + 
                                                                            ' -1 IDR-SPP-'  + label + '-Rep2PseudoReps 0 F signal.value')


        CommandDict['PlotsRep12'].append('Rscript batch-consistency-plot.r 1 ' + 'IDR-SPP-'  + label + ' IDR-SPP-'  + label)

        CommandDict['PlotsPooledPseudoRep'].append('Rscript batch-consistency-plot.r 1 ' + 'IDR-SPP-'  + label + '-PooledPseudoReps' + ' IDR-SPP-'  + label + '-PooledPseudoReps')

        CommandDict['PlotsIndividualPseudoRep'].append('Rscript batch-consistency-plot.r 1 ' + 'IDR-SPP-'  + label + '-Rep1PseudoReps' + ' IDR-SPP-'  + label + '-Rep1PseudoReps')
        CommandDict['PlotsIndividualPseudoRep'].append('Rscript batch-consistency-plot.r 1 ' + 'IDR-SPP-'  + label + '-Rep2PseudoReps' + ' IDR-SPP-'  + label + '-Rep2PseudoReps')

    print '# make SPP output folders for individual replicates:'
    for command in CommandDict['mkdirRep']:
        print command
    print '\n# make SPP output folders for pooled runs:'
    for command in CommandDict['mkdirPooled']:
        print command
    print '\n# make SPP output folders for pooled pseudoreplicates:'
    for command in CommandDict['mkdirPooledPseudoRep']:
        print command
    print '\n# make SPP output folders for individual pseudoreplicates:'
    for command in CommandDict['mkdirIndidividualPseudoRep']:
        print command
    print '\n# merge replicate BAM files:'
    for command in CommandDict['merge']:
        print command
    print '\n# sort merged BAM files:'
    for command in CommandDict['sort']:
        print command
    print '\n# index sorted BAM files:'
    for command in CommandDict['index']:
        print command
    print '\n# Generate pseudoreplicates:'
    for command in CommandDict['BAMPseudoRep']:
        print command
    print '\n# Peak calls for individual replicates:'
    for command in CommandDict['PeakCallRep12']:
        print command
    print '\n# Peak calls for for pooled datasets:'
    for command in CommandDict['PeakCallPooled']:
        print command
    print '\n# Peak calls for for pooled pseudoreplicates:'
    for command in CommandDict['PeakCallPooledPseudoRep']:
        print command
    print '\n# Peak calls for for individual pseudoreplicates:'
    for command in CommandDict['PeakCallIndividualPseudoRep']:
        print command
    print '\n# Clean up pseudoreplicate BAM files:'
    print 'rm *.pooled.bam'
    print 'rm *.sorted.bam'
    print 'rm *.pseudoRep1.bam'
    print 'rm *.pseudoRep2.bam'
    print 'rm *.sorted.bamb.bai'
    print 'rm *.pseudoRep1.bam.bai'
    print 'rm *.pseudoRep2.bam.bai'
    print '\n# unzip peak call files:'
    print 'gunzip SPP-300K*/*gz'
    print '\n# IDR for individual replicates:'
    for command in CommandDict['IDRRep12']:
        print command
    print '\n# IDR for pooled pseudoreplicates:'
    for command in CommandDict['IDRPooledPseudoRep']:
        print command
    print '\n# IDR for individual pseudoreplicates:'
    for command in CommandDict['IDRIndividualPseudoRep']:
        print command
    print '\n# IDR plots for individual replicates:'
    for command in CommandDict['PlotsRep12']:
        print command
    print '\n# IDR plots for pooled pseudoreplicates:'
    for command in CommandDict['PlotsPooledPseudoRep']:
        print command
    print '\n# IDR plots individual pseudoreplicates:'
    for command in CommandDict['PlotsIndividualPseudoRep']:
        print command
    print '\n# compress peak call files:'
    print 'gzip SPP-300K*/*eak'

run()
