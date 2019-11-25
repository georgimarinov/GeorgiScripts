##################################
#                                #
# Last modified 2017/09/22       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s TSSprofile radius sideposition [-subtractMin smoothradius]' % sys.argv[0]
        print '\tMoothing will be done over windows twice as big as the radius'
        sys.exit(1)

    TSSprofile = sys.argv[1]
    radius = int(sys.argv[2])
    sidePos = int(sys.argv[3])

    doSubMin = False
    if '-subtractMin' in sys.argv:
        doSubMin = True
        smooth = int(sys.argv[sys.argv.index('-subtractMin') + 1])

    CovDict = {}
    AllScores = []
    linelist = open(TSSprofile)
    for line in linelist:
        fields = line.strip().split('\t')
        pos = int(fields[0])
        score = float(fields[1])
        CovDict[pos] = score
        AllScores.append(score)

    if doSubMin:
        AllScoresSmoothened = []
        for i in range(smooth):
            smoothened = sum(AllScores[0:2*smooth])/(2*smooth + 0.0)
            AllScoresSmoothened.append(smoothened)
        for i in range(smooth,len(AllScores) - smooth):
            smoothened = sum(AllScores[i - smooth:i + smooth])/(2*smooth + 0.0)
            AllScoresSmoothened.append(smoothened)
        for i in range(len(AllScores) - smooth,len(AllScores)):
            smoothened = sum(AllScoresSmoothened[len(AllScores) - 2*smooth:len(AllScores)])/(2*smooth + 0.0)
            AllScoresSmoothened.append(smoothened)

#        print len(AllScores), len(AllScoresSmoothened)

        minScore = min(AllScoresSmoothened)

    TSS = 0
    for i in range(-radius,radius):
        if doSubMin:
            TSS += (CovDict[i] - minScore)
        else:
            TSS += CovDict[i]

    Sides = 0
    for i in range(-sidePos-radius,-sidePos + radius):
        if doSubMin:
            Sides += (CovDict[i] - minScore)
        else:
            Sides += CovDict[i]
    for i in range(sidePos-radius,sidePos + radius):
        if doSubMin:
            Sides += (CovDict[i] - minScore)
        else:
            Sides += CovDict[i]
    Sides = Sides/2

    print TSSprofile + '\tTSS to flanks ratio:\t' + str(TSS) + '\t' + str(Sides) + '\t' + str(TSS/Sides)

#    if doSubMin:
#        for i in range(len(AllScores)):
#            outline = str(i-len(AllScores)/2) + '\t' + str(AllScores[i]) + '\t' + str(AllScoresSmoothened[i])
#            print outline

run()
