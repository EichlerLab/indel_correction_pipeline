#!/usr/bin/env python


# differs from find_nongenekilling_indels3.py in that it looks for pairs
# of deletions (or pairs of insertions ) where the sum of the 2 events 
# is a multiple of 3

szInput = "for_david_final.repair.bed"
szPairsOfIndels = "pairs_of_indels4.bed"
szPairsOfIndelsRemoved = "pairs_of_indels_removed4.bed"
szTempFile = "sorted.bed"

nDoubleInsOrDoubleDelDistance = 6
nInsDelPairDistance = 6

nNumberOfNonGeneKillingDoubleInsOrDoubleDelEvents = 0
nNumberOfNonGeneKillingInsDelPairedEvents = 0
nNumberOfGeneKillingEventsRemaining = 0

import subprocess

szCommand = "sort -k1,1 -k2,2n  " + szInput + " >sorted.bed"
subprocess.call( szCommand, shell = True )


dictLinesWeDoNotWant = {}
with open( szTempFile, "r" ) as fTempFile, open( szPairsOfIndels, "w" ) as fPairsOfIndels:
    aInputLines = fTempFile.readlines()

    n = 0
    while( (n + 1) < len( aInputLines ) ):
        szLineA = aInputLines[n]
        szLineB = aInputLines[n+1]
        
        aWordsA = szLineA.split()
        aWordsB = szLineB.split()

        # lookes like:
        # 000000F_1_27017830_quiver_pilon 444754  444758  chr13   113325660       113325660       insertion       4       -      
        # 0                                  1      2      3          4               5                6          7

        nEndA = int( aWordsA[2] )
        nStartB = int( aWordsB[1] )

        if ( ( aWordsA[0] == aWordsB[0] ) and ( aWordsA[6] != aWordsB[6] ) and ( ( nStartB - nEndA ) < nInsDelPairDistance ) and ( aWordsA[7] == aWordsB[7] ) ):
            dictLinesWeDoNotWant[n] = 1
            dictLinesWeDoNotWant[n+1] = 1
            fPairsOfIndels.write( szLineA )
            fPairsOfIndels.write( szLineB )
            nNumberOfNonGeneKillingInsDelPairedEvents += 2
            n += 2
        else:
            nSizeOfEventA = int( aWordsA[7] )
            nSizeOfEventB = int( aWordsB[7] )
            
            if  ( ( aWordsA[0] == aWordsB[0] ) and ( aWordsA[6] == aWordsB[6] ) and ( nStartB - nEndA < nDoubleInsOrDoubleDelDistance ) and ( ( nSizeOfEventA + nSizeOfEventB ) % 3 == 0 ) ):
                nNumberOfNonGeneKillingDoubleInsOrDoubleDelEvents += 2
                dictLinesWeDoNotWant[n] = 1
                dictLinesWeDoNotWant[n+1] = 1
                fPairsOfIndels.write( szLineA )
                fPairsOfIndels.write( szLineB )
                n += 2
            else:
                n += 1


with open( szInput, "r" ) as fInput:
    aLines = fInput.readlines()

with open( szPairsOfIndelsRemoved, "w" ) as fPairsOfIndelsRemoved:
    for n in range( 0, len( aLines) ):
        if ( n not in dictLinesWeDoNotWant ):
            fPairsOfIndelsRemoved.write( aLines[n] )
            nNumberOfGeneKillingEventsRemaining += 1


print "nNumberOfNonGeneKillingDoubleInsOrDoubleDelEvents = {:d}".format( nNumberOfNonGeneKillingDoubleInsOrDoubleDelEvents )
print "nNumberOfNonGeneKillingInsDelPairedEvents = {:d}".format( nNumberOfNonGeneKillingInsDelPairedEvents )
print "nNumberOfGeneKillingEventsRemaining = {:d}".format( nNumberOfGeneKillingEventsRemaining )

with open( "summary.txt", "w" ) as fSummary:
    fSummary.write( "nNumberOfNonGeneKillingDoubleInsOrDoubleDelEvents = {:d}\n".format( nNumberOfNonGeneKillingDoubleInsOrDoubleDelEvents ) )
    fSummary.write( "nNumberOfNonGeneKillingInsDelPairedEvents = {:d}\n".format( nNumberOfNonGeneKillingInsDelPairedEvents ) )
    fSummary.write( "nNumberOfGeneKillingEventsRemaining = {:d}\n".format( nNumberOfGeneKillingEventsRemaining ) )

    
