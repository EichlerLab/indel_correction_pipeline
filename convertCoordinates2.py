#!/usr/bin/env python

import argparse                                                                
import subprocess
import os
import sys
import re
import datetime
import gzip


parser = argparse.ArgumentParser()
parser.add_argument("--szInputHighAndLowDepthBedFile", help="can be full path or in current directory", required = True )
parser.add_argument("--szFreebayesVCFFile", help="what was used to convert the old genome file to the new genome file vis vcf-consensus", required = True )
parser.add_argument("--szOutputHighAndLowDepthBedFileWithRespectToNewGenome", required = True )
parser.add_argument( "--szNewGenomeFaiFile", help="just to check that the coordinate change is within range of the new genome", required = True )

args = parser.parse_args()


aContigLengths = {}


aConversionLines = []


def nCompare( szContigA, nPosA, szContigB, nPosB ):
    if ( szContigA < szContigB ):
        return -1
    elif( szContigA > szContigB ):
        return 1
    else:
        if ( nPosA < nPosB ):
            return -1
        elif( nPosA > nPosB ):
            return 1
        else:
            return 0


def nFindMatchOrPredecessor( szContig, nPos ):

    aTuple = aConversionLines[0]

    if ( nCompare( szContig, nPos, aTuple[0], aTuple[1] ) == -1 ):
        # item is below bottom of list
        return -1;

    nMinIndex = 0

    nTopIndex = len( aConversionLines ) - 1
    aTuple = aConversionLines[ nTopIndex ]
    if ( nCompare( aTuple[0], aTuple[1],szContig, nPos ) <= 0 ):
        # item is above top of list
        return nTopIndex

    nTooBigIndex = nTopIndex

    # if reached here, soMatch < operator[]( nTooBigIndex )
    # and operator[](0) <= soMatch
    # Thus nTooBigIndex != 0  and the correct index is somewhere
    # less than nTooBigIndex and correct index >= 0 so 
    # correct index >= 0.
    # Thus bisect this range over and over, making it smaller
    # and smaller until it is 0.

    while( True ):
        if ( nTooBigIndex - nMinIndex <= 1 ):
            return nMinIndex
        else:
           nTestIndex = ( nTooBigIndex + nMinIndex ) / 2;
           aTuple = aConversionLines[ nTestIndex ]
           if ( nCompare( szContig, nPos, aTuple[0], aTuple[1] ) == -1 ):
              nTooBigIndex = nTestIndex;
           else:
              nMinIndex = nTestIndex;

      
    





def nFindNewPosition( szContig, nPos ):
    "converts old position to new position"

    # bDebug = False
    # if ( ( "000001F_1_22489689_quiver_pilon" == szContig) and (nPos == 1)):
    #     bDebug = True


    nIndex = nFindMatchOrPredecessor( szContig, nPos )

    if ( nIndex == -1 ):
        # either there is no entries for the contig
        # or this is before the first entry for this contig
        nCorrection = 0
    else:
        aTuple = aConversionLines[nIndex]
        if ( aTuple[0] != szContig ):
            # this is before the first entry for this contig
            nCorrection = 0
        else:
            # correction of position <= nPos
            nCorrection = aTuple[2]


    # last found nCorrection

    nNewPosition = nPos + nCorrection

    assert szContig in aContigLengths, szContig + " not found in " + args.szNewGenomeFaiFile + " perhaps the vcf is for a different genome?"
    nLargestPositionThisContig = aContigLengths[ szContig ]

    assert 1 <= nNewPosition, " coordinates should start at 1 but is " + str( nNewPosition ) + " converted from " + szContig + " " + str( nPos )
    
    assert nNewPosition <= nLargestPositionThisContig, szContig + " is length " + str( nLargestPositionThisContig ) + " but " + str( nPos ) + " converted to " + str( nNewPosition )

#    print szContig, nPos, nPos + nCorrection
    return nNewPosition


############################  main ###########################

with open( args.szNewGenomeFaiFile, "r" ) as fFai:
    for szLine in fFai.readlines():
        aWords = szLine.split()
        szContig = aWords[0]
        nLength = int( aWords[1] )
        aContigLengths[ szContig ] = nLength


with gzip.open( args.szFreebayesVCFFile, "r" ) as fVCF:
    for szLine in fVCF:
        if ( szLine.startswith( "#" ) ):
            continue
        aWords = szLine.split()
        # looks like:
        # 000000F_1_80365438_quiver_pilon 20442   .       GAAAAAAAAAAAAAC GAAAAAAAAAAAAC  20102.1 . 
        # or
        # 000000F_1_80365438_quiver_pilon 63212   .       ATG     ATCG,ATTG       28726.2 .   
        # 0                                       1        2       3
        szContig = aWords[0]
        nContigPos = int( aWords[1] )

        szOriginal = aWords[3]
        szNew = aWords[4]

        aNewWords = szNew.split(",")
        szNewFirst = aNewWords[0]

        nDiff = len( szNewFirst ) - len( szOriginal )
        if ( nDiff == 0 ):
            continue

        aTriple = [ szContig, nContigPos,  nDiff ]
        aConversionLines.append( aTriple )

sys.stderr.write( "done reading/parsing vcf file\n" )


# now let's sort and accumulate the differences



# sort by secondary key first
aConversionLines = sorted( aConversionLines, key=lambda student: student[1]) 
aConversionLines = sorted( aConversionLines, key=lambda student: student[0])

sys.stderr.write( "done sorting\n" )


szPreviousContig = ""
nPreviousCorrection = -666
for n in range( 0, len( aConversionLines ) ):
    aTuple = aConversionLines[n]
    if ( szPreviousContig == aTuple[0] ):
        nCorrection = nPreviousCorrection + aTuple[2]
        aConversionLines[n] = [ aTuple[0], aTuple[1], nCorrection ]
        nPreviousCorrection = nCorrection
    else:
        szPreviousContig = aTuple[0]
        nPreviousCorrection = aTuple[2]


sys.stderr.write( "completed setup of conversion table\n" )


nLinesRead = 0
with open( args.szInputHighAndLowDepthBedFile, "r" ) as fInputBed, open( args.szOutputHighAndLowDepthBedFileWithRespectToNewGenome, "w" ) as fOutputBed:
    aLines = fInputBed.readlines()
    for szLine in aLines:
        nLinesRead += 1
        if ( ( nLinesRead % 10000 == 0 ) or ( nLinesRead == ( len( aLines) - 1 ) ) ):
            sys.stderr.write( "{:,} lines read out of {:,}\n".format( nLinesRead , len( aLines  ) ) )
        aWords = szLine.split()
        szContig = aWords[0]
        nStart = int( aWords[1] )
        nEnd = int( aWords[2] )

        nNewStart = nFindNewPosition( szContig, nStart + 1 ) - 1
        nNewEnd = nFindNewPosition( szContig, nEnd )

        aWords[1] = str( nNewStart )
        aWords[2] = str( nNewEnd )
        
        szNewLine = "\t".join( aWords )
        
        fOutputBed.write( szNewLine + "\n" )

        
