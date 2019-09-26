#!/usr/bin/env python

#module load numpy/1.7.0
#module load biopython/1.63


import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--szInputFasta", required = True )
parser.add_argument("--szOutputFasta", required = True )
args = parser.parse_args()


from Bio import SeqIO


with open( args.szInputFasta, "r" ) as fInput, open( args.szOutputFasta, "w" ) as fOutput:

    
    for record in SeqIO.parse( fInput, "fasta"):
        szOldName = record.id

        szNewName = re.sub( "_1_", "_", szOldName )
        szNewName = re.sub( "_quiver_pilon", "_qpd", szNewName )
        szNewName = re.sub( "_arrow_pilon",  "_qpd", szNewName )
        szNewName = re.sub( "_qp", "_qpd", szNewName )
    
        record.id = szNewName
        record.description = ''

        SeqIO.write( record, fOutput, "fasta")


