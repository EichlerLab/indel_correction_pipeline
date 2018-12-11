#!/usr/bin/env python

import argparse                                                                
import subprocess
import os
import sys
import re
import datetime




nRunFreebayesPlusMinusThisDistanceFromGenekillingIndels = 100
 
nHowFarToLookForCorrelatedIndels = 6

szSmartiePipelineDir = "/net/eichler/vol20/projects/whole_genome_assembly/nobackups/yoruban/fix_pilon/zev_pipe/smartie-sv/pipeline"

szHg38CDS = "/net/eichler/vol2/eee_shared/assemblies/hg38/genes/refGene.CDS_exons.bed"

szHg38LowConfidenceRegions = "/net/eichler/vol18/zevk/great_apes/smartie/data/hg38_low_confidence_regions.bed"

szHg38WgacSuperDup = "/net/eichler/vol2/eee_shared/assemblies/hg38/wgac/genomicSuperDup.sort.bed.gz"

nHowCloseToEnds = 10000

nMinimumSizeContig = 500000


parser = argparse.ArgumentParser()
parser.add_argument("--szInputGenome", help="can be full path or in current directory", required = True )
parser.add_argument("--szInputFreebayesVCFFile", help="must be with respect to szInputGenome", required = True )
parser.add_argument("--szPrefixForSmartie", help="just so it is unique across genomes", required = True )
parser.add_argument( "--nProcessors", help="number of cpus on a node that has 60G available", required = True, type = int )
parser.add_argument( "--szIlluminaPiecesForward", help="should end with pieces_forward   Replacing pieces_forward by pieces_reverse should give the reverse directory", required = True )
parser.add_argument( "--szHighAndLowDepthRegionsBed", help="I guess in pilon coordinates so a problem", required = True )

args = parser.parse_args()


assert os.path.exists( args.szHighAndLowDepthRegionsBed ), "parameter --szHighAndLowDepthRegionsBed " + args.szHighAndLowDepthRegionsBed + " does not exist"
# I will be using this in a subdirectory so I want it as an absolute path

szCommand = "readlink -f --no-newline " + args.szHighAndLowDepthRegionsBed 
print "about to execute: " + szCommand
szHighAndLowDepthRegionsBedFullPath = subprocess.check_output( szCommand, shell = True )
assert os.path.exists( szHighAndLowDepthRegionsBedFullPath ), "full path " + szHighAndLowDepthRegionsBedFullPath + " didn't exist"


szCurrentDirectory = os.getcwd()

szSmartieOnCorrectedGenome1Subdirectory  = "smartie_on_correctedGenome1"

szSmartieIndelFileCorrectedGenome1 = szCurrentDirectory  + "/" + szSmartieOnCorrectedGenome1Subdirectory + "/variants/hg38-" + args.szPrefixForSmartie + ".1.indel.bed"
if ( os.path.isfile( szSmartieIndelFileCorrectedGenome1 ) ):
    sys.exit( szSmartieIndelFileCorrectedGenome1 + " already exists.  Delete it first." )

szSmartieOnCorrectedGenome2Subdirectory  = "smartie_on_correctedGenome2"

szSmartieIndelFileCorrectedGenome2 = szCurrentDirectory + "/" + szSmartieOnCorrectedGenome2Subdirectory + "/variants/hg38-" + args.szPrefixForSmartie + ".2.indel.bed"
if ( os.path.isfile( szSmartieIndelFileCorrectedGenome2 ) ):
    sys.exit( szSmartieIndelFileCorrectedGenome2 + " already exists.  Delete it first." )

####################################################################
#####  filter existing freebayes file to yield correctFreebayes1.fa genome

szFreebayesCorrection1 = "freebayes_correction1"
szCommand = "mkdir " + szFreebayesCorrection1
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

os.chdir( szFreebayesCorrection1 )

szCommand = "gunzip -c " + args.szInputFreebayesVCFFile + " | grep \"^#\" >header.txt"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand = "gunzip -c " + args.szInputFreebayesVCFFile + " | grep -v \"^#\" | awk '{ if ( length( $5 ) != length( $4 ) ) { {if ( ( $2 - x ) > " + str( nHowFarToLookForCorrelatedIndels ) + " ) print $0; } x = $2 } }' >body.txt"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )


szCurrentDirectory = os.getcwd()
# want this full path since we will use it later in a different directory
# for converting coordinates 
szFilteredFreebayes1VCF = szCurrentDirectory + "/" + "filtered_freebayes1.vcf.gz"

szCommand = "module load tabix/0.2.6  && cat header.txt body.txt | bgzip -c >" + szFilteredFreebayes1VCF
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )


szFreebayesCorrectedGenome1 = szCurrentDirectory + "/" + os.path.basename( args.szInputGenome ) + ".freebayes1.fa"

szCommand = "module load zlib/1.2.11 VCFtools/0.1.12b && module load tabix/0.2.6 && tabix " + szFilteredFreebayes1VCF + " && cat " + args.szInputGenome + " | vcf-consensus " + szFilteredFreebayes1VCF + " >" + szFreebayesCorrectedGenome1
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand = "samtools faidx " + szFreebayesCorrectedGenome1
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szFreebayesCorrectedGenome1Fai = szFreebayesCorrectedGenome1 + ".fai"
assert os.path.isfile( szFreebayesCorrectedGenome1Fai ), "file " + szFreebayesCorrectedGenome1Fai + " should exist at this point"

os.chdir( ".." )

########## aligning reads ###############


# in prepartion for running freebayes, align illumina reads to szFreebayesCorrectedGenome1

szAlignIlluminaSubdirectory  = "alignIlluminaToFreebayesCorrectedGenome1"
szCommand = "mkdir " + szAlignIlluminaSubdirectory
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

os.chdir( szAlignIlluminaSubdirectory )

szCommand = "ln -s " + args.szIlluminaPiecesForward
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szIlluminaPiecesReverse = re.sub( r"pieces_forward", "pieces_reverse", args.szIlluminaPiecesForward )

szCommand = "ln -s " + szIlluminaPiecesReverse
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand = "ln -s " + szFreebayesCorrectedGenome1
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand = "ln -s " + szFreebayesCorrectedGenome1Fai
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

# now this is part of the git repository
#szCommand = "cp ~dgordon/pipelines/align_illumina_against_reference/align_illumina_against_reference.snake ."
#print "about to execute: " + szCommand
#subprocess.call( szCommand, shell = True )

szConfigJson = "align_illumina_against_reference_config.json"
with open( szConfigJson, "w" ) as fConfig:
    fConfig.write( "{\n     \"align_reads_to_this_fasta\": " + "\"" + szFreebayesCorrectedGenome1 + "\",\n     \"tmp_dir\": \"/var/tmp/\" \n}" )




szCommand = "mkdir -p log"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand = "source source_this_first.sh && snakemake -s align_illumina_against_reference.snake --drmaa \" -q eichler-short.q -l h_rt=35:00:00 -V -cwd -e ./log -o ./log {params.sge_opts}  -S /bin/bash\"  -w 300 --jobs 100 -p"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )


szCurrentDirectory = os.getcwd()
szIlluminaVsFreebayesCorrected1Bam = szCurrentDirectory + "/illumina_vs_assembly.sorted.bam"
szIlluminaVsFreebayesCorrected1Bai = szIlluminaVsFreebayesCorrected1Bam + ".bai"

assert os.path.exists( szIlluminaVsFreebayesCorrected1Bam ), szIlluminaVsFreebayesCorrected1Bam + "ust exist at this point but doesn't"
assert os.path.exists( szIlluminaVsFreebayesCorrected1Bai ), szIlluminaVsFreebayesCorrected1Bai + " must exist at this point but doesn't"


os.chdir( ".." )


########### end of aligning reads ##############

########### run smartie pipeline on corrected genome1 ##############

szCommand = "mkdir " + szSmartieOnCorrectedGenome1Subdirectory
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

os.chdir( szSmartieOnCorrectedGenome1Subdirectory )

# next run zev's pipeline

szCommand = "cp ~dgordon/smartie/smartie.snake ."
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand = "cp ~dgordon/smartie/smartie-sv/pipeline/config.sh ."
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )


szSmartiePipelineConfigCorrectedGenome1 = "config.json"


with open( szSmartiePipelineConfigCorrectedGenome1, "w" ) as fSmartiePipelineConfig:
    szFixedPart1 = \
"""\
{
        "install" :"/net/eichler/vol20/projects/whole_genome_assembly/nobackups/yoruban/fix_pilon/zev_pipe/smartie-sv",
        "targets" : {
                  "hg38" : "/net/eichler/vol2/eee_shared/assemblies/hg38/indexes/blasr/ucsc.hg38.no_alts.fasta"
                  },
        "queries" : {
"""
    fSmartiePipelineConfig.write( szFixedPart1 )
    szQueryLine = "                 "
    szQueryLine += "\"" 
    szQueryLine += args.szPrefixForSmartie 
    szQueryLine += ".1"
    szQueryLine += "\" : \"" 
    szQueryLine += szFreebayesCorrectedGenome1 
    szQueryLine += "\"\n"
    fSmartiePipelineConfig.write( szQueryLine )
    szFixedPart2 = \
"""           },
        "blasr_mem" : "10G",
        "processors" : """
    fSmartiePipelineConfig.write( szFixedPart2 )
    fSmartiePipelineConfig.write( "\"" + str( args.nProcessors ) + "\"\n}\n" )


szCommand = "mkdir -p log && snakemake -j 22 --drmaa \" -q eichler-short.q -l h_rt=90:00:00 -V  {params.sge_opts} -cwd -e ./log -o ./log -S /bin/bash\" -s smartie.snake --verbose -p"

print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

# get the indel file from "variants"

assert os.path.isfile( szSmartieIndelFileCorrectedGenome1 ), szSmartieIndelFileCorrectedGenome1 + " must exist at this point but doesn't"

os.chdir( ".." )


########### end of running smartie pipeline on corrected genome1 ##############


###########  filtering indels in corrected genome 1 ##################

# prepare to run freebayes on the problem locations of corrected1 genome
# to get the problem locations, filter the smartie output

szRunFreebayesOnCorrection1Genome = "run_freebayes_on_corrected1_genome"

szCommand = "mkdir " + szRunFreebayesOnCorrection1Genome
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

os.chdir( szRunFreebayesOnCorrection1Genome )

szCurrentDirectory = os.getcwd()
print "current directory: " + szCurrentDirectory


# filter the smartie analysis just as in
# /net/eichler/vol20/projects/whole_genome_assembly/nobackups/yoruban/fix_pilon/freebayes_polishing4/analysis
# to get a bed file of the regions to run freebayes on which is done
# in /net/eichler/vol20/projects/whole_genome_assembly/nobackups/yoruban/fix_pilon/freebayes_polishing5
# by run_freebayes_in_these_regions.sh

# finds indels in human coding regions
szCorrected1GenomeCDS = "hg38-Correct1Genome.cds.bed"
szCommand = "bedtools intersect -wa -u -a " + szSmartieIndelFileCorrectedGenome1 + " -b " + szHg38CDS + " >" + szCorrected1GenomeCDS
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )


# find indels (with respect to human) that are gene-killing (not a
# multiple of 3) in human coding regions and other filters.  I don't
# believe there is any guarantee that freebayes will find any Illumina
# reads that would correct such indels.  But we will try to find such
# by running freebayes at these locations.

szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigs = "indel.cds.noseg-nolowc-noends-nosmallcontigs.bed"

szCommand = "perl -lane 'print if $F[4] % 3 != 0' " + szCorrected1GenomeCDS + " | sort -k1,1 -k2,2n | bedtools subtract -A -a - -b " + szHg38LowConfidenceRegions + " | bedtools subtract -A -a - -b " + szHg38WgacSuperDup + " | perl -lane 'print if ($F[7] > " + str( nHowCloseToEnds ) +  " && ($F[9] - $F[7]) > " + str( nHowCloseToEnds ) + " && $F[9] > " + str( nMinimumSizeContig ) + " ) ' | sort -k1,1 -k2,2n > " + szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigs
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace = "indel.cds.noseg-nolowc-noends-nosmallcontigs.falconspace.bed"

# convert to falcon coordinates:
# (can do this easily because the smartie pipeline gave both human and falcon coordinates
# for each indel)
szCommand = "awk 'BEGIN {OFS = \"\t\" } {print $7,$8,$9,$0}' " + szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigs + " | sort -k1,1 -k2,2n >" +  szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

# remove low and high depth regions

# convert high and low depth regions into coordinates of the 
# freebayes_corrected1 genome

szCommand = "cp ~dgordon/pipelines/freebayes_polishing/convertCoordinates2.py ."
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

# change high/low depth from falcon coordinates *before* this pipeline ran to 
# falcon coordinate of correctedGenome1 coordinates.  Does this using the vcf file
# used to make correctedGenome1

szOutputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1 = szCurrentDirectory + "/outputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1.bed"
szCommand = "./convertCoordinates2.py --szInputHighAndLowDepthBedFile " + szHighAndLowDepthRegionsBedFullPath  + " --szFreebayesVCFFile " + szFilteredFreebayes1VCF + " --szOutputHighAndLowDepthBedFileWithRespectToNewGenome " + szOutputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1 + " --szNewGenomeFaiFile " + szFreebayesCorrectedGenome1Fai
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )


szFilteredGeneKillingIndels = "for_david_final.repair.bed"

szCommand = "bedtools subtract -a " + szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace + " -b " + szOutputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1 + " -A  > " + szFilteredGeneKillingIndels
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

# remove neighboring events in which ins/del cancel each other or in which ins/ins or 
# del/del sum to multiple of 3

szCommand = "cp ~dgordon/pipelines/freebayes_polishing/find_nongenekilling_indels4.py ."
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

# input to find_nongenekilling_indels4.py is for_david_final.repair.bed
# output is pairs_of_indels_removed4.bed

szCommand = "./find_nongenekilling_indels4.py"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szRemainingIndelsCorrectedGenome1 = "pairs_of_indels_removed4.bed"

assert os.path.isfile( szRemainingIndelsCorrectedGenome1 ), szRemainingIndelsCorrectedGenome1 + " must exist at this point but doesn't"

print "number of filtered gene-killing indels:"
szCommand = "wc -l " + szRemainingIndelsCorrectedGenome1
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )






szCommand = "ln -s " + szFreebayesCorrectedGenome1Fai
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szRunFreebayesInTheseRegions = "run_freebayes_in_these_regions.bed"

szCommand = "bedtools slop -b " + str( nRunFreebayesPlusMinusThisDistanceFromGenekillingIndels ) + " -g " + szFreebayesCorrectedGenome1Fai + " -i " + szRemainingIndelsCorrectedGenome1 + " | sort -k1,1 -k2,2n | bedtools merge -i stdin >" + szRunFreebayesInTheseRegions
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

assert os.path.isfile( szRunFreebayesInTheseRegions ), szRunFreebayesInTheseRegions + " should exist at this point but doesn't"


###########  end filtering indels in corrected genome 1 ##################

########### run freebayes at pin-point locations in corrected genome 1 ############



szConfig = "config.yaml"
with open( szConfig, "w" ) as fConfig:
    fConfig.write( "---\n\n" )
    fConfig.write( "populations: populations.txt\n" )
    fConfig.write( "references: [correctedGenome1]\n" )
    fConfig.write( "correctedGenome1:\n" )
    fConfig.write( "    fasta: " + szFreebayesCorrectedGenome1 + "\n" )
    fConfig.write( "    bamlist: correctedGenome1.bamlist.txt\n" )
    fConfig.write( "    exclude: []\n" )


with open( "correctedGenome1.bamlist.txt", "w" ) as fBamList:
    fBamList.write( szIlluminaVsFreebayesCorrected1Bam + "\n" )

# create data/regions.txt

szRegionsFile = "data/regions.txt"

szCommand = "mkdir data"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )
  
szCommand =  "awk '{print $1\":\"$2\"-\"$3}' run_freebayes_in_these_regions.bed >" + szRegionsFile
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand =  "cp ~dgordon/pipelines/running_freebayes_on_regions/env.cfg ."
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand =  "cp ~dgordon/pipelines/running_freebayes_on_regions/freebayes_on_regions.snake ."
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand =  "cp ~dgordon/pipelines/running_freebayes_on_regions/run_snakemake.sh ."
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand = "module unload python && run_snakemake.sh"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )




szCurrentDirectory = os.getcwd()
print "current directory: " + szCurrentDirectory

szCommand = "module unload python && module list && run_snakemake.sh"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szFreebayesVCFOnCorrectedGenome1 = szCurrentDirectory + "/final/merged.correctedGenome1.vcf.gz"
assert os.path.isfile( szFreebayesVCFOnCorrectedGenome1 ), szFreebayesVCFOnCorrectedGenome1 + " should exist at this point but doesn't"


os.chdir( ".." )
########### end run freebayes at pin-point locations in corrected genome 1 ############


####################   create correctedGenome2.fa, our final genome ##########################

szCorrectedGenome2Dir = "freebayes_correction2"
szCommand = "mkdir " + szCorrectedGenome2Dir
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

os.chdir( szCorrectedGenome2Dir )

# apply these indels to create correctedGenome2.fa  (In the future, should exclude snps.)

szCurrentDirectory = os.getcwd()
print "current directory: " + szCurrentDirectory

szCorrectedGenome2 = szCurrentDirectory + "/correctedGenome2.fa"

szCommand = "module load VCFtools/0.1.12b && module load tabix/0.2.6 && cat " + szFreebayesCorrectedGenome1 + " | vcf-consensus " + szFreebayesVCFOnCorrectedGenome1 + " >" + szCorrectedGenome2
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

# check that correctedGenome2.fa exists



assert os.path.exists( szCorrectedGenome2 ), szCorrectedGenome2 + " should exist at this point but doesn't"



szCommand = "samtools faidx " + szCorrectedGenome2
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

os.chdir( ".." )

szCommand = "ln -s " + szCorrectedGenome2 + " almost_done.fa"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand = "module load numpy/1.7.0 && module load biopython/1.63 &&  rename_sequences.py --szInputFasta almost_done.fa --szOutputFasta indel_corrected.fa"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand = "samtools faidx indel_corrected.fa"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )



################ done creating correctedGenome2.fa  (In the future, should exclude snps.) ####################


################   run smartie pipeline again to see how many gene-killing indels remain ##########################################################

szCommand = "mkdir " + szSmartieOnCorrectedGenome2Subdirectory
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )


os.chdir( szSmartieOnCorrectedGenome2Subdirectory )
print "chdir to "  + szSmartieOnCorrectedGenome2Subdirectory


szCommand = "cp ~dgordon/smartie/smartie.snake ."
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )


szCommand = "cp ~dgordon/smartie/smartie-sv/pipeline/config.sh ."
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )



szSmartiePipelineConfigCorrectedGenome2 = "config.json"

with open( szSmartiePipelineConfigCorrectedGenome2, "w" ) as fSmartiePipelineConfig:
    szFixedPart1 = \
"""\
{
        "install" :"/net/eichler/vol20/projects/whole_genome_assembly/nobackups/yoruban/fix_pilon/zev_pipe/smartie-sv",
        "targets" : {
                  "hg38" : "/net/eichler/vol2/eee_shared/assemblies/hg38/indexes/blasr/ucsc.hg38.no_alts.fasta"
                  },
        "queries" : {
"""
    fSmartiePipelineConfig.write( szFixedPart1 )
    szQueryLine = "                 "
    szQueryLine += "\"" 
    szQueryLine += args.szPrefixForSmartie 
    szQueryLine += ".2"
    szQueryLine += "\" : \"" 
    szQueryLine += szCorrectedGenome2
    szQueryLine += "\"\n"
    fSmartiePipelineConfig.write( szQueryLine )
    szFixedPart2 = \
"""           },
        "blasr_mem" : "10G",
        "processors" : """
    fSmartiePipelineConfig.write( szFixedPart2 )
    fSmartiePipelineConfig.write( "\"" + str( args.nProcessors ) + "\"\n}\n" )


szCommand = "mkdir -p log && snakemake -j 22 --drmaa \" -q eichler-short.q -l h_rt=90:00:00 -V  {params.sge_opts} -cwd -e ./log -o ./log -S /bin/bash\" -s smartie.snake --verbose -p"

print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

assert os.path.isfile( szSmartieIndelFileCorrectedGenome2 ), szSmartieIndelFileCorrectedGenome2 + " must exist at this point but doesn't"

os.chdir( ".." )
    
# end run smartie pipeline again to see how many gene-killing indels remain
##########################################################################

# now filter the smartie indels and see how many remain

szCorrected2GenomeCDS = "hg38-Correct2Genome.cds.bed"
szCommand = "bedtools intersect -wa -u -a " + szSmartieIndelFileCorrectedGenome2 + " -b " + szHg38CDS + " >" + szCorrected2GenomeCDS
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigs = "indel.cds.noseg-nolowc-noends-nosmallcontigs.bed"

szCommand = "perl -lane 'print if $F[4] % 3 != 0' " + szCorrected2GenomeCDS + " | sort -k1,1 -k2,2n | bedtools subtract -A -a - -b " + szHg38LowConfidenceRegions + " | bedtools subtract -A -a - -b " + szHg38WgacSuperDup + " | perl -lane 'print if ($F[7] > " + str( nHowCloseToEnds ) + " && ($F[9] - $F[7]) > " + str( nHowCloseToEnds ) + " && $F[9] > " + str( nMinimumSizeContig ) + " ) ' | sort -k1,1 -k2,2n > " + szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigs
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace = "indel.cds.noseg-nolowc-noends-nosmallcontigs.falconspace.bed"

# convert to falcon coordinates:
szCommand = "awk 'BEGIN {OFS = \"\t\" } {print $7,$8,$9,$0}' " + szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigs + " | sort -k1,1 -k2,2n >" +  szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

# remove low and high depth regions

szFilteredGeneKillingIndels = "for_david_final.repair.bed"

szCommand = "bedtools subtract -a " + szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace + " -b " + szOutputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1 + " -A  > " + szFilteredGeneKillingIndels
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

# remove neighboring events in which ins/del cancel each other or in which ins/ins or 
# del/del sum to multiple of 3

szCommand = "cp ~dgordon/pipelines/freebayes_polishing/find_nongenekilling_indels4.py ."
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )

szCommand = "./find_nongenekilling_indels4.py"
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )


szRemainingIndelsCorrectedGenome2 = "pairs_of_indels_removed4.bed"

assert os.path.isfile( szRemainingIndelsCorrectedGenome2 ), szRemainingIndelsCorrectedGenome2 + " must exist at this point but doesn't"

print "number of filtered gene-killing indels:"
szCommand = "wc -l " + szRemainingIndelsCorrectedGenome2
print "about to execute: " + szCommand
subprocess.call( szCommand, shell = True )





    
    


                                      
    

