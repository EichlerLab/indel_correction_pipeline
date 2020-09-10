#!/usr/bin/env python

import argparse                                                                
import subprocess
import os
import sys
import re
import datetime




nRunFreebayesPlusMinusThisDistanceFromGenekillingIndels = 100
 
nHowFarToLookForCorrelatedIndels = 6

szHg38CDS = "/net/eichler/vol26/eee_shared/assemblies/hg38/legacy/genes/refGene.CDS_exons.bed"

szHg38LowConfidenceRegions = "~dgordon/pipelines/indel_correction_pipeline/hg38_low_confidence_regions.bed"

szHg38WgacSuperDup = "/net/eichler/vol26/eee_shared/assemblies/hg38/legacy/wgac/genomicSuperDup.sort.bed.gz"

nHowCloseToEnds = 10000

nMinimumSizeContig = 500000


parser = argparse.ArgumentParser()
parser.add_argument("--szInputGenome", help="must be full path for read_depth pipeline", required = True )
parser.add_argument("--szPrefixForSmartie", help="just so it is unique across genomes", required = True )
parser.add_argument( "--nProcessors", help="number of cpus on a node that has 60G available", required = True, type = int )
parser.add_argument( "--szIlluminaPiecesForward", help="should end with pieces_forward   Replacing pieces_forward by pieces_reverse should give the reverse directory", required = True )
parser.add_argument( "--szPacBioReadsFof", help="to calculate read depth profile and find too high and too low read depth", required = True )

args = parser.parse_args()


# assert os.path.exists( args.szHighAndLowDepthRegionsBed ), "parameter --szHighAndLowDepthRegionsBed " + args.szHighAndLowDepthRegionsBed + " does not exist"
# # I will be using this in a subdirectory so I want it as an absolute path

# szCommand = "readlink -f --no-newline " + args.szHighAndLowDepthRegionsBed 
# print "about to execute: " + szCommand
# szHighAndLowDepthRegionsBedFullPath = subprocess.check_output( szCommand, shell = True )
# assert os.path.exists( szHighAndLowDepthRegionsBedFullPath ), "full path " + szHighAndLowDepthRegionsBedFullPath + " didn't exist"





szCurrentDirectory = os.getcwd()
szTopLevelDirectory = szCurrentDirectory

szSmartieOnCorrectedGenome1Subdirectory  = "smartie_on_correctedGenome1"

szSmartieIndelFileCorrectedGenome1 = szCurrentDirectory  + "/" + szSmartieOnCorrectedGenome1Subdirectory + "/variants/hg38-" + args.szPrefixForSmartie + ".1.indel.bed"
#if ( os.path.isfile( szSmartieIndelFileCorrectedGenome1 ) ):
#    sys.exit( szSmartieIndelFileCorrectedGenome1 + " already exists.  Delete it first." )

szSmartieOnCorrectedGenome2Subdirectory  = "smartie_on_correctedGenome2"

szSmartieIndelFileCorrectedGenome2 = szCurrentDirectory + "/" + szSmartieOnCorrectedGenome2Subdirectory + "/variants/hg38-" + args.szPrefixForSmartie + ".2.indel.bed"
#if ( os.path.isfile( szSmartieIndelFileCorrectedGenome2 ) ):
#    sys.exit( szSmartieIndelFileCorrectedGenome2 + " already exists.  Delete it first." )



####################################################################
#####  align pacbio reads against input genome
print "\nAligning pacbio reads against input genome\n"

dirPacBioReadDepthByPosition = "read_depth_by_position"
szCommand = "mkdir -p " + dirPacBioReadDepthByPosition
print "about to execute: " + szCommand
subprocess.check_call( szCommand, shell = True )


szPacBioReadDepthByPositionDone = dirPacBioReadDepthByPosition + "/read_depth_by_position_done"
if ( not os.path.exists( szPacBioReadDepthByPositionDone )):

    os.chdir( dirPacBioReadDepthByPosition )

    szCommand = "git clone git@github.com:EichlerLab/read_depth_by_position.git"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    os.chdir( "read_depth_by_position" )

    with open( "config.json", "w" ) as fConfig:
        fConfig.write( "{\n" )
        fConfig.write( "    \"fasta_to_align_reads_to\": \"" + args.szInputGenome + "\",\n" )
        fConfig.write( "    \"mmi_index_to_align_reads_to\" : \"" + args.szInputGenome + ".mmi\",\n" )
        fConfig.write( "    \"fasta_fofn\": \"" + args.szPacBioReadsFof + "\",\n" )
        fConfig.write( "    \"cluster_settings\" : {\n" )
        fConfig.write( "        \"heavy\" : \"-l h_rt=48:00:00 -l mfree=20G -q eichler-short.q -l disk_free=20G\"\n" )
        fConfig.write( "         }\n" )
        fConfig.write( "}\n" )



    szCommand = "run_calculate_read_depth.sh"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    os.chdir( szTopLevelDirectory )

    szCommand = "touch " + szPacBioReadDepthByPositionDone
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )


szHighAndLowDepthRegionsFalconCoordsBedFullPath = szTopLevelDirectory + "/" + dirPacBioReadDepthByPosition + "/read_depth_by_position/high_and_low_depth_regions.bed"
assert os.path.exists( szHighAndLowDepthRegionsFalconCoordsBedFullPath ), "looking for " + szHighAndLowDepthRegionsFalconCoordsBedFullPath
    




dirBwaVsInputGenome = szTopLevelDirectory + "/bwa_vs_input_genome"
szCommand = "mkdir -p " + dirBwaVsInputGenome
print "about to execute: " + szCommand
subprocess.check_call( szCommand, shell = True )

####################################################################
#####  bwa align Illumina reads against input genome
print "\nbwa align Illumina reads against input genome\n"


# must match name as specified in align_illumina_against_reference pipeline
szBamOfBwaVsInputGenome = dirBwaVsInputGenome + "/illumina_vs_assembly.sorted.bam"

szBwaVsInputGenomeDoneFlag = dirBwaVsInputGenome + "/bwa_vs_input_genome_done"
if ( not os.path.exists( szBwaVsInputGenomeDoneFlag )):

    os.chdir( dirBwaVsInputGenome )

    # cloning into the current directory
    szCommand = "module list && rm -rf log && rm -rf .git && rm -f * && git clone git@github.com:dgordon562/align_illumina_against_reference.git ."
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    
    # set up links to pieces_forward and pieces_reverse of Illumina reads

    szCommand = "ln -s -f " + args.szIlluminaPiecesForward
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szIlluminaPiecesReverse = re.sub( r"pieces_forward", "pieces_reverse", args.szIlluminaPiecesForward )

    szCommand = "ln -s -f " + szIlluminaPiecesReverse
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szConfigJson = "align_illumina_against_reference_config.json"
    with open( szConfigJson, "w" ) as fConfig:
        fConfig.write( "{\n     \"align_reads_to_this_fasta\": " + "\"" + args.szInputGenome + "\",\n     \"tmp_dir\": \"/var/tmp/\" \n}" )


    szCommand = "mkdir -p log"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "source source_this_first.sh && snakemake -s align_illumina_against_reference.snake --drmaa \" -w n -q eichler-short.q -l h_rt=4:00:00:00 -V -cwd -e ./log -o ./log {params.sge_opts} -S /bin/bash\" -w 300 --jobs 100 -p -k --rerun"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szBamBaiOfBwaVsInputGenome =  szBamOfBwaVsInputGenome + ".bai"
    assert os.path.exists( szBamOfBwaVsInputGenome ), szBamOfBwaVsInputGenome + " must exist at this point but doesn't"
    assert os.path.exists( szBamBaiOfBwaVsInputGenome ), szBamBaiOfBwaVsInputGenome + " must exist at this point but doesn't"

    szCommand = "touch " + szBwaVsInputGenomeDoneFlag
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    os.chdir( szTopLevelDirectory )

####################################################################
#####  run freebayes using bwa alignments of illumina reads against input genome

print "\nrun freebayes using bwa alignments of illumina reads against input genome\n"



dirRunFreebayesOnInputGenome = "run_freebayes_on_input_genome"

szCommand = "mkdir -p " + dirRunFreebayesOnInputGenome
print "about to execute: " + szCommand
subprocess.check_call( szCommand, shell = True )
szRunFreebayesOnInputGenomeDoneFlag = dirRunFreebayesOnInputGenome + "/run_freebayes_on_input_genome_done"

# debug
print "checking for: " + szRunFreebayesOnInputGenomeDoneFlag  
# end debug

if ( not os.path.exists( szRunFreebayesOnInputGenomeDoneFlag ) ):
    os.chdir( dirRunFreebayesOnInputGenome ) 

    szCommand = "cp ~/pipelines/freebayes_on_regions/{freebayes_on_regions1.snake,freebayes_on_regions2.snake,env.cfg,run_snakemake.sh} ."
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    with open( "bamlist.fofn", "w" ) as fBamList:
        fBamList.write( szBamOfBwaVsInputGenome + "\n" )


    with open( "config.yaml", "w" ) as fConfigYaml:
        fConfigYaml.write( "---\n" )
        fConfigYaml.write( "\n" )
        fConfigYaml.write( "fasta: " + args.szInputGenome + "\n" )
        fConfigYaml.write( "bamlist: bamlist.fofn\n" )
        fConfigYaml.write( "exclude: []\n" )

    szCommand = "run_snakemake.sh"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "touch run_freebayes_on_input_genome_done"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    os.chdir( szTopLevelDirectory )

szFreebayesOnInputGenomeFullPath = szTopLevelDirectory + "/" + dirRunFreebayesOnInputGenome + "/final/merged.vcf.gz"

assert( os.path.exists( szFreebayesOnInputGenomeFullPath ) )


####################################################################
#####  filter existing freebayes file to yield correctFreebayes1.fa genome

print "\nrun freebayes using bwa alignments of illumina reads against input genome\n"

szFreebayesCorrection1 = "freebayes_correction1"
szCommand = "mkdir -p " + szFreebayesCorrection1
print "about to execute: " + szCommand
subprocess.check_call( szCommand, shell = True )


szFreebayesCorrectedGenome1 = szCurrentDirectory + "/" + szFreebayesCorrection1 + "/" + os.path.basename( args.szInputGenome ) + ".freebayes1.fa"
szFreebayesCorrectedGenome1Fai = szFreebayesCorrectedGenome1 + ".fai"

# want this full path since we will use it later in a different directory
# for converting coordinates.  Moved here so it has global scope
szFilteredFreebayes1VCF = szCurrentDirectory + "/" + szFreebayesCorrection1 + "/" + "filtered_freebayes1.vcf.gz"












szDoneFlag = szFreebayesCorrection1 + "/freebayes_correction1_done"

if ( not os.path.exists( szDoneFlag ) ):

    os.chdir( szFreebayesCorrection1 )

    szCommand = "set -eo pipefail && gunzip -c " + szFreebayesOnInputGenomeFullPath + " | grep \"^#\" >header.txt"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "set -eo pipefail && gunzip -c " + szFreebayesOnInputGenomeFullPath + " | grep -v \"^#\" | awk '{ if ( length( $5 ) != length( $4 ) ) { {if ( ( $2 - x ) > " + str( nHowFarToLookForCorrelatedIndels ) + " ) print $0; } x = $2 } }' >body.txt"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )


    #szCurrentDirectory = os.getcwd()

    szCommand = "set -eo pipefail && module load tabix/0.2.6  && cat header.txt body.txt | bgzip -c >" + szFilteredFreebayes1VCF
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )


    szCommand = "set -eo pipefail && module load zlib/1.2.11 vcftools/0.1.13 && module load tabix/0.2.6 && tabix " + szFilteredFreebayes1VCF + " && cat " + args.szInputGenome + " | vcf-consensus " + szFilteredFreebayes1VCF + " >" + szFreebayesCorrectedGenome1
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "samtools faidx " + szFreebayesCorrectedGenome1
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    assert os.path.isfile( szFreebayesCorrectedGenome1Fai ), "file " + szFreebayesCorrectedGenome1Fai + " should exist at this point"

    szCommand = "touch freebayes_correction1_done"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )
    

    os.chdir( ".." )

########## aligning reads ###############
# in prepartion for running freebayes, align illumina reads to szFreebayesCorrectedGenome1

print "\nin prepartion for running freebayes, align illumina reads to szFreebayesCorrectedGenome1\n"

szAlignIlluminaSubdirectory  = "alignIlluminaToFreebayesCorrectedGenome1"
szAlignDoneFlag =  "alignIlluminaToFreebayesCorrectedGenome1/alignment_done"

szCurrentDirectory = os.getcwd()
szIlluminaVsFreebayesCorrected1Bam = szCurrentDirectory + "/" + szAlignIlluminaSubdirectory + "/illumina_vs_assembly.sorted.bam"


if ( not os.path.isfile( szAlignDoneFlag ) ):

    szCommand = "mkdir -p " + szAlignIlluminaSubdirectory
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    os.chdir( szAlignIlluminaSubdirectory )

    szCommand = "ln -s -f " + args.szIlluminaPiecesForward
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szIlluminaPiecesReverse = re.sub( r"pieces_forward", "pieces_reverse", args.szIlluminaPiecesForward )

    szCommand = "ln -s -f " + szIlluminaPiecesReverse
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "ln -s -f " + szFreebayesCorrectedGenome1
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "ln -s -f " + szFreebayesCorrectedGenome1Fai
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    # now this is part of the git repository
    #szCommand = "cp ~dgordon/pipelines/align_illumina_against_reference/align_illumina_against_reference.snake ."
    #print "about to execute: " + szCommand
    #subprocess.check_call( szCommand, shell = True )

    szConfigJson = "align_illumina_against_reference_config.json"
    with open( szConfigJson, "w" ) as fConfig:
        fConfig.write( "{\n     \"align_reads_to_this_fasta\": " + "\"" + szFreebayesCorrectedGenome1 + "\",\n     \"tmp_dir\": \"/var/tmp/\" \n}" )


    szCommand = "mkdir -p log"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "source source_this_first.sh && snakemake -s align_illumina_against_reference.snake --drmaa \" -w n -q eichler-short.q -l h_rt=35:00:00 -V -cwd -e ./log -o ./log {params.sge_opts}  -S /bin/bash\"  -w 300 --jobs 100 -p -k --rerun"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )


    szIlluminaVsFreebayesCorrected1Bai = szIlluminaVsFreebayesCorrected1Bam + ".bai"

    assert os.path.exists( szIlluminaVsFreebayesCorrected1Bam ), szIlluminaVsFreebayesCorrected1Bam + "ust exist at this point but doesn't"
    assert os.path.exists( szIlluminaVsFreebayesCorrected1Bai ), szIlluminaVsFreebayesCorrected1Bai + " must exist at this point but doesn't"


    os.chdir( ".." )

    szCommand = "touch " + szAlignDoneFlag
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

########### end of aligning reads ##############

########### run smartie pipeline on corrected genome1 ##############

print "\nrun smartie pipeline on corrected genome1\n"


szSmartie1DoneFlag = szSmartieOnCorrectedGenome1Subdirectory + "/smartie1_done_flag"
if ( not os.path.isfile( szSmartie1DoneFlag ) ):

    szCommand = "mkdir -p " + szSmartieOnCorrectedGenome1Subdirectory
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    os.chdir( szSmartieOnCorrectedGenome1Subdirectory )
    print "chdir to "  + szSmartieOnCorrectedGenome1Subdirectory


    szSmartiePipelineConfigCorrectedGenome1 = "config.json"


    with open( szSmartiePipelineConfigCorrectedGenome1, "w" ) as fSmartiePipelineConfig:
        szFixedPart1 = \
    """\
    {
            "targets" : {
                      "hg38" : "/net/eichler/vol26/eee_shared/assemblies/hg38/legacy/indexes/blasr/ucsc.hg38.no_alts.fasta"
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
            "n_proc" : """
        fSmartiePipelineConfig.write( szFixedPart2 )
        fSmartiePipelineConfig.write( "\"" + str( args.nProcessors ) + "\"\n}\n" )


    szCommand = "SMARTIE_DIR=/net/eichler/vol26/7200/software/pipelines/smartie_sv/201909 && module purge && module load modules modules-init modules-gs/prod modules-eichler/prod miniconda/4.5.12 && module list && mkdir -p log && snakemake -s ${SMARTIE_DIR}/Snakefile -j 10  --jobname \"{rulename}.{jobid}\" --cluster-config ${SMARTIE_DIR}/cluster.config.sge.json --drmaa \" -w n -V -cwd -j y -o ./log -l {cluster.h_rt} -l {cluster.mfree} -pe {cluster.pe} -q {cluster.q} -w n -S /bin/bash\" --verbose -p -k --rerun"

    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )


    # get the indel file from "variants"

    assert os.path.isfile( szSmartieIndelFileCorrectedGenome1 ), szSmartieIndelFileCorrectedGenome1 + " must exist at this point but doesn't"

    os.chdir( ".." )
    szCommand = "touch " + szSmartie1DoneFlag
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )
    

########### end of running smartie pipeline on corrected genome1 ##############


###########  filtering indels in corrected genome 1 ##################

# prepare to run freebayes on the problem locations of corrected1 genome
# to get the problem locations, filter the smartie output

print "\nprepare to run freebayes on the problem locations of corrected1 genome to get the problem locations, filter the smartie output\n"

szRunFreebayesOnCorrection1Genome = "run_freebayes_on_corrected1_genome"

szCommand = "mkdir -p " + szRunFreebayesOnCorrection1Genome
print "about to execute: " + szCommand
subprocess.check_call( szCommand, shell = True )

os.chdir( szRunFreebayesOnCorrection1Genome )

szCurrentDirectory = os.getcwd()
print "current directory: " + szCurrentDirectory

szOutputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1 = szCurrentDirectory + "/outputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1.bed"





szWhereToRunFreebayesDoneFlag = "where_to_run_freebayes_done_flag"
if ( not os.path.isfile( szWhereToRunFreebayesDoneFlag ) ):



    # filter the smartie analysis just as in
    # /net/eichler/vol20/projects/whole_genome_assembly/nobackups/yoruban/fix_pilon/freebayes_polishing4/analysis
    # to get a bed file of the regions to run freebayes on which is done
    # in /net/eichler/vol20/projects/whole_genome_assembly/nobackups/yoruban/fix_pilon/freebayes_polishing5
    # by run_freebayes_in_these_regions.sh

    # finds indels in human coding regions
    szCorrected1GenomeCDS = "hg38-Correct1Genome.cds.bed"
    szCommand = "bedtools intersect -wa -u -a " + szSmartieIndelFileCorrectedGenome1 + " -b " + szHg38CDS + " >" + szCorrected1GenomeCDS
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )


    # find indels (with respect to human) that are gene-killing (not a
    # multiple of 3) in human coding regions and other filters.  I don't
    # believe there is any guarantee that freebayes will find any Illumina
    # reads that would correct such indels.  But we will try to find such
    # by running freebayes at these locations.

    szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigs = "indel.cds.noseg-nolowc-noends-nosmallcontigs.bed"

    szCommand = "set -eo pipefail && perl -lane 'print if $F[4] % 3 != 0' " + szCorrected1GenomeCDS + " | sort -k1,1 -k2,2n | bedtools subtract -A -a - -b " + szHg38LowConfidenceRegions + " | bedtools subtract -A -a - -b " + szHg38WgacSuperDup + " | perl -lane 'print if ($F[7] > " + str( nHowCloseToEnds ) +  " && ($F[9] - $F[7]) > " + str( nHowCloseToEnds ) + " && $F[9] > " + str( nMinimumSizeContig ) + " ) ' | sort -k1,1 -k2,2n > " + szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigs
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace = "indel.cds.noseg-nolowc-noends-nosmallcontigs.falconspace.bed"

    # convert to falcon coordinates:
    # (can do this easily because the smartie pipeline gave both human and falcon coordinates
    # for each indel)
    szCommand = "awk 'BEGIN {OFS = \"\t\" } {print $7,$8,$9,$0}' " + szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigs + " | sort -k1,1 -k2,2n >" +  szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    # remove low and high depth regions

    # convert high and low depth regions into coordinates of the 
    # freebayes_corrected1 genome

    szCommand = "ln -s ../convertCoordinates2.py "
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    # change high/low depth from falcon coordinates *before* this pipeline ran to 
    # falcon coordinate of correctedGenome1 coordinates.  Does this using the vcf file
    # used to make correctedGenome1

    # szOutputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1 = szCurrentDirectory + "/outputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1.bed"
    szCommand = "./convertCoordinates2.py --szInputHighAndLowDepthBedFile " + szHighAndLowDepthRegionsFalconCoordsBedFullPath  + " --szFreebayesVCFFile " + szFilteredFreebayes1VCF + " --szOutputHighAndLowDepthBedFileWithRespectToNewGenome " + szOutputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1 + " --szNewGenomeFaiFile " + szFreebayesCorrectedGenome1Fai
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szFilteredGeneKillingIndels = "for_david_final.repair.bed"

    szCommand = "bedtools subtract -a " + szGeneKillingIndelsInCorrectedGenome1NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace + " -b " + szOutputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1 + " -A  > " + szFilteredGeneKillingIndels
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    # remove neighboring events in which ins/del cancel each other or in which ins/ins or 
    # del/del sum to multiple of 3

    szCommand = "ln -s ../find_nongenekilling_indels4.py"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    # input to find_nongenekilling_indels4.py is for_david_final.repair.bed
    # output is pairs_of_indels_removed4.bed

    szCommand = "./find_nongenekilling_indels4.py"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szRemainingIndelsCorrectedGenome1 = "pairs_of_indels_removed4.bed"

    assert os.path.isfile( szRemainingIndelsCorrectedGenome1 ), szRemainingIndelsCorrectedGenome1 + " must exist at this point but doesn't"

    print "number of filtered gene-killing indels:"
    szCommand = "wc -l " + szRemainingIndelsCorrectedGenome1
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "ln -s -f " + szFreebayesCorrectedGenome1Fai
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szRunFreebayesInTheseRegions = "run_freebayes_in_these_regions.bed"

    szCommand = "set -eo pipefail && bedtools slop -b " + str( nRunFreebayesPlusMinusThisDistanceFromGenekillingIndels ) + " -g " + szFreebayesCorrectedGenome1Fai + " -i " + szRemainingIndelsCorrectedGenome1 + " | sort -k1,1 -k2,2n | bedtools merge -i stdin >" + szRunFreebayesInTheseRegions
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    assert os.path.isfile( szRunFreebayesInTheseRegions ), szRunFreebayesInTheseRegions + " should exist at this point but doesn't"

    szCommand = "touch " + szWhereToRunFreebayesDoneFlag
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

###########  end filtering indels in corrected genome 1 ##################

########### run freebayes at pin-point locations in corrected genome 1 ############

print "\nrun freebayes at pin-point locations in corrected genome 1\n"

szFreebayesVCFOnCorrectedGenome1 = szCurrentDirectory + "/final/merged.vcf.gz"

szRunFreebayesDoneFlag = "run_freebayes_done_flag"
if ( not os.path.isfile( szRunFreebayesDoneFlag ) ):

    szConfig = "config.yaml"
    with open( szConfig, "w" ) as fConfig:
        fConfig.write( "---\n\n" )
        fConfig.write( "fasta: " + szFreebayesCorrectedGenome1 + "\n" )
        fConfig.write( "bamlist: correctedGenome1.bamlist.txt\n" )
        fConfig.write( "exclude: []\n" )


    with open( "correctedGenome1.bamlist.txt", "w" ) as fBamList:
        fBamList.write( szIlluminaVsFreebayesCorrected1Bam + "\n" )

    # create data/regions.txt

    szRegionsFile = "data/regions.txt"

    szCommand = "mkdir -p data"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )
      
    szCommand =  "awk '{print $1\":\"$2\"-\"$3}' run_freebayes_in_these_regions.bed >" + szRegionsFile
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand =  "cp ~dgordon/pipelines/running_freebayes_on_regions/env.cfg ."
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand =  "cp ~dgordon/pipelines/running_freebayes_on_regions/freebayes_on_regions*.snake ."
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand =  "cp ~dgordon/pipelines/running_freebayes_on_regions/run_snakemake.sh ."
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCurrentDirectory = os.getcwd()
    print "current directory: " + szCurrentDirectory

    szCommand = "module unload python && module list && run_snakemake.sh"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    #szFreebayesVCFOnCorrectedGenome1 = szCurrentDirectory + "/final/merged.vcf.gz"
    assert os.path.isfile( szFreebayesVCFOnCorrectedGenome1 ), szFreebayesVCFOnCorrectedGenome1 + " should exist at this point but doesn't"

    szCommand = "touch " + szRunFreebayesDoneFlag
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )


os.chdir( ".." )
########### end run freebayes at pin-point locations in corrected genome 1 ############


####################   create correctedGenome2.fa, our final genome ##########################
print "\ncreate correctedGenome2.fa, our final genome\n"

szCorrectedGenome2Dir = "freebayes_correction2"


szCurrentDirectory = os.getcwd()
print "current directory: " + szCurrentDirectory

szCorrectedGenome2 = szCurrentDirectory + "/" + szCorrectedGenome2Dir + "/correctedGenome2.fa"


szMakeCorrectedGenome2 = szCorrectedGenome2Dir + "/make_corrected_genome2_done"

if ( not os.path.isfile( szMakeCorrectedGenome2 ) ):

    szCommand = "mkdir -p " + szCorrectedGenome2Dir
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    os.chdir( szCorrectedGenome2Dir )

    # apply these indels to create correctedGenome2.fa  (In the future, should exclude snps.)


    szCommand = "set -eo pipefail && module load zlib/1.2.11 vcftools/0.1.13 && module load tabix/0.2.6 && cat " + szFreebayesCorrectedGenome1 + " | vcf-consensus " + szFreebayesVCFOnCorrectedGenome1 + " >" + szCorrectedGenome2
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    # check that correctedGenome2.fa exists

    assert os.path.exists( szCorrectedGenome2 ), szCorrectedGenome2 + " should exist at this point but doesn't"

    szCommand = "samtools faidx " + szCorrectedGenome2
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    os.chdir( ".." )

    szCommand = "touch " + szMakeCorrectedGenome2
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

szIndelCorrectedDoneFlag = "indel_corrected_done"
if ( not os.path.isfile( szIndelCorrectedDoneFlag ) ):

    szCommand = "ln -s -f " + szCorrectedGenome2 + " almost_done.fa"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "module load python/2.7.13 numpy/1.16.6 biopython/1.71 &&  ./rename_sequences.py --szInputFasta almost_done.fa --szOutputFasta indel_corrected.fa"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "chmod a-w indel_corrected.fa"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )


    szCommand = "samtools faidx indel_corrected.fa"
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    szCommand = "touch " + szIndelCorrectedDoneFlag
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )


################ done creating correctedGenome2.fa  (In the future, should exclude snps.) ####################


################   run smartie pipeline again to see how many gene-killing indels remain ##########################################################

print "\nrun smartie pipeline again to see how many gene-killing indels remain\n"



szSmartieOnCorrectedGenome2Done = szSmartieOnCorrectedGenome2Subdirectory + "/smartie2_done"
if ( not os.path.isfile( szSmartieOnCorrectedGenome2Done ) ):

    szCommand = "mkdir -p " + szSmartieOnCorrectedGenome2Subdirectory
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )


    os.chdir( szSmartieOnCorrectedGenome2Subdirectory )
    print "chdir to "  + szSmartieOnCorrectedGenome2Subdirectory


    szSmartiePipelineConfigCorrectedGenome2 = "config.json"

    with open( szSmartiePipelineConfigCorrectedGenome2, "w" ) as fSmartiePipelineConfig:
        szFixedPart1 = \
    """\
    {
            "targets" : {
                      "hg38" : "/net/eichler/vol26/eee_shared/assemblies/hg38/legacy/indexes/blasr/ucsc.hg38.no_alts.fasta"
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
            "n_proc" : """
        fSmartiePipelineConfig.write( szFixedPart2 )
        fSmartiePipelineConfig.write( "\"" + str( args.nProcessors ) + "\"\n}\n" )


    szCommand = "SMARTIE_DIR=/net/eichler/vol26/7200/software/pipelines/smartie_sv/201909 && module purge && module load modules modules-init modules-gs/prod modules-eichler/prod miniconda/4.5.12 && module list && mkdir -p log && snakemake -s ${SMARTIE_DIR}/Snakefile -j 10  --jobname \"{rulename}.{jobid}\" --cluster-config ${SMARTIE_DIR}/cluster.config.sge.json --drmaa \" -w n -V -cwd -j y -o ./log -l {cluster.h_rt} -l {cluster.mfree} -pe {cluster.pe} -q {cluster.q} -w n -S /bin/bash\" --verbose -p -k --rerun"

    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )

    assert os.path.isfile( szSmartieIndelFileCorrectedGenome2 ), szSmartieIndelFileCorrectedGenome2 + " must exist at this point but doesn't"

    os.chdir( ".." )


    szCommand = "touch " + szSmartieOnCorrectedGenome2Done
    print "about to execute: " + szCommand
    subprocess.check_call( szCommand, shell = True )
    
# end run smartie pipeline again to see how many gene-killing indels remain
##########################################################################



# now filter the smartie indels and see how many remain
print "\nnow filter the smartie indels and see how many remain\n"

szCorrected2GenomeCDS = "hg38-Correct2Genome.cds.bed"
szCommand = "bedtools intersect -wa -u -a " + szSmartieIndelFileCorrectedGenome2 + " -b " + szHg38CDS + " >" + szCorrected2GenomeCDS
print "about to execute: " + szCommand
subprocess.check_call( szCommand, shell = True )

szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigs = "indel.cds.noseg-nolowc-noends-nosmallcontigs.bed"

szCommand = "set -eo pipefail && perl -lane 'print if $F[4] % 3 != 0' " + szCorrected2GenomeCDS + " | sort -k1,1 -k2,2n | bedtools subtract -A -a - -b " + szHg38LowConfidenceRegions + " | bedtools subtract -A -a - -b " + szHg38WgacSuperDup + " | perl -lane 'print if ($F[7] > " + str( nHowCloseToEnds ) + " && ($F[9] - $F[7]) > " + str( nHowCloseToEnds ) + " && $F[9] > " + str( nMinimumSizeContig ) + " ) ' | sort -k1,1 -k2,2n > " + szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigs
print "about to execute: " + szCommand
subprocess.check_call( szCommand, shell = True )

szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace = "indel.cds.noseg-nolowc-noends-nosmallcontigs.falconspace.bed"

# convert to falcon coordinates:
szCommand = "set -eo pipefail && awk 'BEGIN {OFS = \"\t\" } {print $7,$8,$9,$0}' " + szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigs + " | sort -k1,1 -k2,2n >" +  szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace
print "about to execute: " + szCommand
subprocess.check_call( szCommand, shell = True )

# remove low and high depth regions

szFilteredGeneKillingIndels = "for_david_final.repair.bed"

szCommand = "bedtools subtract -a " + szGeneKillingIndelsInCorrectedGenome2NoSegNoLowConfNoContigEndsNoSmallContigsFalconSpace + " -b " + szOutputHighAndLowDepthBedFileWithRespectToFreebayesCorrectedGenome1 + " -A  > " + szFilteredGeneKillingIndels
print "about to execute: " + szCommand
subprocess.check_call( szCommand, shell = True )

# remove neighboring events in which ins/del cancel each other or in which ins/ins or 
# del/del sum to multiple of 3

# szCommand = "cp ~dgordon/pipelines/freebayes_polishing/find_nongenekilling_indels4.py ."
# print "about to execute: " + szCommand
# subprocess.check_call( szCommand, shell = True )

szCommand = "./find_nongenekilling_indels4.py"
print "about to execute: " + szCommand
subprocess.check_call( szCommand, shell = True )


szRemainingIndelsCorrectedGenome2 = "pairs_of_indels_removed4.bed"

assert os.path.isfile( szRemainingIndelsCorrectedGenome2 ), szRemainingIndelsCorrectedGenome2 + " must exist at this point but doesn't"

# this is now subsumed by /find_nongenekilling_indels4.py
# print "number of filtered gene-killing indels:"
# szCommand = "wc -l " + szRemainingIndelsCorrectedGenome2
# print "about to execute: " + szCommand
# subprocess.check_call( szCommand, shell = True )

# szCommand = "echo  \"number of filtered gene-killing indels:\" >gene_killing_indels.txt"
# print "about to execute: " + szCommand
# subprocess.check_call( szCommand, shell = True )

# szCommand = "wc -l " + szRemainingIndelsCorrectedGenome2 + " >>gene_killing_indels.txt"
# print "about to execute: " + szCommand
# subprocess.check_call( szCommand, shell = True )






    
    


                                      
    

