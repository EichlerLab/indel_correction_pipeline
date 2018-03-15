see wiki:

https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=users:indel_correction_pipeline_using_freebayes

which has more detailed information, especially for the previous
steps.


This pipeline does the following:


1.Run FreeBayes (Garrison E, et al. 2012. arXiv preprint) of genome*
Illumina reads aligned against genome reference giving a list of
variants 

2.In cases in which there are 2 indels ¡Ü6 bp apart, remove the 2nd
from the list and remove non-indels from list.

3.Stage 1 correction: Change reference to the variants in list from #2
(genome-wide) using vcf-consensus (The Variant Call Format and
VCFtools, Petr Danecek, Adam Auton, Goncalo Abecasis, Cornelis
A. Albers, Eric Banks, Mark A. DePristo, Robert Handsaker, Gerton
Lunter, Gabor Marth, Stephen T. Sherry, Gilean McVean, Richard Durbin
and 1000 Genomes Project Analysis Group, Bioinformatics, 2011) giving
new reference: ¡°Corrected FreeBayes 1

4.Find remaining high-confidence gene-killing indels by applying
filters that eliminated lower confidence (indels in regions of
segmental duplication, high depth regions of the assembly, telomeric,
centromeric, small contigs (less than 500K), and mapping within 10kb
of the ends of contigs). 
Empirically we found most of the remaining gene killing indels in
CorrectedFreeBayes1 were heterozygous. To change them to the benign
variant that restored frame, we follow the following steps:

5.Align Illumina reads against CorrectedFreeBayes1

6.Run FreeBayes (again) with the genome Illumina reads against
CorrectedFreeBayes1 but just at the sites (+/-100 bp) of indels from
#4

7.Stage 2 correction: Change the reference to the variants in #6 (not
just indels, 215 in chimp) using vcf-consensus giving
¡°CorrectedFreeBayes2¡±--the output of this pipeline.

*Yoruban, Clint-chimpanzee, CHM13 or Susie-orangutan



> less run_fix_indels.sh
module load python/2.7.3
fix_indels.py --szInputGenome
~/eee_shared/assemblies/clint_panTro/polished_quivered_masked/clint.contigs.quiver.pilon.masked.fasta
--szInputFreebayesVCFFile
/net/eichler/vol25/projects/whole_genome_assembly/nobackups/chimp/Drafts/V3.0a-falcon/freebayes_polishing/clint.vcf.gz
--szPrefixForSmartie ylint170309 
--nProcessors 22
--szIlluminaPiecesForward
/net/eichler/vol25/projects/whole_genome_assembly/nobackups/chimp/Drafts/V3.0a-falcon/freebayes_polishing5/illumina_reads/pieces_forward
--szHighAndLowDepthRegionsBed
high_low_depth_2.7_170_falcon_coordinates.bed

where: 
--szInputGenome
is the fasta file with indel errors

--szInputFreebayesVCFFile
is the output of freebayes of Illumina reads run against the fasta
file above

--szPrefixForSmartie 
is (I believe) an arbitrary prefix that the smartie
pipeline needs to output and intermediate files

--nProcessors
is the number of processors used by blasr of the corrected genomes
(stage 1 and stage 2) against human as part of the smartie pipeline

--szIlluminaPiecesForward
is the directory of fastq files each with lots of forward Illumina
reads.  the path must end with "pieces_forward" which is replaced by
"pieces_reverse" to get to the directory with the corresponding files
containing Illumina reverse files.  The purpose is to parallelize
running bwa of these reads against the stage1 and stage2 freebayes
corrected genomes.  Since the input of this script is a freebayes vcf
file, there is no need to run bwa against the original (uncorrected)
file.

--szHighAndLowDepthRegionsBed
is a bed file of high and low depth regions in the coordinates of the
uncorrected genome.  As such it is useless unless converted to
coordinates of the stage1 correction genome.  This is done by
convertCoordinates2.py 







