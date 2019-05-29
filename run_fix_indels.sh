export PATH=.:$PATH
source source_this_first.sh
fix_indels_start_earlier.py --szInputGenome /net/eichler/vol27/projects/autism_genome_assembly/nobackups/assemblies/proband/squashed_pilon/squashed_arrowed_piloned.fa \
--szPrefixForSmartie mother190515 \
--nProcessors 22 \
--szIlluminaPiecesForward /net/eichler/vol27/projects/autism_genome_assembly/nobackups/assemblies/proband/illumina_reads/cram_to_fasta/fastqs/pieces_forward \
--szHighAndLowDepthRegionsBed /net/eichler/vol27/projects/autism_genome_assembly/nobackups/assemblies/proband/squashed_indel/read_depth_by_position/high_and_low_depth_regions.bed
