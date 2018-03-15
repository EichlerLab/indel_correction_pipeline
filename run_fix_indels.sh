export PATH=.:$PATH
source source_this_first.sh
fix_indels.py --szInputGenome /net/eichler/vol27/projects/assemblies/nobackups/macaque/download_stlouis/macaque.fa \
--szInputFreebayesVCFFile /net/eichler/vol27/projects/assemblies/nobackups/macaque/initial_freebayes/final/merged.macaque.vcf.gz \
--szPrefixForSmartie mac180309 \
--nProcessors 22 \
--szIlluminaPiecesForward /net/eichler/vol27/projects/assemblies/nobackups/macaque/bwa_illumina_vs_pilond/pieces_forward \
--szHighAndLowDepthRegionsBed high_low_depth_2.7_99_falcon_coordinates.bed
