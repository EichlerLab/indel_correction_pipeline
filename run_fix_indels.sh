export PATH=.:$PATH
source source_this_first.sh
fix_indels_start_earlier.py --szInputGenome /net/eichler/vol27/projects/assemblies/nobackups/marmoset/pilon/marmoset_arrow_pilon.fa \
--szPrefixForSmartie marm \
--nProcessors 15 \
--szIlluminaPiecesForward /net/eichler/vol27/projects/assemblies/nobackups/marmoset/illumina_reads/split/pieces_forward \
--szPacBioReadsFof /net/eichler/vol27/projects/assemblies/nobackups/marmoset/pacbio_reads/input_pacbio.fofn
