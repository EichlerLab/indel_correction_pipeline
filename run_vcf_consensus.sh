set -v

module load VCFtools/0.1.12b
module load tabix/0.2.6

cat quivered.pilon.fa | vcf-consensus altered.vcf.gz >quivered.pilon.freebayes.fa
samtools faidx quivered.pilon.freebayes.fa
