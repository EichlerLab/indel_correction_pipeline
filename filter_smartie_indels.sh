set -ve

#szInput="hg38-YRI_fix2.indel.bed"
szInput=$1

bedtools intersect -wa -u -a $szInput -b ~/eee_shared/assemblies/hg38/genes/refGene.CDS_exons.bed  >hg38-yoruban_freebayes5.indel.cds.bed


perl -lane 'print if $F[4] % 3 != 0' hg38-yoruban_freebayes5.indel.cds.bed | sort -k1,1 -k2,2n | bedtools subtract -A -a - -b hg38_low_confidence_regions.bed | bedtools subtract -A -a - -b ~/eee_shared/assemblies/hg38/wgac/genomicSuperDup.sort.bed.gz | perl -lane 'print if ($F[7] > 10000 && ($F[9] - $F[7]) > 10000 && $F[9] > 500000) ' | sort -k1,1 -k2,2n > yri_david_fix2.no-seg-no-lowc.200K.bed

# this takes more than 15 minutes
#bedtools coverage -a hg38.windows.bed -b yri_david_fix2.no-seg-no-lowc.200K.bed  > yri_david_fix2.no-seg-no-lowc.200K.1Mb.cov

# convert to falcon coordinates:
awk 'BEGIN {OFS = "\t" } {print $7,$8,$9,$0}' yri_david_fix2.no-seg-no-lowc.200K.bed | sort -k1,1 -k2,2n >yri_david_fix2.no-seg-no-lowc.200K.yri.space.bed


perl -lane 'print if $F[-1] < 10 || $F[-1] > 150' /net/eichler/vol20/projects/whole_genome_assembly/nobackups/yoruban/aligned_read_depth3/merged_smooth/1000/all.smooth.depth.1000.bed > yri-lt5-gt199.bed


bedtools subtract -a yri_david_fix2.no-seg-no-lowc.200K.yri.space.bed -b yri-lt5-gt199.bed -A  > for_david_final.repair.bed

echo "number of indels not multiple of 3 in yoruban"
wc for_david_final.repair.bed

# convert to hg38 coordinates

#awk 'BEGIN {OFS = "\t" } { print $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17 }' for_david_final.repair.bed > for_david_final.repair.falconspace.bed


# how many are in coding sequence

#bedtools intersect -wa -u -a for_david_final.repair.falconspace.bed -b ~/eee_shared/assemblies/hg38/genes/refGene.CDS_exons.bed  >for_david_final.repair.cds.bed

echo "number of indels not multiple 3 in coding region:"
wc for_david_final.repair.bed
