set -v
awk 'BEGIN {OFS = "\t" } {if ( ($5 < 2.7) || (170 < $5 ) ) print $1, $2, $3, $5 }' /net/eichler/vol18/zevk/great_apes/contig_breaking/clint/depth/merged_smooth/1000/all.smooth.depth.1000.txt >high_low_depth_2.7_170_falcon_coordinates.bed
