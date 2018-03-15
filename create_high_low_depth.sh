set -v
awk 'BEGIN {OFS = "\t" } {if ( ($4 < 2.7) || (98.66 < $4 ) ) print $1, $2, $3, $4 }' /net/eichler/vol27/projects/assemblies/nobackups/macaque/read_depth_by_position/merged_smooth/1000/all.smooth.depth.1000.bed >high_low_depth_2.7_99_falcon_coordinates.bed
