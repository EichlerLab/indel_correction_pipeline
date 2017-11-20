set -v

module load tabix/0.2.6
szInputFile="clint.vcf.gz"

gunzip -c $szInputFile | grep "^#"  >header.txt
gunzip -c $szInputFile | grep -v "^#"  | awk '{if ( length($4) != length( $5) ) print $0 }' | awk '{ if ( $2 - x > 6 ) print $0; x = $2 }' >body.txt
cat header.txt body.txt | bgzip -c >altered.vcf.gz
