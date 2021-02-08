#!/bin/bash
#
# gff2faa.sh
#

GFF=$1
FNA=$2
Name=$3
TAG=$4
gCode=$5

grep -v ^"#"  ${GFF}  > ${GFF}_gff2faa.gff

 Rscript $(which gff2faa.R) ${GFF}_gff2faa.gff ${FNA} ${TAG}  ${gCode:=11}

cat gff2.faa > ${Name}.faa
cat gff2RNA.fna > ${Name}_RNAs.fna

rm ${GFF}_gff2faa.gff
rm gff2.faa
rm gff2RNA.fna
