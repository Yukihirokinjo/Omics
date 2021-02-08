#!/bin/bash

##
## RNAout2gff.sh
##
## USAGE: RNAout2gff.sh  <infernal output file(fm2)>  <tRENAScan-SE output file>  <RNAmmer output file>
##
##

InfernalOUT=$1
tRNAscanOUT=$2
RNAmmerOUT=$3
Name=$4

## Infernal(fm2)
grep -v "#\|tRNA\|rRNA\|?" ${InfernalOUT} | \
awk 'BEGIN {OFS="\t"} {if ($2 == "tmRNA") print $4, "Infernal", "tmRNA", $10, $11, $16, $12, ".", $2; \
     else print $4, "Infernal", "ncRNA", $10, $11, $16, $12, ".", $2}' > tmp_pInfernalOUT.gff
awk 'BEGIN {OFS="\t"} {if ($7 == "+") print  $0; \
     else print $1, $2, $3, $5, $4, $6, $7, $8, $9}' tmp_pInfernalOUT.gff > tmp_InfernalOUT.gff

## tRNAScan-SE
# NOTE: need to reverse "start($3)" and "end($4)" when tRNAs are coded in reverse strand; for proper gff format
grep "0" ${tRNAscanOUT} | \
awk 'BEGIN {OFS="\t"} {if ($3 < $4) print $1, "tRNAScan-SE", "tRNA", $3, $4, $9, "+", ".", "tRNA_type=", $5, "Anti_codon=", $6; \
       else if ($3 > $4) print $1, "tRNAScan-SE", "tRNA", $4, $3, $9, "-", ".", "tRNA_type=", $5, "Anti_codon=", $6}'  > tmp1_tRNAscanOUT.txt
sed 's/=\t/=/g' tmp1_tRNAscanOUT.txt | sed 's/\tAnti_codon/;Anti_codon/g' > tmp_tRNAscanOUT.gff

## RNAmmer
grep -v "#" ${RNAmmerOUT} | sed 's/ /_/g' | sed 's/:/-/g' > tmp_RNAmmerOUT.gff


## Merge
cat  tmp_InfernalOUT.gff tmp_tRNAscanOUT.gff tmp_RNAmmerOUT.gff |  sort -k 1,1 -k 4n,4 > ${Name}_RNAoutMerge.gff

rm ./tmp_*OUT.gff
rm ./tmp1_*OUT.txt


