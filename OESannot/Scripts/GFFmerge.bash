#!/bin/bash

gff1=$1
gff2=$2
Name=$3

## priority: Pseudo > ShortORF > CDS

cat ${gff1}  ${gff2} | grep -v "#" | sed 's/:/-/g' | awk '$3 ~ /CDS|pseudo|PPgene/ {OFS="\t"; print $0}' | sort -k 1,1 -k 4n,5 | sed 's/  /\t/g' | uniq > tmp_premerge.gff
cat ${gff1}  ${gff2} | grep -v "#" | sed 's/:/-/g' | awk '$3 ~ /RNA/ {OFS="\t"; print $0}' | sort -k 1,1 -k 4n,5  | sed 's/  /\t/g' | uniq > tmp_RNA.gff


 Rscript  $(which GFFmerge.R)

if [ -e tmp_insideORF.gff ]; then
  cat tmp_merged.gff tmp_RNA.gff | uniq | sort -k 1,1 -k 4n,5  >  ${Name}_merge.gff
else
  cat tmp_premerge.gff tmp_RNA.gff | uniq | sort -k 1,1 -k 4n,5  >  ${Name}_merge.gff
  echo "No ORF duplications were found."
fi
rm  tmp_*.gff
