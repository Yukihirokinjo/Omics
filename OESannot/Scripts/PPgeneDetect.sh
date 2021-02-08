#!/bin/bash

##
## PPgeneDetect.sh
##
## USAGE: PPgeneDetect.sh  <query_genome.gff> <query_genome.fna> <ref.proteome.faa>  <output_prefix>

if [ -e "PPgene" ]; then
 rm -rf PPgene
fi
mkdir  PPgene

cat $1 | awk '$3 ~ /CDS/ {print $0}'  > ./PPgene/tmp_input.gff
cat $2 > ./PPgene/tmp_input.fna
cat $3 > ./PPgene/tmp_ref.faa
Name=$4
numCPU=$5
COGinfo=$6
gCode=$7

cd ./PPgene

Rscript $(which gff2faa.R)  tmp_input.gff  tmp_input.fna  ${Name}  ${gCode:=11}
cat gff2.faa > tmp_input.faa

cp ${COGinfo} ./tmp_COGinfo.csv

### 1. serch >2 continuous imcomplete fragments which hit to same regerence protein

diamond makedb --in tmp_ref.faa -d tmp_ref.faa 1> dmdDB.log
diamond blastp --sensitive  -d tmp_ref.faa -q tmp_input.faa -k 1 --max-hsps 1 -e 0.1  -o blast.out -p ${numCPU:=11} 1> dmd.log

#makeblastdb -in tmp_ref.faa  -dbtype prot -hash_index
#blastp -db tmp_ref.faa -query tmp_input.faa -outfmt 6 -out blast.out -max_target_seqs 1 -max_hsps 1 -evalue 0.1 -num_threads ${numCPU:=1}
rm ./tmp_ref.faa.*
awk '{print $1, $2, $11, $9, $10}' blast.out | uniq > tmp_bla.out

Rscript $(which BlastMultiHitExtract.R)  tmp_bla.out tmp_COGinfo.csv


Rscript $(which PPgeneDetect.R)   tmp_BlastMultiHit_Continuous.txt  tmp_input.fna ${gCode:=11}


if [ -e "PPrecode.faa" ]; then
  mv PPrecode.fna  ${Name}_PPrecode.fna
  mv PPrecode.faa  ${Name}_PPrecode.faa
  cp  ${Name}_PPrecode.faa  ..
fi

mv PPgene.gff  ${Name}_PPgene.gff
#mv PPgene.txt  ${Name}_PPgene.txt
mv PseudoGeneCandi.gff ${Name}_PseudoGeneCandi.gff
cat ${Name}_PseudoGeneCandi.gff ${Name}_PPgene.gff  > PPgeneDetect.gff
mv PPgeneDetect.gff ${Name}_PPgeneDetect.gff
cp ${Name}_PPgeneDetect.gff  ..

mkdir ./tmp
mv ./tmp_*  ./tmp
rm  ./tmp/tmp_ref.faa

cd ../

##END
