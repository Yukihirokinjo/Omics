#!/bin/bash
#
#IGSextract.bash
#
#.gff ファイルからの遺伝子間領域の抜きだし
#
#
#2015_0304
#
#Usage: IGSextract.bash  <input.gff>  <queryGenome.fna>  <outname>

infile=$1
query=$2
outname=$3

mkdir IGS
cat ${infile} > ./IGS/tmp_input.gff
cat ${query}  > ./IGS/tmp_query_genome.fna
cd ./IGS

awk '$3 ~ /RNA|CDS|pseudo|PPgene/ {OFS="\t"; print $1,$4,$5}' tmp_input.gff > tmp_genepos.tsv

#R------------------------------------------------------

Rscript $(which IGSeq.R)

#R------------------------------------------------------

mv ./IGS_Seq.fna  ../IGS_Seq_${outname}.fna
mv ./IGSeq.gff    ../IGS_${outname}.gff
rm ./tmp_genepos.tsv
rm ./tmp_query_genome.fna

cd ..
