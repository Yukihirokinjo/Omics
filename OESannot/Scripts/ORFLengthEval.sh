#!/bin/bash
##
## ORFLengthEval.sh
##
## USAGE: ORFLengthEval.sh  [input.gff] [input.fna] [blastDB] [DB_length.txt] [output_prefix] [Num.CPU]

INgff=$1
INfna=$2
blaDB=$3
refInfo=$4
refLen=$5
Name=$6
numCPU=$7
Prop=$8
gCode=$9

if [ -e "ORFLen_${Name}" ]; then
 rm -rf ORFLen_${Name}
fi

mkdir ORFLen_${Name}
cat ${INgff}  >  ORFLen_${Name}/tmp_input.gff
cat ${INfna}  > ORFLen_${Name}/tmp_input.fna
cat ${blaDB}  >  ORFLen_${Name}/tmp_DB.faa
cat ${refInfo} > ORFLen_${Name}/tmp_COGinfo.csv
cat ${refLen}  > ORFLen_${Name}/tmp_COGLength.txt
cd ORFLen_${Name}

## R: Obtain sequence length information of reference DB

#  R --vanilla --slave  --args  tmp_DB.faa  RefLen.txt  <  `which GetDBseqLength.R`


## Generate amino acid sequence file from input gff file

gff2faa.sh  tmp_input.gff tmp_input.fna tmp_input  tmp  ${gCode:=11}

diamond makedb --in tmp_DB.faa -d tmp_DB 1> dmdDB.log
diamond blastp --sensitive  -d tmp_DB -q tmp_input.faa -k 1 --max-hsps 1 -e 0.1  -o tmp_ORFLen_bla.out -p ${numCPU:=1} 1> dmd.log
#makeblastdb -in tmp_DB.faa  -dbtype prot -hash_index  -out tmp_DB
#blastp -db tmp_DB  -query tmp_input.faa -outfmt 6 -out tmp_ORFLen_bla.out  -max_target_seqs 1 -max_hsps 1 -evalue 0.1 -num_threads ${numCPU:=1}

## R: Evaluate sequence length

  Rscript $(which ORFLengthEval.R) tmp_ORFLen_bla.out  tmp_input.faa  tmp_COGinfo.csv   tmp_COGLength.txt  tmp_DB.faa  ${Prop}

##

#cat BlastHitLengthCheck.txt > ../BlastHitLengthCheck_${Name}.txt
cat  tmp_ShortEval.gff > ../ShortEval_${Name}.gff
rm ./tmp_DB.faa
rm ./tmp_COGinfo.csv

cd ..

## END OF FILE
