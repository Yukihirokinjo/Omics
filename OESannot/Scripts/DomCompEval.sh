#!/bin/bash

## DomCompEval.sh
##
## USAGE: DomCompEval.sh [ input.faa ] [PSSM_database] [CDDid file] [output prefix]
##
## Note: using only Pfam and CDD_ncbi PSSM database. (constracted as "CDDncbi_Pfam.pn" in this script)
##
## /db/cdd/CDDncbi_Pfam

INgff=$1
INfna=$2
CDD=$3
CDDid=$4
Name=$5
numCPU=$6
Prop=$7
#COG=$8

if [ -e "DomComp_${Name}" ]; then
 rm -rf DomComp_${Name}
fi

mkdir DomComp_${Name}
cat ${INgff} >  DomComp_${Name}/tmp_input.gff
cat ${INfna} > DomComp_${Name}/tmp_input.fna


## prerequisite:
#
# 1. construct profile db (cdd_ncbi + pfam)
#
#   cat Cdd_NCBI.pn Pfam.pn > CDDncbi_Pfam.pn
#   makeprofiledb  -in CDDncbi_Pfam.pn -out CDDncbi_Pfam  -threshold 9.82 -scale 100.0 -dbtype rps -index true
#
# 2. extract id and length information from cddid.tbl file
#
   awk '{print $1,$2,$3,$NF}' ${CDDid} > ./DomComp_${Name}/cddid_extract.tbl
#

cd ./DomComp_${Name}

gff2faa.sh  tmp_input.gff tmp_input.fna tmp_input  tmp

# cd-search by rps-blast

rpsblast -db ${CDD} -query tmp_input.faa -evalue 1e-5  -outfmt 6 -out RPS_custom6  -num_threads ${numCPU:=1} 1> rpsblast.log
#rpsblast -db ${COG} -query tmp_input.faa -evalue 1e-5  -outfmt 6 -out RPS_custom6  -num_threads ${numCPU:=1}

## R

  Rscript $(which DomCompEval.R)  RPS_custom6  tmp_input.faa  cddid_extract.tbl ${Prop}  

##


# GFFmerge.bash  tmp_input.gff  DomEval.gff  ${Name}_LengthDomEval

cp DomEval.gff  ../${Name}_DomEval.gff
cd ..


## omit
## pfam: hmmscan
#   hmmscan --cpu 8 --domE 1e-5 --domtblout pfam.domtblout2 /db/pfam/Pfam-A.hmm  PPtest.faa
#   cat  pfam.domtblout2  | grep -v ^"#" | awk '{for(i=1;i<22;i++){printf $i "\t"} print""}' >  tmp_pfamhmscn.txt
#
