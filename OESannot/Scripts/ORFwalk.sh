#!/bin/bash

## ORFwalk.sh
## Extend possible longer flames for existing ORFs
##
## USAGE: ORFwalk <orfs.gff>  <genome.fna> <overlap(nt)> <5/3 Edge(5 or 3)>
##


Gff=$1
Fna=$2
Ext=$3
Name=$4
gCode=$5


  Rscript $(which ORFwalk.R)  "${Gff}" "${Fna}" ${Ext} ${gCode:=11} 


cat ORFw.gff  >  ${Name}_ORFw.gff

rm  ORFw.gff
