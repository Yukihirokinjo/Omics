#!/bin/bash

##
## PseudoGeneDetect.sh
##
## USAGE: PseudoGeneDetect.sh  <input.gff> <Name>
#
#
# - SingleHit:Short		->
# - MultiHitConti:notPPcandi	-> PseudoGeneCandi.gff(1st)
# - MultiHitConti:PPcandi	-> PseudoGeneCandi.gff(1st+2nd)
#
# -> PseudoGeneCandi_ALL.gff

InGFF=$1
Name=$2

mkdir Pseudo
cp ${InGFF} Pseudo/input.gff
cd Pseudo

R --vanilla --slave --args  input.gff < $(which PseudoGeneDetect.R) 

GFFmerge.bash   Pseudogene.gff  CandiFunctional.gff   Abnormal

GFFmerge.bash   Abnormal_merge.gff  NoHit.gff  Abnormal_NoHit

GFFmerge.bash   Abnormal_NoHit_merge.gff Normal.gff   PseudoEval

cp PseudoEval_merge.gff ../${Name}_PseudoEval_merge.gff

cd ..

#EOF