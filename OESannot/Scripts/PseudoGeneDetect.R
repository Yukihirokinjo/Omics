##
## PseudoGeneDetect.sh 
##
## USAGE: 

chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)
if(!require(stringr)) install.packages("stringr")
library(stringr)

file.create("Pseudogene.gff", "CandiFunctional.gff", "Normal.gff", "NoHit.gff", "NoHitRemove.gff")

args       <- commandArgs(trailingOnly = T)
Ingff  <- args[1]

gff        <- as.matrix(read.table(Ingff))
gffMHC     <- subset(gff, str_detect(gff[,9], "(MultiHitConti)"))    # MultiHit-Continuous
gffSH      <- subset(gff, ! str_detect(gff[,9], "(MultiHitConti)"))  # SingleHit
gffSHS     <- subset(gffSH, str_detect(gffSH[,6], "Short") )   # Short_Hit
gffOther   <- subset(gffSH, ! str_detect(gffSH[,6], "Short") ) # Normal or NoHit
gffNormal1 <- subset(gffOther, str_detect(gffOther[,6], "Normal") )     # Normal1: Normal_Hit
gffOther2  <- subset(gffOther, ! str_detect(gffOther[,6], "Normal") )   # NoHit
gffNoHit   <- subset(gffOther2, str_detect(gffOther2[,6], "NoDomHit") ) # NoHit_and_NoDomHit
gffNormal2 <- subset(gffOther2, ! str_detect(gffOther2[,6], "NoDomHit") ) # Normal2: RNAs + NoHit_but_DomHit

# write normal
if(nrow(gffNormal1) > 0 || nrow(gffNormal2) > 0 ){
  gffNormal  <- rbind(gffNormal1, gffNormal2)
  write.table(gffNormal, file="Normal.gff", quote=F, col.names=F, row.names=F, sep="\t", append=T)
}else{
  write("#No Normal Genes", file="Normal.gff",  ncolumns=9, sep="\t", append=T)
}
# no hit (for both COG and Domain)

if(nrow(gffNoHit) > 0 ){
  for(n in 1:nrow(gffNoHit)){
    orgScore  <- as.numeric(strsplit(gffNoHit[n,6], "\\|")[[1]][1])
    if(orgScore > 10){
      write(gffNoHit[n,], file="NoHit.gff", ncolumns=9, sep="\t", append=T)
    }else{
      write(gffNoHit[n,], file="NoHitRemove.gff", ncolumns=9, sep="\t", append=T)
    }
  }
}else{
  write("#No NoHit Genes", file="NoHit.gff",  ncolumns=9, sep="\t", append=T)
  write("#No NoHit Genes", file="NoHitRemove.gff",  ncolumns=9, sep="\t", append=T)
}

# Single Hit
if(nrow(gffSHS) > 0){
  for(i in 1:nrow(gffSHS)){                        #      [1]           [2]                [3]                           [4]      [5]
    LDScore  <- strsplit(gffSHS[i,6], "\\|")[[1]]  # OriginalScore|Normal/Short|RelativeLength(to COGstandard-minimum)|Evalue|DomainStatus
    LenComp  <- as.numeric(LDScore[3])
    DomComp  <- LDScore[5]
    HitEval  <- as.numeric(LDScore[4])
    DomHitP  <- charmatch("DomPartial", DomComp)
    DomHitC  <- charmatch("DomHit", DomComp)

    if(!is.na(DomHitP)){
      Pseudo <- paste(gffSHS[i,1], gffSHS[i,2], "pseudo", gffSHS[i,4], gffSHS[i,5], gffSHS[i,6], gffSHS[i,7], gffSHS[i,8], "Truncated_DomPartial", sep="\t")
      write(Pseudo, file="Pseudogene.gff", ncolumns=9, sep="\t", append=T)
    }else if( !is.na(DomHitC) || LenComp > 0.6){ # Complete Domain or >60% RelativeLength of COG-hit(standard-minimum)
      write(gffSHS[i,], file="CandiFunctional.gff", sep="\t", ncolumns=9, append=T)
    }else{
      if(HitEval < 1e-5 ){ # No Domain-hit but significant(<1e-5) COG-hit
        PseudoND <- paste(gffSHS[i,1], gffSHS[i,2], "pseudo", gffSHS[i,4], gffSHS[i,5], gffSHS[i,6], gffSHS[i,7], gffSHS[i,8], "Truncated_NoDomHit", sep="\t")
        write(PseudoND, file="Pseudogene.gff", ncolumns=9, sep="\t", append=T)
      }else{  # No Domain-hit and non-significant(>1e-5) COG-hit
        PseudoNDL <- paste(gffSHS[i,1], gffSHS[i,2], "pseudo", gffSHS[i,4], gffSHS[i,5], gffSHS[i,6], gffSHS[i,7], gffSHS[i,8], "Candi_MisHit", sep="\t")
        write(PseudoNDL, file="NoHitRemove.gff", ncolumns=9, sep="\t", append=T)  
      }
    }
  } # other ORFs will be removed
}else{
  write("#No SHS Genes", file="Pseudogene.gff",  ncolumns=9, sep="\t", append=T)
  write("#No SHS Genes", file="CandiFunctional.gff",  ncolumns=9, sep="\t", append=T)
}

# MultiHit (Continuous) <- no candidata duplication
# <1stORF/2ndORF>
# Normal/Normal: DomHit/DomHit:Normal/Normal;, DomHit/DomPartial:Normal/remove, DomPartial/DomPartial:<pseudo> 
# Normal/Short:
# Short/Short:
if(nrow(gffMHC) > 0){
  NoteMHC  <- sapply(strsplit(gffMHC[,9], ":"), "[" , 1)
  refMHC   <- unique(NoteMHC)
  for(i in 1:length(refMHC)){
    MHCs     <- subset(gffMHC, gffMHC[,9] == refMHC[i])
    numHit   <- nrow(MHCs)
    MHCScore <- strsplit(MHCs[,6], "\\|")
    LenComp  <- as.numeric(sapply(MHCScore, "[", 3))
    DomComp  <- sapply(MHCScore, "[", 5)
    DomHitC  <- charmatch("DomHit", DomComp)
    if(is.na(DomHitC)){
      pStart  <- MHCs[1,4]; pEnd <- MHCs[numHit,5]
      Pseudo <- paste(MHCs[1,1], MHCs[1,2], "pseudo", pStart, pEnd, MHCs[1,6], MHCs[1,7], MHCs[1,8], "Fragmented", sep="\t")
      write(Pseudo, file="Pseudogene.gff", ncolumns=9, sep="\t", append=T)
    }else{
      for(k in DomHitC ){
        write(MHCs[k,], file="CandiFunctional.gff", ncolumns=9, sep="\t", append=T)
      }
    }
  }
}else{
  write("#No MHC Genes", file="Pseudogene.gff",  ncolumns=9, sep="\t", append=T)
  write("#No MHC Genes", file="CandiFunctional.gff",  ncolumns=9, sep="\t", append=T)
}
      #LDStatPP <- sapply(pp2, "[", [7])

## ---------------------------------------Evaluate extracted Pseudogene candidates
# -> regarding their length and/or Domain completeness 
# Score colum:   NA,  Short:0.x:NoDomHit,  Short:0.x:DomPartial,  Short:0.x:DomComplete

##END OF FILE
