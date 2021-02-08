##
## DomCompEval.R 
##
## USAGE: DomCompEval.R [ RPSblast_output ] [ input.faa ] [ cddid_extract.tbl ]
##
##
## Strategy: 
## 1. split hits into each independent region (no need; rep mode can automatically do that)
## 2. select longest match in each region (no need; rep mode can automatically do that)
## 3. select lowest E-value  (no need; rep mode can automatically do that)
## 4. evaluate PSSMs length completeness  (if longest PSSM has insufficient length but some shorters have sufficient/complete length?)
## 5. evaluate whether insufficient length is due to CDS fragmentation or not

chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("IRanges")
library(IRanges)

args        <- commandArgs(trailingOnly = T)
RPSbla.out  <- args[1]
input.faa   <- args[2]
CDDid.tbl   <- args[3]
Prop        <- as.numeric(args[4])

# 0. Add length information into RPSblast output

RPS  <- as.matrix(read.table(RPSbla.out))
faa  <- read.fasta(input.faa)
tbl  <- read.table(CDDid.tbl)

Head <- c("#QryID", "QryLen", "CDDid", "CDDLen", "Ident", "AlignLen", "MisMatch", "GapOpen", "QhitFrom", "QhitTo", "CDhitFrom", "CDhitTo", "Evalue", "BitScore")
write(Head, file="RPS_custom6mod", ncolumns=14, sep="\t", append=F)

for(i in 1:nrow(RPS)){
  CDDid   <- strsplit(RPS[i,2], ":")[[1]][2]  # to remove "CDD:" from rpsblast output
  Qryid   <- RPS[i,1]
  CDDLen  <- as.numeric(subset(tbl, tbl$V1 == CDDid)[4])
  Seq     <- subset(faa, names(faa) == Qryid)                          # extract sequence which match to the Qryid
  QryLen  <- length(Seq[[1]])
  RPSmod  <- c(Qryid, QryLen, CDDid, CDDLen, RPS[i,3:12])
  write(RPSmod, file="RPS_custom6mod", ncolumns=14, sep="\t", append=T)
}

RPS2  <- as.matrix(read.table("RPS_custom6mod"),head=T)

# 1. Remove redundancy  (remove overlapped and low Evalue)


write(Head, file="RPS_custom6rep", ncolumns=14, sep="\t", append=F)

uniID  <- unique(RPS2[,1])

for(i in 1:length(uniID)){
  DupHit    <- subset(RPS2, RPS2[,1] == uniID[i])
  Ranges    <- IRanges(start=c(as.numeric(DupHit[,9])),end=c(as.numeric(DupHit[,10])))
  rRanges   <- reduce(Ranges)

  # Extract tophit in each range
  for(k in 1:length(rRanges)){
    oRanges   <- findOverlaps(Ranges, rRanges[k])
    RangeBest <- as.matrix(oRanges)[1,1]    # extract besthit in each range blast output is arleady sorted by Evalue
    write(DupHit[RangeBest,], file="RPS_custom6rep", ncolumns=14, sep="\t", append=T)  
  }
}


# 2. Evaluate Domain completeness

RPSrep     <- read.table("RPS_custom6rep", head=T)

Head2 <- c("#QryID", "QryLen", "CDDid", "CDDLen", "Ident", "AlignLen", "MisMatch", "GapOpen", "QhitFrom", "QhitTo", "CDhitFrom", "CDhitTo", "Evalue", "BitScore", "CompStat", "Partial")
write(Head2, file="RPS_custom6fin", ncolumns=16, sep="\t", append=F)

for(i in 1:length(RPSrep[[1]])){
  DomStart   <- RPSrep[i,11]
  DomEnd     <- RPSrep[i,12]
  HitDomLen  <- (DomEnd - DomStart)
  OrgDomLen  <- RPSrep[i,4]

  AlnStart   <- RPSrep[i,9]
  AlnEnd     <- RPSrep[i,10]
  OrgSeqLen  <- RPSrep[i,2]

  Partial <- "-"

  if( HitDomLen > ( OrgDomLen * Prop ) ){
    CompStat <- "-"
    Partial <- "-"
  }else{
    if( (DomStart > (OrgDomLen * 0.2)) & ((DomEnd / OrgDomLen) < Prop) ){
      CompStat <- "NC"

      ## 3. evaluate whether incomplete domain is due to incomplete ORF
      if((DomStart * Prop) > AlnStart){ Partial <- "yes" # due to 5' deletion
      }else if(((OrgDomLen - DomEnd) * Prop) > (OrgSeqLen - AlnEnd)){ Partial <- "yes" #  due to 3' deletion
      }else{ Partial <- "no" } # not due to incomplete ORF
    }else if( DomStart > (OrgDomLen * 0.2) ){
      CompStat <- "N"
      if((DomStart * Prop) > AlnStart){ Partial <- "yes" }else{ Partial <- "no" }
    }else if( (DomEnd / OrgDomLen) < Prop ){
      CompStat <- "C"
      if(((OrgDomLen - DomEnd) * Prop) > (OrgSeqLen - AlnEnd)){ Partial <- "yes" }else{ Partial <- "no" }
    }
  }

  # "Partial = yes": partial domain hit due to incomplete ORF (5'/3' deletion)
  # "Partial = no" : partial domain hit but ORF is likely complete
  # "Partial = -"  : complete domain hit (CompStat is also "-")
  RPSfin <- c(as.matrix(RPSrep[i,1:14]), CompStat, Partial)
  write(RPSfin, file="RPS_custom6fin", ncolumns=16, sep="\t", append=T)
}

## write gff file for Partial ORF

Fin <- as.matrix(read.table("RPS_custom6fin"))
file.create("DomEval.gff")
Faa <- names(faa)

for(i in 1: length(Faa)){
  name   <- strsplit(Faa[i], ":")[[1]]
  DomHit <- subset(Fin, Fin[,1]==Faa[i])  # matrix
  if( length(DomHit) == 0 ){
    Stat <- paste(name[6], "NoDomHit",    sep="|")
  }else{
    numHit  <- nrow(DomHit)
    DomStatEdge <- c(DomHit[1,16], DomHit[numHit,16])
    DomStatTotal <- as.character(length(grep("yes", DomStatEdge)))
    Stat=switch(DomStatTotal,
      "2"  =  paste(name[6], "DomPartial",  sep="|"),
      "1"  =  paste(name[6], "DomPartial",  sep="|"),
      "0"  =  paste(name[6], "DomHit", sep="|"),
    )
  }
  gff  <- paste(name[1], name[2], name[3], name[4], name[5], Stat, name[7], name[8], name[9] , sep="\t")
  write(gff, file="DomEval.gff", ncolumns=9, append=T)
}

##END OF FILE
