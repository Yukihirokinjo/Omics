#
# BlastMultiHitExtract.R
#
# Implemented in PPgeneDetect.sh


chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)


args       <- commandArgs(trailingOnly = T)
bla.out    <- args[1]
cogids     <- args[2]

###
###--------------  PPAP blast dual-hit based mode  --------------
###

## prepare output files
file.create("tmp_blaCOG.out", "tmp_BlastMultiHit_Continuous.txt")


## covert protein-id into COGid
bla     <- as.matrix(read.table(bla.out,head=F))  # blast query(.faa) moust be prepared by "gff2faa.sh" ; need its unique id structure
COGinfo <- as.matrix(read.table(cogids, sep=","))
Refs    <- bla[,2]
Pids    <- sapply(strsplit(as.character(Refs), "\\|"), "[", 2)
for(i in 1:length(Pids)){
  COGid   <- COGinfo[match(Pids[i], COGinfo[,3]),7]
  if( is.na(COGid) ){
    blaCOG  <- paste(bla[i,1], bla[i,2], bla[i,3], bla[i,4], bla[i,5], sep="\t")
    write(blaCOG, ncolumns = 5, file="tmp_blaCOG.out", append=T)
  }else{
    blaCOG  <- paste(bla[i,1], COGid, bla[i,3],  bla[i,4], bla[i,5], sep="\t")
    write(blaCOG, ncolumns = 5, file="tmp_blaCOG.out", append=T)
  }
}



## extract continuous ORFs which hit to same COGid -> # continuous in order of blasthit list -> OK!(ref[i]==ref[i+1])

blaC <- as.matrix(read.table("tmp_blaCOG.out"))

preHit <- 0
HitRgn <- 0 
for(i in 1:nrow(blaC)){
  if(i == nrow(blaC)){ # for final record
    if(preHit >=1){
      HitRgn <- c(HitRgn, gsub(" ", "", paste(blaC[i,4], blaC[i,5], sep="|")))
      write.table(cbind( blaC[(i -preHit):i,1:2], rep(preHit+1, preHit+1) , HitRgn), file="tmp_BlastMultiHit_Continuous.txt", row.names=F, col.names=F,  quote=F, append=T )
    }
  }else{ # untill the final record
    if(blaC[i,2] == blaC[i+1,2]){
      if(preHit == 0){
        HitRgn <- gsub(" ", "", paste(blaC[i,4], blaC[i,5], sep="|"))
      }else if(preHit >= 1){
        HitRgn <- c(HitRgn, gsub(" ", "", paste(blaC[i,4], blaC[i,5], sep="|")))
      }
      preHit <- preHit + 1
    }else{
      if(preHit >= 1){ # final record of each continuous-hit
        HitRgn <- c(HitRgn, gsub(" ", "", paste(blaC[i,4], blaC[i,5], sep="|")))
        write.table(cbind( blaC[(i -preHit):i,1:2], rep(preHit+1, preHit+1), HitRgn), file="tmp_BlastMultiHit_Continuous.txt", row.names=F, col.names=F,  quote=F, append=T )
      }
      preHit <- 0
    }
  }
}


#rownames(blaC) <- c(1:nrow(blaC))
#dupRef <- unique(blaC[duplicated(blaC[,2]),2])

#for(i in 1:length(dupRef)){
#  Qids  <- subset(blaC, blaC[,2] == dupRef[i], select = V1)
#  rQids <- as.numeric(rownames(Qids))
#  if( (rQids[length(rQids)] - rQids[1]) == (length(rQids) - 1) ){                #mod
#    for(k in 1:length(Qids)){
#      write(c(Qids[k], dupRef[i], length(Qids)), file="tmp_BlastMultiHit_Continuous.txt", ncolumns = 3, append=T )
##      QidSep   <- strsplit(Qids[k], ":")[[1]] # change to vector
##      LDStat     <- QidSep[6]
##      gffMulti <- paste(QidSep[1], "BlastMultiHitConti", "CDS", QidSep[4], QidSep[5], LDStat, QidSep[6], "0", dupRef[i] )  # Score colum indicates number of hits
##      write(gffMulti, file="tmp_BlastMultiHit_Continuous.gff", ncolumns = 9, append=T)
#    }
#  }else{
#    for(k in 1:length(Qids)){
#      write(c(Qids[k], dupRef[i], length(Qids)), file="tmp_BlastMultiHit_nonContinuous.txt", ncolumns = 3, append=T)
#    }
#  } 
#}

## extract ORFs which hit to same protein (<300nt apart, <3,000nt apart): for large insertion (e.g. mobile element)
# need to distinguish from gene duplication <- check ORF length

