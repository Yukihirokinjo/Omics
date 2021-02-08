#
# PPgeneDetect.R
#
# Implemented in PPgeneDetect.sh

chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)
if(!require(stringr)) install.packages("stringr")
library(stringr)


args       <- commandArgs(trailingOnly = T)
blaMHC.out <- args[1]
input.fna  <- args[2]
gCode      <- as.numeric(args[3])


## ---------------------------------------Extract candidate pseudo-pseudogenes
# check1: same strand
# check2: diff. frame

blaMHC       <- as.matrix(read.table(blaMHC.out,head=F)) # Qid <\t> Refid <\t> NumHits

file.create("tmp_PP.txt", "tmp_TP.txt", "PPgene.gff", "PseudoGeneCandi.gff", "PPrecode.faa", "PPrecode.fna") 


Qinfo  <- matrix("NA", nc=9, nr=nrow(blaMHC))
for(n in 1:nrow(blaMHC)){
  Qinfo[n,] <- as.vector(strsplit(blaMHC[n,1], ":")[[1]])
}


i <- 1
while( i < nrow(blaMHC)){
  numHit  <- as.numeric(blaMHC[i,3])

  ## Eval. duplication
  # get blast-alignment region
  HitRgn      <- strsplit(blaMHC[i:(i+(numHit-1)),4], "\\|")
  HitRgnStart <- sort( as.numeric(sapply(HitRgn, "[", 1)) ) # need sort for reverce hit (i.e. 5'-Hit(i+1)-Hit(i)-3' )
  HitRgnEnd   <- sort( as.numeric(sapply(HitRgn, "[", 2)) )
  # get blast-alignment overlap
  for(k in 1:(length(HitRgn)-1)){
    HitRgnOVL <- HitRgnEnd[k] - HitRgnStart[k+1]
  }

  # length check: 
  LDstatL  <- strsplit(Qinfo[i:(i+(numHit-1)),6], "\\|")
  LDstats  <- as.numeric(sapply(LDstatL, "[", 3))
  #sumLD    <- sum(LDstats)
  QryLen   <- (as.numeric(Qinfo[i,5]) - as.numeric(Qinfo[i,4]) +1) /3
  RefLen   <- QryLen/LDstats[1]


  # same strand -> -> if there is even one fragment corded on different strand, the fragments set is assigned as pseudo 
  strand1 <- as.character(Qinfo[i,7])
  strands <- as.vector(Qinfo[i:(i+(numHit-1)),7])
  sameS   <-  length(subset(strands, strands == strand1) ) == numHit

  # whether all continuous fragments are on different frame -> if there is even one continuous same frames, the fragments set is assigned as pseudo 
  frame1  <- as.numeric(Qinfo[i,4]) %% 3  # first frame
  frames  <- as.numeric(Qinfo[i:(i+(numHit-1)),4]) %% 3 # all frames (in the fragments set)
  Fdif <- rep(1, length=(length(frames)-1))
  for(k in 1:(length(frames)-1)) { Fdif <- frames[i] - frames[i+1] }
  diffF   <- is.na(match(0, Fdif))

  if(HitRgnOVL < RefLen * 0.5){  # Fragmented ORFs should have short overlap: set as shorter than 50% of the reference protein
  if( sameS && diffF ){                                                   # Condition for pseudo-pseudo genne ORFs (on same strand and diff. frame)
    write.table(blaMHC[i:(i+(numHit-1)),], file="tmp_PP.txt", row.names=F, col.names=F,  quote=F, append=T)
  } else {                                                                # Others should be split gene or true pseudogene
    write.table(blaMHC[i:(i+(numHit-1)),], file="tmp_TP.txt", row.names=F, col.names=F,  quote=F, append=T)
    tag2 <- paste("PseudogeneCndi(MultiHitConti)", blaMHC[i,2], sep="-")
    for(n in i:i:(i+(numHit-1))){
      tpgff2 <- paste(Qinfo[n,1], "PPgeneDetect", "CDS", Qinfo[n,4], Qinfo[n,5], Qinfo[i,6], Qinfo[i,7], Qinfo[i,8], tag2, sep="\t")
      write(tpgff2, file="PseudoGeneCandi.gff", ncolumns = 9, append=T)         # gff for ORFLengthEval.sh
    }
  }
  }  
  i <- i + numHit
}


## ------------------------------------------------------------------------------



## ------------------------------------------------------------------------------

###
###--------------------------------   2 Function definition   ------ Start
###

# function1: HPTrecode 
HPTrecode <- function(ORF, DNA, Type, HPTpos){ 

  POSs   <- as.numeric(strsplit(as.character(HPTpos), ":")[[1]])
  ADLs   <- strsplit(Type, ":")[[1]]
  adls   <- rep(0, length(POSs))
  DNAs   <- strsplit(DNA, ":")[[1]]

  for(i in 1:length(POSs)){
    switch(ADLs[i],
      "ad"  =  adls[i] <-  1,
      "dl"  =  adls[i] <- -1,
    )
    HPT    <- c2s(rep(DNAs[i], 8))
    HPTadl <- c2s(rep(DNAs[i], 8 + adls[i]))

    POSmod   <- POSs[i] + sum(adls[0:(i-1)])
    HPTstart <- POSmod
    HPTend   <- POSmod + 7
    substring(ORF, HPTstart, HPTend)  <- HPT  # Replase "a{8}" to "A{8}" to distiguish from other hits (All other sequence characters in PPorf are "lower letter") 
    ORF    <- sub(HPT, HPTadl, ORF)         # insert/delete 1nt 
  }
  return(ORF)

} 


# function2: Check_tPP
Check_tPP <- function(ad, dl, Name, FileOut){

  tPP_1ad   <- translate(ad, sens="F", numcode=gCode)
  tPP_1dl   <- translate(dl, sens="F", numcode=gCode)

  # get stop codon positions
  ePP_1ad   <- gregexpr("\\*", c2s(tPP_1ad))  # need duble"\" to escape meta-charactor "*"
  names(ePP_1ad)<-"position"
  ePP_1dl   <- gregexpr("\\*", c2s(tPP_1dl))
  names(ePP_1dl)<-"position"

  outN      <- paste(FileOut, "fna", sep=".")
  outA      <- paste(FileOut, "faa", sep=".")
  name_1ad  <- paste(Name,"1ad", sep="")
  name_1dl  <- paste(Name,"1dl", sep="")
  
  # whether insertion/eletion eliminates other stop codons (so whether new ORF has "single stop codon" in the "final position") 
  logi_Single_ad <- length(ePP_1ad$position) == 1            # whether new ORF has "single stop codon"
  logi_final_ad  <- ePP_1ad$position[1] == length(tPP_1ad)   # whether the SingleStopC is in "final position"
  logi_Single_dl <- length(ePP_1dl$position) == 1
  logi_final_dl  <- ePP_1dl$position[1] == length(tPP_1dl)

  if( logi_Single_ad  &&  logi_final_ad ){
    return("ad")
  } else if ( logi_Single_dl  &&  logi_final_dl ){
    return("dl")
  } else {
    return("no")
  }
}



# function3: PP_Check

PP_Check <- function(PPorf){
  Slip  <- rep("no", 4)   # <- (add_T, del_T, add_A, del_A)
  # "T"
  FindT    <- "no"
  # get HPT pos. in each ORF
  posT        <- gregexpr("t{8}", PPorf)
  names(posT) <- "position"
  # make slippage frame for all HPTs
  if(posT$position[1] != -1){
  for( hitT in length(posT$position):1){         # Start from most-anterior hit
    PositionT  <- as.character(posT$position[hitT])
    PPorf_1ad <- HPTrecode(PPorf, "T", "ad", PositionT)
    PPorf_1dl <- HPTrecode(PPorf, "T", "dl", PositionT)
    # check PP gene completeness
    FindT <- Check_tPP( s2c(PPorf_1ad), s2c(PPorf_1dl), Name="CheckPP", FileOut="PPrecode")
    if(FindT != "no"){
      switch(FindT,
      "ad"  =   Slip[1] <- PositionT,
      "dl"  =   Slip[2] <- PositionT,
      )
      break  # break roop when PPrecode-frame was found (forT)
    }
  } # end of T roop
  } # end of if (posT$position[1] != -1)

  # "A"
  FindA   <- "no"
  posA        <- gregexpr("a{8}", PPorf)
  names(posA) <- "position"
  if(posA$position[1] != -1){
  for(hitA in length(posA$position):1){
    PositionA  <- as.character(posA$position[hitA])
    PPorf_1ad <- HPTrecode(PPorf, "A", "ad", PositionA)
    PPorf_1dl <- HPTrecode(PPorf, "A", "dl", PositionA)
    FindA     <- Check_tPP( s2c(PPorf_1ad), s2c(PPorf_1dl), Name="CheckPP", FileOut="PPrecode")
    if(FindA != "no"){
      switch(FindA,
      "ad"  =   Slip[3] <- PositionA,
      "dl"  =   Slip[4] <- PositionA,
      )
      break
    } 
  } # end of A roop
  } # end of if (posA$position[1] != -1)

  # analyse type of Slippage (A/T, position, ad/dl)
  Slip2    <- as.numeric(Slip)
  Slippage <-  subset(Slip2, !is.na(Slip2))
  if( length(Slippage) != 0 ){
    PPmax  <- max(Slippage)
    Type   <- grep(PPmax, Slip2)
    switch(Type,
      "1"  =  PPout <- paste("ad", "T", PPmax, sep=":"),
      "2"  =  PPout <- paste("dl", "T", PPmax, sep=":"),
      "3"  =  PPout <- paste("ad", "A", PPmax, sep=":"),
      "4"  =  PPout <- paste("dl", "A", PPmax, sep=":"),
    )
  }else{
    PPout  <- 0
  }
  return(PPout) 
} # end of PP_Check

###
###--------------------------------  2 Function definition   ------ End
###





###
###--------------------------------  3 search PP gene (main procedure)-------Start
###
# FrameShift ORF by insert or delete all detected HPTs; and if ..

if(file.info("tmp_PP.txt")$size != 0 ){  # if any PP candi exists

  fna <- read.fasta(input.fna)
  pp  <- as.matrix(read.table("tmp_PP.txt", head=F))
  Ref <- unique(pp[,2])

  # roop for each set of MultiHit ORFs
  for(k in 1:length(Ref)){
    pp1      <- subset(pp, pp[,2] == Ref[k], select=V1)  # pp1 <- matrix
    pp2      <- strsplit(pp1, ":")                       # pp2 <- list
    NumORF   <- length(pp2)
    NumSplit <- NumORF -1

    PPget    <- rep(0, NumSplit)
    HPTpos   <- rep(0, NumSplit)
    SlipType <- rep(0, NumSplit)

    seqID    <- pp2[[1]][1]
    Contig   <- match(seqID, names(fna))
    Strand   <- pp2[[1]][7]

    StartA   <- as.numeric(pp2[[1]][4])
    EndA     <- as.numeric(pp2[[NumORF]][5])

    # modify sequence and positions to ignore strand information
    Start  <- rep(0, length(pp2))
    End    <- rep(0, length(pp2))
    for(n in 1:length(pp2)){
      if(Strand == "+"){
        Fna      <- fna[[Contig]]
        Start[n] <- as.numeric(pp2[[n]][4])
        End[n]   <- as.numeric(pp2[[n]][5])
      }else if(Strand == "-"){                                  ## If reverse strand,,,
        Fna      <- rev(comp(fna[[Contig]]))                    # change sequence into RevComp
        Start[n] <- length(Fna) - as.numeric(pp2[[n]][5]) + 1   # change positions to fit changed sequence
        End[n]   <- length(Fna) - as.numeric(pp2[[n]][4]) + 1
      }
    }
    StartS <- sort(Start)
    EndS   <- sort(End)

    PPfit <- "yes"
    for( i in 1:NumSplit ){
      ## prepare row sequence (combined)
      PPorf <- c2s(Fna[StartS[i]:EndS[i+1]])
      PPget[i]  <- PP_Check(PPorf)                           # PPget: "type:DNA:Pos"
      if(PPget[i] == "0" ){
        PPfit <- "no"
      }

    } # end of for( i in 1:(length(pp2)-1) ) : PPfit

    if(PPfit != "no"){
        HPTpos <- rep(0, NumSplit)
        Type   <- rep(0, NumSplit)
        DNA    <- rep(0, NumSplit)
        for( h in 1:NumSplit ){
          Type[h]   <- strsplit(PPget[h], ":")[[1]][1]
          DNA[h]    <- strsplit(PPget[h], ":")[[1]][2]
          HPTpos[h] <- as.numeric(strsplit(PPget[h], ":")[[1]][3])
          HPTpos[h] <- HPTpos[h] + (StartS[h] - StartS[1])
          ORF       <- c2s(Fna[StartS[1]:EndS[NumORF]])
          
          HPTs   <- paste(HPTpos, collapse=":")
          Types  <- paste(Type, collapse=":")
          DNAs   <- paste(DNA, collapse=":")
        }
        ORFrecod <-  HPTrecode(ORF, DNAs, Types, HPTs)

    # write candidate ORFs (Slipped frames)
      PPstat   <- paste(DNAs, Types, HPTs, sep="|")
      orfname  <- paste(seqID,  pp2[[1]][2], pp2[[1]][3], pp2[[1]][4], pp2[[NumORF]][5], "NA", Strand, PPstat, "Pseudo-pseudo", sep=":")
      write.fasta(s2c(ORFrecod), orfname, file.out="PPrecode.fna", open="a")
      write.fasta(translate(s2c(ORFrecod), numcode=gCode), orfname, file.out="PPrecode.faa", open="a")
      PPgff    <- paste(seqID, "PPgeneDetect", "PPgene", StartA, EndA, "NA", Strand, PPstat, "Pseudo-pseudo", sep="\t")
      write(PPgff , file="PPgene.gff", ncolumns=9, append=T)
    }else{
      pNote <- paste("PseudogeneCndi(MultiHitConti)", Ref[k], sep="-")
      for(x in 1:length(pp2)){
        TPgff  <- paste(seqID, "PPgeneDetect", "CDS", pp2[[x]][4], pp2[[x]][5], pp2[[x]][6], Strand, "0", pNote, sep="\t")
        write(TPgff, file="PseudoGeneCandi.gff", ncolumns=9, append=T)
      }
    } # end of PPfit eval
  } # end of for(k in 1:length(Ref))


}else{
  NoPP <- "#No PPgene was detected"
  write(NoPP , file="PPgene.gff", ncolumns=1, append=F)
} # end of tmp_PP.txt roop







###
###--------------  PPAP exhaustive search mode  --------------
###




###
###--------------  PPAP exhaustive search mode  --------------
###

