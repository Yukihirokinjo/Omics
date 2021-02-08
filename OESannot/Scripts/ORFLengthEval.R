##
## ORFLengthEval.R 
##

chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)
if(!require(stringr)) install.packages("stringr")
library(stringr)


args       <- commandArgs(trailingOnly = T)
bla.out    <- args[1]
input.faa  <- args[2]
cogids     <- args[3]
coglen     <- args[4]
ref.faa    <- args[5]
Prop       <- as.numeric(args[6])

file.create("tmp_blaCOG.out", "tmp_ShortEval.gff", "BlastHitLengthCheck.txt")


## covert protein-id into COGid
bla    <- read.table(bla.out, head=F)
qryseq <- read.fasta(input.faa)
refseq <- read.fasta(ref.faa)
COGinfo <- as.matrix(read.table(cogids, sep=","))
COGLen  <- as.matrix(read.table(coglen))

Refs    <- bla[,2]
Pids    <- sapply(strsplit(as.character(Refs), "\\|"), "[", 2)  # extract Pids from referencce sequence in cogxxx-xxx.fasta
# convert Pid (of reference sequence in cog) to COG-id.
for(i in 1:length(Pids)){
  COGid   <- COGinfo[match(Pids[i], COGinfo[,3]),7] 
  if( is.na(COGid) ){                               # because of cogDB was used as Blast reference, all hit have a COGid, but just in case
    blaCOG  <- paste(bla[i,1], bla[i,2], bla[i,11], sep="\t")
    write(blaCOG, ncolumns = 3, file="tmp_blaCOG.out", append=T)
  }else{
    blaCOG  <- paste(bla[i,1], COGid, bla[i,11], sep="\t")
    write(blaCOG, ncolumns = 3, file="tmp_blaCOG.out", append=T)
  }
}

#blaCOG   <- as.vector(unique(bla[,1]))                          # collapse alignment redundancy
blaCOG  <- as.matrix(read.table("tmp_blaCOG.out"))

Head   <- c("Qry_Name", "Qry_Len", "Mean_HitLen",  "Acc.Len", "Status", "Propotion")
write(Head, file="BlastHitLengthCheck.txt", ncolumns=6, append=F)

for (i in 1:nrow(blaCOG)){
  HitCOGidLenInfo  <- subset(COGLen, COGLen[,1] == blaCOG[i,2])     # output of subset is a list
  if(length(HitCOGidLenInfo) == 0){
    RefLen     <- length(subset(refseq, names(refseq) == blaCOG[i,2])[[1]])
    HitCOGidLens_M   <- RefLen
  }else{
    HitCOGidLens_M   <- as.numeric(HitCOGidLenInfo[5])
  }

  if(is.na(HitCOGidLenInfo[3])){
    HitCOGidLens_SD  <- 0
  }else{
    HitCOGidLens_SD  <- as.numeric(HitCOGidLenInfo[3])
  }

  AcLen      <- HitCOGidLens_M
  QryLen     <- length(subset(qryseq, names(qryseq) == blaCOG[i,1])[[1]])
  QLenProp   <- format( (QryLen / HitCOGidLens_M), digits = 2)
  QryName    <- subset(names(qryseq), names(qryseq) == blaCOG[i,1])

  Evalue     <- as.numeric(blaCOG[i,3])

  splitName  <- strsplit(blaCOG[i,1], ":")[[1]]
  Score      <- as.numeric(strsplit(splitName[6], "\\|")[[1]])

  if(QryLen > (HitCOGidLens_M * Prop)){
    Status <- paste(Score[1], "Normal", QLenProp, Evalue, sep="|")
  } else {
    Status <- paste(Score[1], "Short", QLenProp, Evalue, sep="|")
  }

  info       <- paste(QryName, QryLen, HitCOGidLens_M, AcLen, Status, QLenProp)
  write(info, file="BlastHitLengthCheck.txt", ncolumns = 5, append=T)

  outS        <- paste(splitName[1], splitName[2], splitName[3], splitName[4], splitName[5], Status, splitName[7], splitName[8], splitName[9], sep="\t")
  write(outS, file="tmp_ShortEval.gff", ncolumns = 9, append=T)

}

##END OF FILE
