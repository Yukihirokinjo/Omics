#
# GetDBseqLength.R
#

chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)



args <- commandArgs(trailingOnly = T)
Ref.faa  <- args[1]
FileOut  <- args[2]


Head  <- c("#RefID", "SeqLength")
write(Head, file=FileOut, ncolum=2, append=F)

refDB <- read.fasta(Ref.faa)
for( i in 1:length(refDB) ){
  Name    <- names(refDB[i])
  Length  <- length(refDB[[i]])
  NLinfo  <- paste(Name, Length, sep="\t") 
  write(NLinfo, file=FileOut, ncolum=2, append=T)
}

