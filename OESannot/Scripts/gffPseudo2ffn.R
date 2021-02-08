#
# gffPseudo2ffn.R
#


chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)
if(!require(stringr)) install.packages("stringr")
library(stringr)


args       <- commandArgs(trailingOnly = T)
input.gff  <- args[1]
input.fna  <- args[2]
N_pos      <- as.numeric(args[3]) # colum position of gff, to be the name of sequences


fna <- read.fasta(input.fna)
gff <- as.matrix(read.table(input.gff))
#cds <- subset(gff, grepl("gff", gff[,3]))

file.create("gffPseudo2ffn.ffn")

for(i in 1:nrow(gff)){
#  Num   <- formatC(i, width=4, flag="0")
#  LocusTag <- paste(N_tag, Num, sep="_")
  Seqname  <- names(fna)
  Contig   <- gff[i,1]
  TContig  <- match(Contig, Seqname) # to teat multiple Contigs

  Method   <- gff[i,2]
  ORFtype  <- gff[i,3]
  Start    <- as.numeric(gff[i,4])
  End      <- as.numeric(gff[i,5])
  Score    <- as.character(gff[i,6])     # must be a character, to treat "short" status ("Short|xx")
  Strand   <- gff[i,7]
  Frame    <- gff[i,8]
  Note     <- gff[i,9]
  Name     <- gsub(" ", "", paste(Contig, Method, ORFtype, Start, End, Score, Strand, Frame, Note, sep=":" ))

  if ( gff[i,3] == "pseudo"  ){   
    Pseudo   <- fna[[TContig]][gff[i,4]:gff[i,5]]
    pName    <- gff[i,N_pos]
    write.fasta(toupper(Pseudo), pName, file.out="gffPseudo2ffn.ffn", open="a")
  }  

}
