#
# gff2faa.R
#

chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)
if(!require(stringr)) install.packages("stringr")
library(stringr)


args       <- commandArgs(trailingOnly = T)
input.gff  <- args[1]
input.fna  <- args[2]
N_tag      <- args[3]
gCode      <- as.numeric(args[4])


###################################### function1: HPTrecode
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
######################################


fna <- read.fasta(input.fna)
gff <- as.matrix(read.table(input.gff))
#cds <- subset(gff, grepl("gff", gff[,3]))

file.create("gff2.faa", "gff2RNA.fna")

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

  if(gff[i,3] == "CDS"){
    switch(Strand,
      "+"  =   tORF    <- translate(fna[[TContig]][gff[i,4]:gff[i,5]], sens="F", frame=0, numcode=gCode),
      "-"  =   tORF    <- translate(rev(comp(fna[[TContig]][gff[i,4]:gff[i,5]])), sens="F", frame=0, numcode=gCode),
    )
      write.fasta(tORF, Name, file.out="gff2.faa", open="a")
  } else if ( str_detect(gff[i,3], "RNA") ){                     # for RNA
    Rna   <-  fna[[TContig]][gff[i,4]:gff[i,5]]
    write.fasta(toupper(Rna), Name, file.out="gff2RNA.fna", open="a")
  } else if ( gff[i,3] == "PPgene"  ){                     # for RNA
    switch(Strand,
      "+" =  ORF    <- c2s(fna[[TContig]][gff[i,4]:gff[i,5]]),
      "-" =  ORF    <- c2s(rev(comp(fna[[TContig]][gff[i,4]:gff[i,5]]))),
    )
    DNAs  <- strsplit(strsplit(gff[i,8], "\\|")[[1]][1], ":")[[1]]  # colum8 of PPgene is frameshift information 
    Types <- strsplit(strsplit(gff[i,8], "\\|")[[1]][2], ":")[[1]]
    HPTs  <- strsplit(strsplit(gff[i,8], "\\|")[[1]][3], ":")[[1]]
    ORFrecod <-  HPTrecode(ORF, DNAs, Types, HPTs)
    write.fasta(translate(s2c(ORFrecod)), Name, file.out="gff2.faa", open="a")
  } else if ( gff[i,3] == "pseudo"  ){   
    Pseudo   <- fna[[TContig]][gff[i,4]:gff[i,5]]
    write.fasta(toupper(Pseudo), Name, file.out="gff2Pseudo.fna", open="a")
  }  

}
