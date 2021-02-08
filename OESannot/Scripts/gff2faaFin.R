#
# gff2faaFin.R
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
Gcode      <- as.numeric(args[4])
InfoTag    <- as.numeric(args[5])

print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])

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
preGff <- as.matrix(read.table(input.gff))
#cds <- subset(gff, grepl("gff", gff[,3]))


FileOutProt <- paste(N_tag, "faa", sep=".")
FileOutCDS  <- paste(N_tag, "ffn", sep=".")
FileOutRNA  <- paste(N_tag, "frn", sep=".")

file.create(FileOutProt, FileOutCDS, FileOutRNA)

Contigs <- unique(preGff[,1])
for(c in 1:length(Contigs)){
  gff <- subset(preGff, preGff[,1] == Contigs[c])

  for(i in 1:nrow(gff)){
    Contig   <- gff[i,1] 
    Num      <- formatC((i * 10), width=5, flag="0")
    Seq_ID   <- paste(as.roman(c), Num, sep="")
    LocusTag <- paste(N_tag, Seq_ID, sep="_")
    Seqname  <- names(fna)
    Strand   <- gff[i,7]
    TContig  <- match(Contig, Seqname) # to teat multiple Contigs
    InfoType     <- gff[i,InfoTag]

    Name     <- gsub(" ", "", paste(InfoType, LocusTag, sep=":" ))

    if(gff[i,3] == "CDS"){
      switch(Strand,
        "+"  =   ORF    <- fna[[TContig]][gff[i,4]:gff[i,5]],
        "-"  =   ORF    <- rev(comp(fna[[TContig]][gff[i,4]:gff[i,5]])),
      )
      write.fasta(toupper(ORF), Name, file.out=FileOutCDS,  open="a")
      switch(Strand,
        "+"  =   tORF    <- translate(fna[[TContig]][gff[i,4]:gff[i,5]], sens="F", frame=0, numcode = Gcode),
        "-"  =   tORF    <- translate(rev(comp(fna[[TContig]][gff[i,4]:gff[i,5]])), sens="F", frame=0, numcode = Gcode),
      )
      write.fasta(tORF, Name, file.out=FileOutProt,  open="a")
    } else if ( str_detect(gff[i,3], "RNA") ){                     # for RNA
      Rna   <-  fna[[TContig]][gff[i,4]:gff[i,5]]
    
      write.fasta(toupper(Rna), Name, file.out=FileOutRNA,  open="a")
    } else if ( gff[i,3] == "PPgene"  ){                     # for RNA
      switch(Strand,
        "+" =  pORF    <- c2s(fna[[TContig]][gff[i,4]:gff[i,5]]),
        "-" =  pORF    <- c2s(rev(comp(fna[[TContig]][gff[i,4]:gff[i,5]]))),
      )
      DNAs  <- strsplit(strsplit(gff[i,8], "\\|")[[1]][1], ":")[[1]]  # colum8 of PPgene is frameshift information 
      Types <- strsplit(strsplit(gff[i,8], "\\|")[[1]][2], ":")[[1]]
      HPTs  <- strsplit(strsplit(gff[i,8], "\\|")[[1]][3], ":")[[1]]
      ORFrecod <-  HPTrecode(pORF, DNAs, Types, HPTs)
      write.fasta(toupper(s2c(ORFrecod)), Name, file.out=FileOutCDS,  open="a")
      write.fasta(translate(s2c(ORFrecod), numcode = Gcode), Name, file.out=FileOutProt,  open="a")
    }

  # write new gff : include locus_tag
    tagNote <- paste(LocusTag, gff[i,9], sep=";")
    tagGFF  <- paste(gff[i,1], gff[i,2], gff[i,3], gff[i,4], gff[i,5], gff[i,6], gff[i,7], gff[i,8], tagNote, sep="\t")
    write(tagGFF, file="addLocusTag.gff", sep="\t", ncolumns=9, append=T)

  }

}

# END OF FILE
