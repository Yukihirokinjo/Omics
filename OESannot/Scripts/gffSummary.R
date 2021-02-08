#
# gffSummary.R
#

chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)


args       <- commandArgs(trailingOnly = T)
Ingff  <- args[1]
fna    <- args[2]
Tag    <- args[3]

FileOut <- paste(Tag, "Summary.txt", sep="_" )
file.create(FileOut)

gff <- as.matrix(read.table(Ingff))

Genes    <- 0
CDSs     <- 0
tRNAs    <- 0
rRNAs    <- 0
ncRNAs   <- 0
Pseudos  <- 0
PPgene  <- 0

for( i in 1:nrow(gff)){
  switch(gff[i,3],
    "CDS"     =  CDSs    <- CDSs + 1,
    "PPgene"  =  {PPgene <- PPgene + 1  ; CDSs <- CDSs + 1 },
    "tRNA"    =  tRNAs   <- tRNAs + 1,
    "rRNA"    =  rRNAs   <- rRNAs + 1,
    "ncRNA"   =  ncRNAs <- ncRNAs +1,
    "tmRNA"   =  ncRNAs <- ncRNAs +1,
    "pseudo"  =  Pseudos <- Pseudos + 1,
  )
}
Genes <- CDSs + tRNAs + rRNAs + ncRNAs


seq    <- read.fasta(fna)
Name   <- rep(NA, length(seq))
Length <- rep(0, length(seq))
GCcont <- rep(0, length(seq))
GCLen  <- rep(0, length(seq))
for(s in 1:length(seq)){
  Name[s]   <- names(seq[s])
  Length[s] <- length(seq[[s]])
  GCcont[s] <- GC(seq[[s]])
  GCLen[s]  <- GCcont[s] * Length[s]
}

NumCon <- length(seq)
Gname  <- paste(Name, collapse=";")
Gsize  <- sum(Length)
Ggc    <- round((sum(GCLen) / sum(Length) * 100), digits = 2)

write(c("#GenomeSummary:", Tag), file=FileOut,  ncolumns=2, sep="\t", append=T)
write(c("   SeqNames:", Gname ), file=FileOut, ncolumns=2, sep="\t", append=T)
write(c("    NumSeqs:", NumCon ), file=FileOut, ncolumns=2, sep="\t", append=T)
write(c("  TotalSize:", Gsize ), file=FileOut, ncolumns=2, sep="\t", append=T)
write(c("GenomeGC(%):", Ggc ), file=FileOut, ncolumns=2, sep="\t", append=T)
write(c("       Gene:", Genes ), file=FileOut, ncolumns=2, sep="\t", append=T)
write(c("        CDS:", CDSs ), file=FileOut, ncolumns=2, sep="\t", append=T)
write(c("       rRNA:", rRNAs ), file=FileOut, ncolumns=2, sep="\t", append=T)
write(c("       tRNA:", tRNAs ), file=FileOut, ncolumns=2, sep="\t", append=T)
write(c("      ncRNA:", ncRNAs ), file=FileOut, ncolumns=2, sep="\t", append=T)
write(c("     Pseudo:", Pseudos ), file=FileOut, ncolumns=2, sep="\t", append=T)
write(c("    PPgene*:", PPgene), file=FileOut, ncolumns=2, sep="\t", append=T)
write("#*PPgene is included in CDS", file=FileOut, ncolumns=2, sep="\t", append=T)

#END OF FILE
