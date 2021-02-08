#
# ProdiLCORFremove.R
#


chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)


LCname <- as.matrix(read.table("tmp_prodiLowConf_Names.txt",head=F))
seq    <- read.fasta("tmp_prodi.faa")

for(i in 1:length(seq)){
  if( !is.na(match(names(seq[i]),LCname)) ){
    write.fasta(seq[i],names(seq[i]), file.out="prodi_LowConf.faa", open="a")
  }
} 

