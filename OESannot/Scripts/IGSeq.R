#
#Rで遺伝子間領域の取り出し
 

chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)
if(!require(stringr)) install.packages("stringr")
library(stringr)


args      <- commandArgs(trailingOnly = T)
InFasta   <- args[1]
InGenPos  <- args[2]


fna     <- read.fasta(InFasta)      
genepos <- read.table(InGenPos)          

Contig  <- genepos[,1]
start   <- genepos[,2]
end     <- genepos[,3]

file.create("IGS_Seq.fna", "IGSeq.gff")
file2.out <- "IGS_Seq.fna"

#遺伝子間領域の抜きだし
#intragenic領域長が30nt以上10Knt未満のものを抜きだし
#annotationが配列と逆順の場合も考慮（else if）
#開始・終止コドンを含まないように！

if( start[1] < start[2] ){

  for(i in 1:(nrow(genepos)-1)){
    TContig <- match(Contig[i], names(fna))
    w  <- paste(names(fna[TContig]), end[i], start[i+1],sep = ":")
    if(  ((start[i+1]-end[i]) > 30)&&((start[i+1]-end[i]) < 100000)  ){      
      write.fasta(fna[[TContig]][(end[i]):(start[i+1])], w , file2.out, open = "a")
      igs <- paste(names(fna[TContig]), "IGSeq", "IGS", (end[i]), (start[i+1]), "NA", "+/-", "Draft", sep="\t")
      write(igs, file="IGSeq.gff", ncolumns= 8, append=T)
    }
  }
}else if( start[1] > start[2]){
  for(i in 2:(nrow(genepos)-1)){
    TContig <- match(Contig[i], names(fna))
    w <- paste(names(fna[TContig]), end[i], start[i+1],sep = ":")
    if(  ((start[i-1]-end[i]) > 30)&&((start[i-1]-end[i]) < 100000)  ){
      write.fasta(fna[[TContig]][(start[i+1]):(end[i])], w , file2.out, open = "a")
      igs <- paste(names(fna[TContig]), "IGSeq", "IGS", (start[i+1]), (end[i]), "NA", "+/-", "Draft", sep="\t")
      write(igs, file="IGSeq.gff", ncolumns= 7, append=T)
    }
  }
}


#EOF
