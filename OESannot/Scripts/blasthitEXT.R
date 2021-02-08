#
# blasthitEXT.R [output_prefix] [Sequence_Type]
#

chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)


args    <- commandArgs(trailingOnly = T)
Name    <- as.character(args[1])
Type    <- as.character(args[2])
gCode	<- as.numeric(args[3])


hit <- read.table("tmp_BLASThit.txt", head=F)
TGB <- read.table("tmp_TGBhit.txt", head=F)
fna <- read.fasta("tmp_genome.fna")

OPname    <-  paste(Name, "BLASThit", sep="_")
FnameGFF  <- "tmp_BLASThit.gff"
FnameFNA  <- "tmp_BLASThit.fna"
file.create(FnameGFF, FnameFNA)

if(Type == "igs"){

  for( i in 1:nrow(hit)){
    seqtag  <- strsplit(as.character(hit[i,1]), ":")[[1]] 
    Contig  <- seqtag[1]
    Seqname <- names(fna)
    TContig <- match(Contig, Seqname)

    # to get hit reagion (blast-hit position to original genomic position )
    iStart  <- as.numeric(seqtag[2])
    iEnd    <- as.numeric(seqtag[3])
    hStart  <- as.numeric(hit[i,2])
    hEnd    <- as.numeric(hit[i,3])

    if(hStart < hEnd ){       
      Start  <- (iStart + hStart -1 )
      End    <- (iStart + hEnd -1 )  # not include (true) stop cpdon
      IGSname <- paste(names(fna[TContig]), OPname, "CDS", Start, End, "NA", "+", "0", "candidate", sep="\t")
      w   <- paste(Name, "_", i, ":", Start, "..", End, sep="")
      tORF<- translate(fna[[TContig]][Start:End], frame=0, numcode=gCode)
    } else if (hStart > hEnd){
      Start  <- (iStart + hEnd +1)   # not include (true) stop cpdon
      End    <- (iStart + hStart -1)
      IGSname <- paste(names(fna[TContig]), OPname, "CDS", Start, End, "NA", "-", "0", "candidate", sep="\t")
      w   <- paste(Name, "_", i, ":", "c(", Start, "..", End, ")", sep="")
      tORF<- translate(rev(comp(fna[[TContig]][Start:End])), frame=0, numcode=gCode)
    }

    stopC <- grep("\\*", tORF)

    if(length(stopC) == 0){
      write.fasta(tORF, w , file.out=FnameFNA, open = "a")
      write(IGSname, file=FnameGFF, ncolumns= 9, append=T)
    }
    
  }

}else if (Type == "prot"){

  prot <- read.fasta("tmp_proteome.faa")

  for( i in 1:nrow(hit)){
    prottag <- as.character(hit[i,1])
    Contig  <- strsplit(as.character(prottag), ":")[[1]][1] 
    Seqname <- names(fna)
    TContig <- match(Contig, Seqname)

    pStart   <- as.numeric(TGB[TGB$V1==prottag, 2])
    pEnd     <- as.numeric(TGB[TGB$V1==prottag, 3])
    # to avoid multihit in genome back blast (tblastn) *NOTE in tblastn, num_alignments 1 dose not work!
    Start    <- pStart[1]
    End      <- pEnd[1]
  
    ## prepare GFF file
    # strand information is determined by the order of Start End position in blast output (i.e. "+": Start End, "-": End Start)
    if(Start < End ){
      ProName <- paste(names(fna[TContig]), OPname, "CDS", Start, End, "NA", "+", "0", "candidate", sep="\t")
      w   <- paste(prottag, " ", Start, "..", End, sep="")
    } else if (Start > End){
      ProName <- paste(names(fna[TContig]), OPname, "CDS", End, Start, "NA", "-", "0", "candidate", sep="\t")
      w   <- paste(prottag, " ", "c(", End, "..", Start, ")", sep="")
    }
    
    write(ProName, file=FnameGFF, ncolumns= 9, append=T)

  }

}

