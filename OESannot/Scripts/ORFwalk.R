##
## ORFwalk.R
##

# Algorithm:
#   1. Check stop codon
#   2. If stop codon(s) was found -> Extend 5' (until stop codon will appear, and find new (most anterior) start codon)
#   3. If no stop codon was found -> Extend 5' (until stop codon will appear)

chooseCRANmirror(graphics=FALSE, ind=49)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)
if(!require(stringr)) install.packages("stringr")
library(stringr)


args <- commandArgs(trailingOnly = T)
gff  <- args[1]
fna  <- args[2]
Ext  <- as.numeric(args[3])
gCode<- args[4]
#Edge <- as.numeric(args[4])


## read infut
GFF  <- as.matrix(read.table(gff, head=F))
Orfs <- subset(GFF, GFF[,3] == "CDS")
Fna  <- read.fasta(fna)



## Procedure
for(i in 1:nrow(Orfs)){
  Seqname  <- names(Fna)
  Contig   <- Orfs[i,1]
  TContig  <- match(Contig, Seqname) # to teat multiple Contigs

  ExFin    <- "F"
  Ext2     <- ((Ext %/% 3)*3)
  Strand   <- Orfs[i,7]
  Frame    <- Orfs[i,8]
  Start    <- as.numeric(Orfs[i,4])
  End      <- as.numeric(Orfs[i,5])
  NoteW    <- paste(Orfs[i,9], "ORFw", sep=";")
  NoteFrag <- paste(Orfs[i,9], "ORFw-Fragmented", sep=";")
  while(ExFin == "F"){  

################################## Forward

    if(Strand == "+" ){                                                                            ## for forward strand
      pOrf    <- Fna[[TContig]][Start:End]
      pCodons <- splitseq(pOrf, frame = 0, word = 3)
	if(gCode == "4" ){
      	  pSCs    <- grep("(taa|tag)",  pCodons) # Stop codons of original ORF: TAA, TAG
        }else{
          pSCs     <- grep("(taa|tag|tga)", pCodons) # Stop codons of extended ORF: TAA, TAG, TGA
        }
      if(length(pSCs) != 0 ){                                                                        # If stop codon(s) was found

        if(Start > Ext2){                                                       

          eOrf    <- Fna[[TContig]][(Start - Ext2):End]                                                    # Extend 5' (until stop codon will appear)
          Codons  <- splitseq(eOrf, frame = 0, word = 3)
		  if(gCode == "4" ){
    	    SCs     <- grep("(taa|tag)",  Codons) # Stop codons of extended ORF: TAA, TAG
		  }else{
	  		SCs     <- grep("(taa|tag|tga)",  Codons) # Stop codons of extended ORF: TAA, TAG, TGA
		  }
          if( length(SCs) > 1 ){
            semiSC  <- SCs[length(SCs) - 1]
            eCodons <- Codons[(semiSC +1):length(Codons)]
            ICs     <- grep("(ttg|ctg|att|atc|ata|atg|gtg)", eCodons) # Start codons
            eIC     <- ICs[1]                                                                        # Find new (most anterior) start codon which makes ORF longest
            eStart  <- ((Start - Ext2) + ((semiSC + eIC - 1) * 3) )   
            if(eStart != Start){
              ExOrf     <- Fna[[TContig]][eStart:End]
              egff      <- paste(Orfs[i,1], Orfs[i,2], Orfs[i,3], eStart, End, Orfs[i,6], Strand,  Frame, NoteW, sep="\t")
              Start     <- eStart
            } else if(eStart == Start){
              egff      <- paste(Orfs[i,1], Orfs[i,2], Orfs[i,3], Start, End, Orfs[i,6], Strand,  Frame, NoteW, sep="\t")
              write(egff, file="ORFw.gff", ncolumns= 9, append=T)
              ExFin     <- "T"
            } 
          } else if( length(SCs) == 1 ){
            egff      <- paste(Orfs[i,1], Orfs[i,2], Orfs[i,3], Start, End, Orfs[i,6], Strand, Frame, NoteW, sep="\t")
            write(egff, file="ORFw.gff", ncolumns= 9, append=T)
            ExFin     <- "T"
          }
        
        } else {
          Status <- "PotentialAnteriorStart"
          egff    <- paste(Orfs[i,1], Orfs[i,2], Orfs[i,3], Start, End, Orfs[i,6], Strand, Frame, "ORFw-PotentialAnteriorStart", sep="\t")
          write(egff, file="ORFw.gff", ncolumns= 9, append=T)
          ExFin     <- "T"
        } 

      }else if(length(pSCs) == 0 ){                                                                # If no stop codon was found
 
        if( (End+Ext2) < length(Fna[[TContig]]) ){

          while(length(pSCs) == 0 ){
            eOrf     <- Fna[[TContig]][Start:(End+Ext2)]                                                   # Extend 3' (until stop codon will appear)
            eCodons  <- splitseq(eOrf, frame = 0, word = 3)
            if(gCode == "4" ){
              eSCs     <- grep("(taa|tag)",  eCodons) # Stop codons of extended ORF: TAA, TAG
            }else{
              eSCs     <- grep("(taa|tag|tga)",  eCodons) # Stop codons of extended ORF: TAA, TAG, TGA
            }
            if(length(eSCs) != 0 ){
              eEnd     <- (Start + (eSCs[1] * 3) -1)
              ExOrf    <- Fna[[TContig]][Start:eEnd]
              pSCs     <- eSCs
              End      <- eEnd            
            } else if(length(eSCs) == 0 ){
              pSCs     <- eSCs
              Ext2     <- (Ext2 + 60)
            }
          }

        } else {
          Status  <- "Fragmented"
          egff    <- paste(Orfs[i,1], Orfs[i,2], Orfs[i,3], Start, End, Orfs[i,6], Strand, Frame, NoteFrag, sep="\t")
          write(egff, file="ORFw.gff", ncolumns= 9, append=T)
          ExFin     <- "T"
        } 
      }


################################## Reverse


    } else if (Strand == "-"){                                                                     ## for reverse strand
      pOrf     <- rev(comp(Fna[[TContig]][Start:End]))
      pCodons  <- splitseq(pOrf, frame = 0, word = 3)
      if(gCode == "4" ){
        pSCs     <- grep("(taa|tag)",  pCodons) # Stop codons of extended ORF: TAA, TAG
      }else{
        pSCs     <- grep("(taa|tag|tga)",  pCodons) # Stop codons of extended ORF: TAA, TAG, TGA
      }
    

      if(length(pSCs) != 0 ){ 

        if( (End + Ext2) < length(Fna[[TContig]]) ){           # If stop codon(s) was found

          eOrf    <- rev(comp(Fna[[TContig]][Start:(End + Ext2)]))     
          Codons  <- splitseq(eOrf, frame = 0, word = 3)       # Extend 5' (until stop codon will appear) 
          if(gCode == "4" ){
            SCs     <- grep("(taa|tag)",  Codons) # Stop codons of extended ORF: TAA, TAG
          }else{
            SCs     <- grep("(taa|tag|tga)",  Codons) # Stop codons of extended ORF: TAA, TAG, TGA
          }
          if( length(SCs) > 1 ){
            semiSC  <- SCs[length(SCs) - 1]
            eCodons <- Codons[(semiSC +1):length(Codons)]
            ICs     <- grep("(ttg|ctg|att|atc|ata|atg|gtg)", eCodons) # Start codons
            eIC     <- ICs[1]                                                                        # Find new (most anterior) start codon which makes ORF longest
            eEnd    <- ((End + Ext2) - ((semiSC + eIC - 1) * 3))
            if(eEnd != End){
              ExOrf     <- rev(comp(Fna[[TContig]][Start:eEnd]))
              egff      <- paste(Orfs[i,1], Orfs[i,2], Orfs[i,3], Start, eEnd, Orfs[i,6], Strand, Frame, NoteW, sep="\t")
              End       <- eEnd
            } else {
              egff    <- paste(Orfs[i,1], Orfs[i,2], Orfs[i,3], Start, End, Orfs[i,6], Strand, Frame, NoteW, sep="\t")
              write(egff, file="ORFw.gff", ncolumns= 9, append=T)
              ExFin     <- "T"
            }
          } else {
            egff    <- paste(Orfs[i,1], Orfs[i,2], Orfs[i,3], Start, End, Orfs[i,6], Strand, Frame, NoteW, sep="\t")
            write(egff, file="ORFw.gff", ncolumns= 9, append=T)
            ExFin     <- "T"
          }

        } else {
          Status <- "PotentialAnteriorStart"
          egff    <- paste(Orfs[i,1], Orfs[i,2], Orfs[i,3], Start, End, Orfs[i,6], Strand, Frame, "ORFw-PotentialAnteriorStart", sep="\t")
          write(egff, file="ORFw.gff", ncolumns= 9, append=T)
          ExFin     <- "T"
        } 

      }else if(length(pSCs) == 0){                                                                   # If no stop codon was found

        if(Start > Ext2 ){                                                                

          while(length(pSCs) == 0 ){
            eOrf     <- rev(comp(Fna[[TContig]][(Start - Ext2):End]))                                      # Extend 3' (until stop codon will appear)
            Codons   <- splitseq(eOrf, frame = 0, word = 3)
            if(gCode == "4" ){
              eSCs     <- grep("(taa|tag)",  Codons) # Stop codons of extended ORF: TAA, TAG
            }else{
              eSCs     <- grep("(taa|tag|tga)",  Codons) # Stop codons of extended ORF: TAA, TAG, TGA
            }
            if(length(eSCs) != 0 ){
              SC1st    <- eSCs[1]
              eCodons  <- Codons[1:SC1st]
              eStart   <- (End - (SC1st * 3) + 1)
              ExOrf    <- rev(comp(Fna[[TContig]][eStart:End]))
              pSCs     <- eSCs
              Start    <- eStart
            } else if(length(eSCs) == 0 ){
              pSCs     <- eSCs
              Ext2     <- (Ext2 + 60)
            } 
          }
        } else {
          Status <- "Fragmented"
          egff    <- paste(Orfs[i,1], Orfs[i,2], Orfs[i,3], Start, End, Orfs[i,6], Strand, Frame, NoteFrag, sep="\t")
          write(egff, file="ORFw.gff", ncolumns= 9, append=T)
          ExFin     <- "T"
        } 

      } 

    }
  }
}

#EOF
