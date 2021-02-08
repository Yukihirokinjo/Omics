#
# GFFmerge.R
#
file.create("tmp_insideORF.gff", "tmp_merged.gff")
preGff <- as.matrix(read.table("tmp_premerge.gff", head=F))

Contigs  <- unique(preGff[,1])
for(c in 1:length(Contigs)){
  gff <- as.matrix(subset(preGff, preGff[,1] == Contigs[c]))

  Fpre <- as.numeric(gff[,4])
  Rpre <- as.numeric(gff[,5])

  Fs     <- sort(Fpre)
  Rs     <- sort(Rpre)
  for(i in 1:nrow(gff)){
    F1     <- as.numeric(gff[i,4])
    F1pos  <- as.numeric(match(F1, Fs))
    R1     <- as.numeric(gff[i,5])
    R1pos  <- as.numeric(match(R1, Rs))
    name   <- as.matrix(gff[i,])

    if((R1pos - F1pos) < 0){
      write(name, file="tmp_insideORF.gff", ncolumns= 9, append=T)
      gff[i,1] <- "NA"
    }
  } # Close iteration

  sortlist  <- order(Fpre, pmax(Fpre, Rpre))
  if(length(sortlist) > 1){ gffs <- gff[sortlist,] }else{ gffs <- t(gff[sortlist,]) }
  gff2  <- as.matrix(subset(gffs, gffs[,1] != "NA"))

  if(nrow(gff2) > 1){
    for(n in 2:nrow(gff2)){
      type1   <- as.character(gff2[n-1,3])
      type2   <- as.character(gff2[n,3])
      Score1  <- strsplit(as.character(gff2[n-1,6]), "\\|")[[1]]
      Score2  <- strsplit(as.character(gff2[n,6]), "\\|")[[1]]
      Note1   <- strsplit(as.character(gff2[n-1,9]), ";")[[1]]
      Note2   <- strsplit(as.character(gff2[n,9]), ";")[[1]]
      Method1 <- as.character(gff2[n-1,2])
      Method2 <- as.character(gff2[n,2])
      strand1 <- as.character(gff2[n-1,7])
      strand2 <- as.character(gff2[n,7])
      F1      <- as.numeric(gff2[n-1,4])
      F2      <- as.numeric(gff2[n,4])
      R1      <- as.numeric(gff2[n-1,5])
      R2      <- as.numeric(gff2[n,5])
      lORF1   <- (R1 - F1)
      lORF2   <- (R2 - F2)

      name1   <- as.matrix(gff2[n-1,])
      name2   <- as.matrix(gff2[n,])
 
      if( strand1 == strand2 ) {
        if( (F1 == F2) || (R1 == R2) ){
          if(lORF1 < lORF2){
            write(name1, file="tmp_insideORF.gff", ncolumns= 9, append=T)
            gff2[n-1,1] <- "NA"
          } else if (lORF1 > lORF2){
            write(name2, file="tmp_insideORF.gff", ncolumns= 9, append=T)
            gff2[n,1]   <- "NA"
          } else if ( (F1 == F2) && (R1 == R2) ){
            if(type1 == "pseudo"){
              write(name2, file="tmp_insideORF.gff", ncolumns= 9, append=T)
              gff2[n,1]   <- "NA"
            } else if(type2 == "pseudo"){
              write(name1, file="tmp_insideORF.gff", ncolumns= 9, append=T)
              gff2[n-1,1] <- "NA"
            } else if(Method1 == "PPgeneDetect" & Method2 != "PPgeneDetect"){
              write(name2, file="tmp_insideORF.gff", ncolumns= 9, append=T)
              gff2[n,1]   <- "NA"
            } else if(Method2 == "PPgeneDetect" & Method1 != "PPgeneDetect"){
              write(name1, file="tmp_insideORF.gff", ncolumns= 9, append=T)
              gff2[n-1,1] <- "NA"
            } else if(Method1 == "IGS_BLASThit" & Method2 != "IGS_BLASThit"){
              write(name1, file="tmp_insideORF.gff", ncolumns= 9, append=T)
              gff2[n-1,1]   <- "NA"
            } else if(Method2 == "IGS_BLASThit" & Method1 != "IGS_BLASThit"){
              write(name2, file="tmp_insideORF.gff", ncolumns= 9, append=T)
              gff2[n,1] <- "NA"
            } else if(length(Score1) > length(Score2)){
              write(name2, file="tmp_insideORF.gff", ncolumns= 9, append=T)
              gff2[n,1]   <- "NA"
            } else if(length(Score1) < length(Score2)){
              write(name1, file="tmp_insideORF.gff", ncolumns= 9, append=T)
              gff2[n-1,1] <- "NA"
            } else if(length(Note1) > length(Note2)){
              write(name2, file="tmp_insideORF.gff", ncolumns= 9, append=T)
              gff2[n,1]   <- "NA"
            } else if(length(Note1) < length(Note2)){
              write(name1, file="tmp_insideORF.gff", ncolumns= 9, append=T)
              gff2[n-1,1] <- "NA"
            }
          }
        }
      }
    } # Close iteration

    gff3<-subset(gff2, gff2[,1]!="NA")
    write.table(gff3, "tmp_merged.gff", quote=F, col.names=F, row.names=F, sep="\t", append=T)
  }else{  # nrow(gff2) == 1
    write.table(gff2, "tmp_merged.gff", quote=F, col.names=F, row.names=F, sep="\t", append=T)

  }
}

# END OF FILE
