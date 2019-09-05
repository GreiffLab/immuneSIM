#calculate vdj occurence function
.calc_vdj_occ<-function(repertoire_in,name){
  #make table of v,d,j_calls and calculate frequencies
  occurrence_v <- as.data.frame(table(repertoire_in$v_call))
  occurrence_v$Freq <- occurrence_v$Freq/sum(occurrence_v$Freq)
  occurrence_v$gene <- "V"

  occurrence_d <- as.data.frame(table(repertoire_in$d_call))
  occurrence_d$Freq <- occurrence_d$Freq/sum(occurrence_d$Freq)
  occurrence_d$gene <- "D"

  occurrence_j <- as.data.frame(table(repertoire_in$j_call))
  occurrence_j$Freq <- occurrence_j$Freq/sum(occurrence_j$Freq)
  occurrence_j$gene <- "J"

  vdj_occurrence <- rbind(occurrence_v,occurrence_d,occurrence_j)
  vdj_occurrence$repertoire <- name

  vdj_occurrence_df <- vdj_occurrence
  vdj_occurrence_df$gene <- factor(vdj_occurrence_df$gene,levels=c("V","D","J"))

  return(vdj_occurrence_df)
}
