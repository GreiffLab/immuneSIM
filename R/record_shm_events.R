#function to identify shm events (output as string that can be easily decoded)
.record_shm_events <- function(seq_a, seq_b){
  #during the shm compare pre and post shm sequences
  if(nchar(seq_a)==nchar(seq_b)){
    #split sequences
    split_seqs <- base::strsplit(c(seq_a, seq_b), split = "")
    #find differences
    compare_seqs <- (split_seqs[[1]] != split_seqs[[2]])
    #save positions where seqs are different
    shm_locations<-which(is.na(compare_seqs) | compare_seqs)
    #if shm occurred save nt pre and post shm and save in df together with position
    if(length(shm_locations)!=0){

      per_loc_seq_a_nt<-split_seqs[[1]][shm_locations]
      per_loc_seq_b_nt<-split_seqs[[2]][shm_locations]

      shm_df<-data.frame(position=shm_locations,nt_pre_shm=per_loc_seq_a_nt,nt_post_shm=per_loc_seq_b_nt)
      #write into encoded string
      shm_df$event<-sapply(1:nrow(shm_df), function(x) paste(shm_df[x,"position"],":",shm_df[x,"nt_pre_shm"],",",shm_df[x,"nt_post_shm"],sep=""))
      shm_record<-paste(shm_df$event,collapse="|")

    }else{
      #if no shm occurred record an empty string
      shm_record <- ""
    }
  }else{
    #if for some reasons two sequences are of unequal length notify user and record SHM as empty string
    base::warning("unequal seq lengths")
    shm_record<-""
  }

  return(shm_record)
}
