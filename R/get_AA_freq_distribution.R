.get_AA_freq_distribution<-function(repertoire_seqs,length_seqs=0){
  #for lenghts 6:28 returns AA frequency list
  #each list contains a sublist for each length
  #which in turn contains a list AA frequencies per position.

  #initialize vector of amino acids
  AA_list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*")

  #if length_seqs parameter is specifice (!=0) fix chosen length as only unique length
  #else: evaluate lengths 6-28.
  if(length_seqs!=0){
    unique_lengths<-length_seqs
  }else{
    unique_lengths<-c(6:28)
  }

  #for each length evaluate amino acid frequency occurrence across all positions
  list_AA_freqs_per_position<-list()
  for(i in 1:length(unique_lengths)){
    #set curr length and current sequences (seqs of this length)
    curr_length<-unique_lengths[[i]]
    curr_seqs<-as.character(repertoire_seqs[nchar(as.character(repertoire_seqs))==curr_length])

    list_AA_freqs_per_position_temp<-list()
    #check whether there are seqs of this length if yes calc positional frequences
    if(length(curr_seqs)>0){
      for(j in 1:curr_length){
        #split each sequence and unlist amino acids
        pos_ <- sapply(1:length(curr_seqs), function(x) unlist(strsplit(curr_seqs[x],""))[j])
        #evaluate ferquencies per position
        pos_freqs <- table(factor(pos_,levels = AA_list))/length(pos_)
        list_AA_freqs_per_position_temp[[j]] <- pos_freqs
      }
    } else{ #if there are no seqs of this length: set all freqs to 0
      for(j in 1:curr_length){
        pos_ <- rep(0,length(AA_list))
        pos_freqs<-table(factor(pos_,levels=AA_list))
        list_AA_freqs_per_position_temp[[j]] <- pos_freqs
      }
    }
    names(list_AA_freqs_per_position_temp)<-c(1:curr_length)

    list_AA_freqs_per_position[[i]]<-list_AA_freqs_per_position_temp

  }
  names(list_AA_freqs_per_position)<-unique_lengths

  #return. first index: length, second index position.
  return(list_AA_freqs_per_position)
}
