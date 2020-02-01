.kmer_occurrence_opt_k_m_gen<-function(sequences_to_kmerize,AA_nt,opt_k,opt_m,dictionary_ref){
  #makes a list of counts/frequencies for all kmer length and gap length combinations
  #outputs list of length '#possible combinations'
  #applicable to both nt and Aa counts

  #load relevant dictionary into dictionary variable depending on whether AA or nt
  #skip the first 3 and last 2 amino acids (usually strongly conserved)
  if(AA_nt=="aa"){
    sequences_input<-as.character(sequences_to_kmerize)
    dictionary<-dictionary_ref[[toupper(AA_nt)]]#c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*")
    start_core<-4
    end_core<-2
  } else if(AA_nt=="nt"){
    sequences_input<-as.character(sequences_to_kmerize)
    dictionary<-dictionary_ref[[AA_nt]]#c("A","C","G","T")
    start_core<-10
    end_core<-6
  }

  #set parameters
  m<-opt_m
  k_length<-opt_k
  list_all_dist<-list()
  #for each gap size evaluate occurrences
  for(j in 1:length(m)){
    k_dist<-m[[j]]
    #cat("k: ",k_length,"\t", "m: ", k_dist,"\n")#

    #choose relevant dictionary and kmer pairs (all-pairs_of_combs)
    nt_pair_dict_gen<-dictionary[[as.character(opt_k)]]
    all_pairs_of_combs<-names(nt_pair_dict_gen)
    valid_combo_length<-nchar(all_pairs_of_combs[[1]])

    #count occurrences for each sequence
    found_valid_subsequences<-list()
    for(p in 1:length(sequences_input)){
      #find core substring of current sequence and choose last index that can be evaluated
      curr_seq<-sequences_input[p]
      curr_seq_cut<-base::substr(curr_seq,start_core,nchar(curr_seq)-end_core)
      last_index_in_range <- nchar(curr_seq_cut)-k_dist-k_length #last index is chosen so that all pairs can occur#
      #go through split sequence and find all kmer pairs
      list_found_pairs<- sapply(1:last_index_in_range,function(s) paste(base::substr(curr_seq_cut,s,(s+(k_length-1))),",", base::substr(curr_seq_cut,(s+k_length+k_dist),(s+k_length+k_dist+k_length-1)),sep=""))
      #save the valid combos (of correct length)
      found_valid_subsequences[[p]]<-toupper(list_found_pairs[nchar(list_found_pairs)==valid_combo_length])

    }
    #for current gap size m make table of found sequences
    found_valid_subsequences_table<-table(factor(unlist(found_valid_subsequences),levels=all_pairs_of_combs))

    #construct dataframe containing all patterns, k, m, counts and in-class frequency
    #if count is zero set frequency to zero otherwise evaluate frequency.
    df_occurrence_curr_params<-data.frame(pairs=names(found_valid_subsequences_table),kmer_length=k_length,kmer_dist=k_dist)
    df_occurrence_curr_params$counts<-found_valid_subsequences_table
    if(sum(df_occurrence_curr_params$counts)!=0){
      df_occurrence_curr_params$freq_counts<-df_occurrence_curr_params$counts/sum(df_occurrence_curr_params$counts)
    } else{
      df_occurrence_curr_params$freq_counts<-0
    }
    list_all_dist[[j]]<-df_occurrence_curr_params#

  }
  names(list_all_dist)<-opt_m
  #bind dataframes for all k,m combos and calculate overall frequency (divide by all counts: shorter gaps weighted more than longer ones)
  kmer_k_m<-base::do.call(base::rbind,list_all_dist)
  kmer_k_m$freq_counts_overall<-kmer_k_m$counts/sum(kmer_k_m$counts)

  kmer_k_m
}
