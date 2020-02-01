#' Simulates an immune repertoire based on user-defined parameters
#'
#' @param number_of_seqs Integer defining the number of sequences that should be simulated
#' @param vdj_list List containing germline genes and their frequencies
#' @param species String defining species for which repertoire should be simulated ("mm": mouse, "hs": human. Default: "mm").
#' @param receptor String defining receptor type ("ig" or "tr". Default: "ig")
#' @param chain String defining chain (for ig: "h","k","l", for tr: "b" or "a". Default: "h")
#' @param insertions_and_deletion_lengths Data.frame containing np1, np2 sequences as well as deletion lengths.
#' (Pooled from murine repertoire data, Greiff,2017)
#' Note: This is a subset of 500000 observations of the dataframe used in the paper. The full dataframe which can
#' be introduced here can be found on: (Git-Link)
#' @param user_defined_alpha Numeric. Scaling parameter used for the simulation of powerlaw distribution
#' (recommended range 2-5. Default: 2, https://en.wikipedia.org/wiki/Power_law)
#' @param name_repertoire String defining chosen repertoire name recorded in the name_repertoire column of the output for identification.
#' @param length_distribution_rand Vector containing lengths of immune receptor sequences based on immune repertoire data (Greiff, 2017).
#' @param random Boolean. If TRUE repertoire will consist of fully random sequences, independent of germline genes.
#' @param shm.mode String defining mode of somatic hypermutation simulation based on AbSim
#' (options: 'none', 'data','poisson', 'naive', 'motif', 'wrc'. Default: 'none'). See AbSim documentation.
#' @param shm.prob Numeric defining probability of a SHM (somatic hypermutation) occurring at each position.
#' @param vdj_noise Numeric between 0,1, setting noise level to be introduced in provided V,D,J germline frequencies. 0 denotes no noise. (Default: 0)
#' @param vdj_dropout Named vector containing entries V,D,J setting the number of germline genes to be dropped out. (Default: c("V"=0,"D"=0,"J"=0))
#' @param ins_del_dropout String determining whether insertions and deletions should occur.
#' Options: "", "no_insertions", "no_insertions_n1", "no_insertions_n2", "no_deletions_v", "no_deletions_d_5",
#' "no_deletions_d_3", "no_deletions_j", "no_deletions_vd", "no_deletions". Default: "")
#' @param equal_cc Boolean that if set TRUE will override user_defined_alpha and generate a clone count distribution that is equal for all sequences.
#' Default: FALSE.
#' @param freq_update_time Numeric determining whether simulated VDJ frequencies agree with input after set amount of sequences to correct for VDJ bias.
#' Default: Update after 50 percent of sequences.
#' @param max_cdr3_length Numeric defining maximal length of cdr3. (Default: 100)
#' @param min_cdr3_length Numeric defining minimal length of cdr3. (Default: 6)
#' @param verbose Boolean toggling printing of progress on and off (Default: FALSE)
#' @param airr_compliant Boolean determining whether output repertoire should be named in an AIRR compliant manner
#' (Default: TRUE). (http://docs.airr-community.org/en/latest/)
#' @return An annotated AIRR-compliant immuneSIM repertoire. (http://docs.airr-community.org/en/latest/)
#' @examples
#' sim_rep <- immuneSIM(number_of_seqs = 10, vdj_list = list_germline_genes_allele_01,
#' species = "mm", receptor = "ig", chain = "h",
#' insertions_and_deletion_lengths = insertions_and_deletion_lengths_df,
#' user_defined_alpha = 2,name_repertoire = "mm_igh_sim",
#' shm.mode = "data",shm.prob=15/350,vdj_noise = 0, vdj_dropout = c(V=0,D=0,J=0),
#' ins_del_dropout = "",min_cdr3_length = 6)


immuneSIM<-function(number_of_seqs=1000,vdj_list=list_germline_genes_allele_01,species="mm",receptor="ig",chain='h',insertions_and_deletion_lengths=insertions_and_deletion_lengths_df,user_defined_alpha=2,name_repertoire="sim_rep",length_distribution_rand=length_dist_simulation,random=FALSE,shm.mode="none",shm.prob=15/350,vdj_noise=0,vdj_dropout=c("V"=0,"D"=0,"J"=0),ins_del_dropout=c(""),equal_cc=FALSE,freq_update_time=round(0.5*number_of_seqs),max_cdr3_length=100,min_cdr3_length=6, verbose=TRUE, airr_compliant=TRUE){
  #vdj_noise determines how much the frequencies get altered: noise_level btw 0,1
  #noise level determines sd of randomly distributed noise term. 1 super noisey. recommended:0.2

  if(random==FALSE){
    #find relevant list of germline genes.
    vdj_list<-vdj_list[[species]][[receptor]][[chain]]

    #introduce noise into VDJ if vdj_noise==TRUE
    if(vdj_noise!=0){
      vdj_list<-.noiseify_VDJ_freqs(vdj_list=vdj_list,vdj_noise_level=vdj_noise)
    }

    #dropout VDJ genes if chosen by user
    vdj_list <- .dropout_VDJ(vdj_list,dropout_vdj=vdj_dropout)

    #modify insertion deletion pool (option to exclude insertions or deletions)
    if(ins_del_dropout[1]!=""){
      insertions_and_deletion_lengths<-.new_insertion_deletion_df(new_insertion_deletion=insertions_and_deletion_lengths,options=ins_del_dropout) #new_insertion_deletion_df(insertions_and_deletion_lengths=insertions_and_deletion_lengths,ins_del_dropout=ins_del_dropout)
    }

    #simulating repertoire
    sim_repertoire<-.immuneSim_standard(number_of_seqs=number_of_seqs,
                                       vdj_list=vdj_list,
                                       species=species,
                                       receptor=receptor,
                                       chain=chain,
                                       insertions_and_deletion_lengths=insertions_and_deletion_lengths,
                                       user_defined_alpha=user_defined_alpha,
                                       name_repertoire=name_repertoire,
                                       freq_update_time=freq_update_time,
                                       max_cdr3_length=max_cdr3_length,
                                       min_cdr3_length=min_cdr3_length,
                                       shm.mode=shm.mode,
                                       shm.prob=shm.prob,
                                       verbose=verbose)

  }else{
    #simulate fully random repertoire
    sim_repertoire<-.immuneSim_random(nb_of_seqs=number_of_seqs,length_distribution_rand=length_dist_simulation,name_repertoire=name_repertoire)
  }

  #overlay clone count
  sim_repertoire<-.add_clone_count(simulated_repertoire=sim_repertoire,user_defined_alpha=user_defined_alpha,equal_cc=equal_cc)

  #if airr_compliant naming is chosen rename columns
  if(airr_compliant==TRUE){
    names(sim_repertoire) <- c("sequence", "sequence_aa", "junction", "junction_aa", "v_call", "d_call",
                               "j_call", "np1", "np2", "del_v", "del_d_5", "del_d_3", "del_j",
                               "v_sequence_alignment", "d_sequence_alignment", "j_sequence_alignment", "freqs","counts","shm_events", "name_repertoire")
  }

  return(sim_repertoire)
}

