# Generates a repertoire containing fully random sequences.
#
# @param nb_of_seqs Defines the number of sequences that should be simulated
# @param length_distribution_rand Desired length distribution for generated repertoire
# @param user_defined_alpha Parameter that for the simulation of powerlaw distribution (recommended range 2-5)
# @param name_repertoire Chosen repertoire name recorded in the name_repertoire column of the output for identification
# @return Random immuneSIM repertoire
# @examples
# immuneSim_random(10, length_dist_simulation, 2, "random_repertoire")


#load("length_dist_simulation") #for lengths of fully immuneSIM_random: length distribution from a 10000seq simualted repertoire

.immuneSim_random<-function(nb_of_seqs,length_distribution_rand=length_dist_simulation,user_defined_alpha=2,name_repertoire){
  #set number of seqs and length distribution that simulation should adhere to.
  #creates a vector of length nb_of seqs with lenghts (integers).
  number_of_simulated_seqs<-nb_of_seqs
  length_distribution_sampled<-sample(length_distribution_rand,number_of_simulated_seqs,replace=TRUE)
  #for each entry in the length_distribution_sampledcreate a random nucleotide and a random amino acid sequence
  sim_full_random_nt<-list()
  sim_full_random_aa<-list()
  for(i in 1:length(length_distribution_sampled)){
    sim_full_random_nt[[i]]<-paste(sample(c("a","c","t","g"),length_distribution_sampled[i],replace=TRUE),collapse="")
    sim_full_random_aa[[i]]<-paste(sample(c('A','V','I','L','M','F','Y','W','C','G','P','S','T','N','Q','R','H','K','D','E'),length_distribution_sampled[i]/3,replace=TRUE),collapse="")
  }
  sim_full_random_nt<-as.character(unlist(sim_full_random_nt))
  sim_full_random_aa<-as.character(unlist(sim_full_random_aa))

  #build data.frame that contains random sequence data and vdj_calls 'random' and name. set the rest to NA
  sim_full_random<-data.frame(sequence=sim_full_random_nt,
                              sequence_aa=sim_full_random_aa,
                              junction="",
                              junction_aa="",
                              v_call="IGHVrandom",
                              d_call="IGHDrandom",
                              j_call="IGHJrandom",
                              np1=NA,
                              np2=NA,
                              del_v=NA,
                              del_d_5=NA,
                              del_d_3=NA,
                              del_j=NA,
                              v_sequence_alignment =NA,
                              d_sequence_alignment=NA,
                              j_sequence_alignment=NA,
                              freqs=NA,
                              counts=NA,
                              shm_events=NA,
                              name_repertoire=name_repertoire)	#set sim to 0 for random, 1 for non random


  return(sim_full_random)

}
