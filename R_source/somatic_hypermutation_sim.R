# Somatically hypermutates full VDJ sequences (modified from AbSim package)
#
# @param vdj_seq An annotated AIRR compliant immuneSIM repertoire.
# @param mut_param Def
# @param v_seq Def
# @param d_seq Def
# @param j_seq Def
# @param SHM.nuc.prob Def

# @return full VDJ sequence that is somatically hypermutated.
# @examples
# motif_implantation(sim_repertoire,list("n"=2,"k"=3,"freq"=c(0.1,0.1)),0)
#

.somatic_hypermutation_sim<-function(vdj_seq, mut_param="data",v_seq,d_seq,j_seq, SHM.nuc.prob=15/350){

  #FROM ABSIM PACKAGE DOC: SHM.method The mode of SHM speciation events. Options are either: "poisson","data","motif","wrc",
  #and "all". Specifying either "poisson" or "naive" will result in mutations that can
  #occur anywhere in the heavy chain region, with each nucleotide having an equal
  #probability for a mutation event. Specifying "data" focuses mutation events during
  #SHM in the CDR regions (based on IMGT), and there will be an increased
  #probability for transitions (and decreased probability for transversions). Specifying
  #"motif" will cause neighbor dependent mutations based on a mutational
  #matrix from high throughput sequencing data sets (Yaari et al., Frontiers in Immunology,
  #2013). "wrc" allows for only the WRC mutational hotspots to be
  #included (where W equals A or T and R equals A or G). Specifying "all" will
  #use all four types of mutations during SHM branching events, where the weights
  #for each can be specified in the "SHM.nuc.prob" parameter.

  #mut param just cdr3: data,
  #mutparam full seq: poisson, naive, motif, wrc
  #require(AbSim)

  base_line_mutations <- 0
  if(mut_param=="naive" || mut_param=="all" ||mut_param=="poisson"){
    holding_mut <- sample(x=c(0,1), nchar(vdj_seq), replace=TRUE,
                          c(SHM.nuc.prob[1], 1-SHM.nuc.prob[1]))
    for (i in 1:nchar(vdj_seq)){
      if(holding_mut[i]==0){
        holding_char <- substr(vdj_seq, i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(rep(1/3,3)))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(rep(1/3,3)))
        #substr(vdj_seq, i, i) <- sample(x=c("A", "T", "G","C"), 1, replace=TRUE, c(rep(.25,4)))
        base_line_mutations <- base_line_mutations + 1
      }
    }
  }
  if(mut_param=="data" || mut_param=="all"){
    index_CDR_start <- nchar(v_seq)-15
    index_CDR_stop <- nchar(vdj_seq)
    CDR_length <- index_CDR_stop - index_CDR_start
    if(mut_param=="data") CDR_prob <- SHM.nuc.prob
    else if(mut_param=="all") CDR_prob <- SHM.nuc.prob[2]
    no_CDR_prob <- 1-CDR_prob
    #CDR_mut <- sample(x=c(0,1), nchar(CDR_length), replace=TRUE, c(no_CDR_prob,CDR_prob))
    for(i in 81:114){
      CDR_mut <- sample(x=c(0,1), 1, replace=TRUE, c(no_CDR_prob,CDR_prob))
      if(CDR_mut==1){
        holding_char <- substr(vdj_seq,i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(15,70,15))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(15,15,70))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(70,15,15))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(15,15,75))
        #substr(vdj_seq,i,i) <- sample(x=c("A","T","G","C"),1,replace=TRUE, c(rep(.25,4)))
      }
    }
    for(i in 168:195){
      CDR_mut <- sample(x=c(0,1), 1, replace=TRUE, c(no_CDR_prob,CDR_prob))
      if(CDR_mut==1){
        holding_char <- substr(vdj_seq,i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(15,70,15))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(15,15,70))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(70,15,15))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(15,15,75))
        #substr(vdj_seq,i,i) <- sample(x=c("A","T","G","C"),1,replace=TRUE, c(rep(.25,4)))
      }
    }
    for(i in index_CDR_start:index_CDR_stop){
      CDR_mut <- sample(x=c(0,1), 1, replace=TRUE, c(no_CDR_prob,CDR_prob))
      if(CDR_mut==1){
        holding_char <- substr(vdj_seq,i,i)
        if(holding_char=="a") substr(vdj_seq, i, i) <- sample(x=c("t", "g","c"), 1, replace=TRUE, c(15,70,15))
        else if(holding_char=="t") substr(vdj_seq, i, i) <- sample(x=c("a", "g","c"), 1, replace=TRUE, c(15,15,70))
        else if(holding_char=="g") substr(vdj_seq, i, i) <- sample(x=c("a", "t","c"), 1, replace=TRUE, c(70,15,15))
        else if(holding_char=="c") substr(vdj_seq, i, i) <- sample(x=c("a", "g","t"), 1, replace=TRUE, c(15,15,75))
        #substr(vdj_seq,i,i) <- sample(x=c("A","T","G","C"),1,replace=TRUE, c(rep(.25,4)))
      }
    }
  }
  if(mut_param == "motif" || mut_param == "all"){
    random_spot <- hotspot_df[sample(nrow(hotspot_df)),]
    current_hot_spots <- 0
    if(mut_param=="motif") hot_spot_limit <- SHM.nuc.prob
    else if(mut_param=="all") hot_spot_limit <- SHM.nuc.prob[3]
    #hot_spot_limit <- SHM.nuc.prob
    for(i in 1:nrow(random_spot)){
      if(grepl(pattern = random_spot$pattern[i],x = vdj_seq)){
        current_hot_spots <- current_hot_spots+1
      }
      if(current_hot_spots>hot_spot_limit) break
      vdj_seq <- base::sub(random_spot$pattern[i],replacement = paste(substring(random_spot$pattern[i],first=1,last=2),sample(x = c("A","C","G","T"),1,replace=TRUE,prob = c(random_spot$toA[i], random_spot$toC[i], random_spot$toG[i], random_spot$toT[i])),substring(random_spot$pattern[i],first=4,last=5),sep=""),x = vdj_seq)

    }
  }
  if(mut_param == "wrc" || mut_param == "all"){
    random_spot <- one_spot_df[sample(nrow(one_spot_df)),]
    current_hot_spots <- 0
    if(mut_param=="wrc") hot_spot_limit <- SHM.nuc.prob
    else if(mut_param=="all") hot_spot_limit <- SHM.nuc.prob[4]
    #hot_spot_limit <- SHM.nuc.prob
    for(i in 1:nrow(random_spot)){
      if(grepl(pattern = random_spot$pattern[i],x = vdj_seq)){
        current_hot_spots <- current_hot_spots+1
      }
      if(current_hot_spots>hot_spot_limit) break
      vdj_seq <- base::sub(random_spot$pattern[i],replacement = paste(substring(random_spot$pattern[i],first=1,last=2),sample(x = c("A","C","G","T"),1,replace=TRUE,prob = c(random_spot$toA[i], random_spot$toC[i], random_spot$toG[i], random_spot$toT[i])),substring(random_spot$pattern[i],first=4,last=5),sep=""),x = vdj_seq)
    }
  }
  return(vdj_seq)

}
