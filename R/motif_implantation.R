#' Implant random or predefined motifs into CDR3
#'
#' @param sim_repertoire An annotated AIRR compliant immuneSIM repertoire.
#' @param motif Either a list that contains number, length and frequencies of motifs or dataframe that contains predefined motifs and their frequencies
#' @param fixed_pos defines position at which motif is to be introduced. if 0 motif will be introduced at random position
#' @return Repertoire with modified sequences containing implanted motifs in CDR3.
#' @examples
#' sim_repertoire <- list_example_repertoires[["example_repertoire_A"]]
#' sim_rep_motifs <- motif_implantation(sim_repertoire,list("n"=2,"k"=3,"freq"=c(0.1,0.1)),0)

motif_implantation <- function(sim_repertoire,motif,fixed_pos=0){
  #library(plyr,quietly=TRUE,warn.conflicts = FALSE)
  #library(stringr,quietly=TRUE,warn.conflicts = FALSE)

  #reformat every sequence to character to preempt errors
  sim_repertoire$junction_aa <- as.character(sim_repertoire$junction_aa)
  sim_repertoire$junction <- as.character(sim_repertoire$junction)
  sim_repertoire$sequence<- as.character(sim_repertoire$sequence)
  sim_repertoire$sequence_aa<- as.character(sim_repertoire$sequence_aa)
  sim_repertoire$np1 <- as.character(sim_repertoire$np1)
  sim_repertoire$np2 <- as.character(sim_repertoire$np2)
  sim_repertoire$v_sequence_alignment  <- as.character(sim_repertoire$v_sequence_alignment)
  sim_repertoire$d_sequence_alignment  <- as.character(sim_repertoire$d_sequence_alignment)
  sim_repertoire$j_sequence_alignment  <- as.character(sim_repertoire$j_sequence_alignment)

  #check if motif is given as a list or dataframe
  #if list: generate motifs randomly
  #if dataframe: take in user defined motifs
  if(class(motif)=="list"){
    #read in variables
    length_of_motif <- motif[['k']]
    number_of_motifs <- motif[['n']]
    motif_freq <- motif[['freq']]

    #generate motif to be injected
    list_motif <- list()
    list_motif[["aa"]] <- list()
    list_motif[["nt"]] <- list()
    for(i in 1:number_of_motifs){
      list_motif[["aa"]][[i]] <- paste(gen_code[sample(1:nrow(gen_code),length_of_motif),"aa"],collapse="")
      list_motif[["nt"]][[i]] <- paste(gen_code[sample(1:nrow(gen_code),length_of_motif),"codon"],collapse="")
    }

    aa_motif<-sapply(1:number_of_motifs,function(z) paste(unlist(list_motif[["aa"]][[z]]),collapse=""))
    nt_motif<-tolower(sapply(1:number_of_motifs,function(z) paste(unlist(list_motif[["nt"]][[z]]),collapse="")))

    #make dataframe with motif information
    motif_df<-data.frame(aa=aa_motif,nt=nt_motif,freq=motif_freq*1/(sum(motif_freq)))

  }else if(class(motif)=="data.frame"){
    #load dataframe with motif information
    motif_df<-motif

    #get numbers relevant
    length_of_motif <- nchar(as.character(motif_df$aa[[1]]))
    number_of_motifs <- length(unique(as.character(motif_df$aa)))
    motif_freq <- motif_df$freq

    #scale freqs up to 1
    motif_df$freq<-motif_df$freq*1/(sum(motif_df$freq))
  }

  #calculate number of sequences that will receive motif injection
  nb_seqs_with_motif<-round(nrow(sim_repertoire)*sum(motif_freq))
  #for each target choose a motif to inject
  motif_per_target<-sample(motif_df$aa,prob=motif_df$freq,nb_seqs_with_motif,replace=TRUE)

  #sample rows that will receive injection and make 'injection guide' df that contains all information
  sample_inj_target <- sample(1:nrow(sim_repertoire),nb_seqs_with_motif)
  sample_inj_target_df<-data.frame(target_seq=sample_inj_target,aa=motif_per_target)
  motif_inj_target_df<-plyr::join(sample_inj_target_df,motif_df)

  #initialize new columns that will receive info on injection if used
  sim_repertoire$motif_aa<-""
  sim_repertoire$motif_nt<-""
  sim_repertoire$motif_pos<-0

  #for each target find an injection position
  for(i in 1:nrow(motif_inj_target_df)){
    #find target and starting point
    target<-sim_repertoire[motif_inj_target_df$target_seq[[i]],]

    #find motifs for this target
    motif_curr_aa<-as.character(motif_inj_target_df$aa[[i]])
    motif_curr_nt<-as.character(motif_inj_target_df$nt[[i]])

    #find injection position (check if fixed_position given by user)
    if(fixed_pos==0){
      starting_point<-sample(1:(nchar(as.character(target$junction_aa))-2-length_of_motif),1)
    }else{
      if(fixed_pos %in% c(1:(nchar(as.character(target$junction_aa))-2-length_of_motif))){
        starting_point <- fixed_pos
      }else{
        starting_point<-sample(1:(nchar(as.character(target$junction_aa))-2-length_of_motif),1)
      }
    }
    starting_point_nt<-(starting_point*3)-2

    #save initial seq for targeting full seq
    initial_sequence<-as.character(target$junction_aa)
    initial_sequence_nt<-as.character(target$junction)
    initial_full_v_aa<-stringr::str_split(as.character(target$sequence_aa), initial_sequence)[[1]][1]
    initial_full_j_aa<-stringr::str_split(as.character(target$sequence_aa), initial_sequence)[[1]][2]
    initial_full_v_nt<-stringr::str_split(as.character(target$sequence), initial_sequence_nt)[[1]][1]
    initial_full_j_nt<-stringr::str_split(as.character(target$sequence), initial_sequence_nt)[[1]][2]

    #inject aa motif
    substr(target$junction_aa,starting_point,starting_point+length_of_motif-1) <- motif_curr_aa
    #inject nt motif
    substr(target$junction,starting_point_nt,starting_point_nt+(3*length_of_motif)-3) <- motif_curr_nt

    #change full length_sequences
    target$sequence_aa<-paste(initial_full_v_aa,target$junction_aa,initial_full_j_aa,sep="")
    target$sequence<-paste(initial_full_v_nt,target$junction,initial_full_j_nt,sep="")

    #change jct subsequences v,d,j,n1,n2
    target$v_sequence_alignment<-substr(target$junction,1,nchar(target$v_sequence_alignment))
    n1_d_n2_j<-substr(target$junction,nchar(target$v_sequence_alignment)+1,nchar(target$junction))
    target$np1<-substr(n1_d_n2_j,1,nchar(target$np1))
    target$d_sequence_alignment<-substr(n1_d_n2_j,nchar(target$np1)+1,nchar(target$np1)+nchar(target$d_sequence_alignment))
    target$np2<-substr(n1_d_n2_j,nchar(target$np1)+nchar(target$d_sequence_alignment)+1,nchar(target$np1)+nchar(target$d_sequence_alignment)+nchar(target$np2))

    target$j_sequence_alignment<-substr(n1_d_n2_j,nchar(target$np1)+nchar(target$d_sequence_alignment)+nchar(target$np2)+1,nchar(n1_d_n2_j))

    #save motif injection metainformation (what was injected and where)
    target$motif_aa<-motif_curr_aa
    target$motif_nt<-motif_curr_nt
    target$motif_pos<-starting_point

    #reintroduce motif back into dataframe
    sim_repertoire[motif_inj_target_df$target_seq[[i]],]<-target
  }

  return(sim_repertoire)

}
