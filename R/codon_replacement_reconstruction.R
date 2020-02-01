#' Decodes immuneSIM repertoire codon replacements events.
#' @param codon_replacement_vec An vector containing strings describing codon replacement events as generated
#' by codon_replacement() function. The string contains information on every replacement event in the form:
#'
#' "initial_codon:replacement_codon:number_of_occurrences"
#'
#' which is combined into: "Replacement1|Replacement2|Replacement3".
#'
#' (For example: "tac,tat:3|agc,agt:1|gtg,gtt:0".)
#' @return List of dataframes. Each entry contains replacement info including count of occurrences for each simulated sequence.
#' @examples
#' codon_replacement_example <- c("tat,tac:3|agt,agc:3|gtt,gtg:0", "tat,tac:1|agt,agc:1|gtt,gtg:1")
#' codon_replacement_list <- codon_replacement_reconstruction(codon_replacement_example)


codon_replacement_reconstruction<-function(codon_replacement_vec){
  #library(reshape2,quietly=TRUE,warn.conflicts = FALSE)

  list_cr_dfs<-list()
  for(x in 1:length(codon_replacement_vec)){
    curr_cr_event<-as.character(codon_replacement_vec[x])
    #check if current sequences had an SHM event
    if(curr_cr_event!=""){
      #decode SHM event string
      curr_df<- reshape2::melt(strsplit(curr_cr_event,"\\|"))
      curr_df$nt<-sapply(1:nrow(curr_df),function(x) strsplit(as.character(curr_df$value),":")[[x]][1])
      curr_df$L1<-sapply(1:nrow(curr_df),function(x) strsplit(as.character(curr_df$value),":")[[x]][2])
      curr_df$initial<-sapply(1:nrow(curr_df),function(x) strsplit(as.character(curr_df$nt),",")[[x]][1])
      curr_df$replacement<-sapply(1:nrow(curr_df),function(x) strsplit(as.character(curr_df$nt),",")[[x]][2])
      #get rid of duplicate columns and rename
      curr_df$value<-NULL
      curr_df$nt<-NULL
      names(curr_df)<-c("occurrences","intial","replacement")
      curr_df$sequence_id<-x
    }else{
      #if no SHM event occurred return NA,"" for current sequences
      curr_df<-data.frame(occurrences=NA,intial=NA,replacement=NA,sequence_id=x)
    }
    list_cr_dfs[[x]]<-curr_df
  }

  return(list_cr_dfs)
}
