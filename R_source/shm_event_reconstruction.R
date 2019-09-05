#' Decodes immuneSIM repertoire shm_events column.
#' @param shm_event_vec An vector containing strings describing SHM events as output in shm_events column of immuneSIM repertoires.
#' The string contains information on every mutation event in the form:
#'
#' "Position:pre_mutation_nucleotide,post_mutation_nucleotide"
#' combined as: "Mutation1|Mutation2|Mutation3". For example: "171:t,a|186:g,a".
#' @return List of dataframes. Each entry contains location and shm mutation info for a simulated sequence
#' @examples
#' shm_events_example<-c("171:t,a|186:g,a|287:g,a|310:t,c","","294:c,g|316:t,c|330:c,t")
#' shm_list<-shm_event_reconstruction(shm_events_example)

shm_event_reconstruction<-function(shm_event_vec){
  #library(reshape2,quietly=TRUE,warn.conflicts = FALSE)

  list_shm_dfs<-list()
  for(x in 1:length(shm_event_vec)){
    curr_shm_event<-as.character(shm_event_vec[x])
    #check if current sequences had an SHM event
    if(curr_shm_event!=""){
      #decode SHM event string
      curr_df<- reshape2::melt(strsplit(curr_shm_event,"\\|"))
      curr_df$L1<-sapply(1:nrow(curr_df),function(x) strsplit(as.character(curr_df$value),":")[[x]][1])
      curr_df$nt<-sapply(1:nrow(curr_df),function(x) strsplit(as.character(curr_df$value),":")[[x]][2])
      curr_df$nt_pre_shm<-sapply(1:nrow(curr_df),function(x) strsplit(as.character(curr_df$nt),",")[[x]][1])
      curr_df$nt_post_shm<-sapply(1:nrow(curr_df),function(x) strsplit(as.character(curr_df$nt),",")[[x]][2])

      curr_df$value<-NULL
      curr_df$nt<-NULL
      names(curr_df)<-c("location","pre_shm_nt","post_shm_nt")
    }else{
      #if no SHM event occurred return NA,"" for current sequences
      curr_df<-data.frame(location=NA,pre_shm_nt="",post_shm_nt="")
    }
    list_shm_dfs[[x]]<-curr_df
  }

  return(list_shm_dfs)
}
