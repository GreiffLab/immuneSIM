#' Deletes top hub sequences from repertoire, changing the network architecture.
#'
#' @param repertoire An annotated AIRR compliant repertoire.
#'
#' (http://docs.airr-community.org/en/latest/)
#' @param top_x Determines what percentage of hub sequences get excluded
#'
#' (Default: 0.005, i.e. Top 0.5 percent)
#' @param report The user can choose to output a report csv file containing the excluded
#' sequences. (Default: FALSE)
#' @param output_dir If user specifies and output directory a csv file containing
#' the excluded sequences is saved at that path, otherwise it will be saved in tempdir().
#' @param verbose Determines whether messages on plot locations are output to user. (Default: TRUE)
#' @return Repertoire reduced by hub sequence (new network architecture)
#' @examples
#' repertoire <- list_example_repertoires[["example_repertoire_A"]]
#' rep_excluded_hubs <- hub_seqs_exclusion(repertoire, top_x = 0.005, output_dir = "")


hub_seqs_exclusion<-function(repertoire, top_x = 0.005, report = FALSE, output_dir = "",verbose = TRUE){

  #set output folder
  if(output_dir == ""){
    wd <- tempdir()
  }else{
    wd <- output_dir
  }

  ###construct similarity network and mutate hub sequences.
  #prep sequences for graph / get unique sequences
  CDR3s<-unique(as.character(repertoire$junction_aa))

  #set similarity threshold here levenshtein distance 1
  similarity<-1

  #calculate levenshtein distance matrix and name rows/columns
  ld_matrix <- stringdist::stringdistmatrix(CDR3s,CDR3s,method="lv")
  rownames(ld_matrix)<-CDR3s
  colnames(ld_matrix)<-CDR3s

  #turn ld_matrix into adjacency matrix
  adj_matrix<-ld_matrix
  adj_matrix[adj_matrix!=similarity]<-0
  adj_matrix[adj_matrix==similarity]<-1

  #construct graph
  CDR3s_graph<-igraph::graph_from_adjacency_matrix(adj_matrix,mode=c('undirected'))

  #evaluate hub scores and find seqs to be excluded
  hub_scores_cdr3s<-igraph::hub_score(CDR3s_graph, weights=NA)$vector
  seqs_to_be_excluded<-names(utils::tail(sort(hub_scores_cdr3s),round(length(hub_scores_cdr3s)*top_x)))

  #find sequences to keep
  kept_CDR3s<-CDR3s[!(CDR3s %in% seqs_to_be_excluded)]

  #subset repertoire to only include kept sequences.
  repertoire_hubs_excluded<-repertoire[repertoire$junction_aa %in% kept_CDR3s,]

  #if report is chosen extract excluded sequences and write as .csv file.
  if(report == TRUE){
    report_excluded<-repertoire[repertoire$junction_aa %in% seqs_to_be_excluded,]
    write.csv(report_excluded,file=file.path(wd, paste(repertoire$name_repertoire[1],'_hub_seqs_excluded.csv',sep="")))

    if(verbose==TRUE){
      cat("csv file containing excluded sequences has been saved to: ",wd,"\n")
    }
  }

  return(repertoire_hubs_excluded)

}

