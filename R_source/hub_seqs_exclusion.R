#' Deletes top hub sequences from repertoire, changing the network architecture.
#'
#' @param repertoire An annotated AIRR compliant repertoire.
#'
#' (http://docs.airr-community.org/en/latest/)
#' @param top_x Determines what percentage of hub sequences get excluded
#'
#' (Default: 0.005, i.e. Top 0.5 percent)
#' @param report If TRUE saves a csv file containing the excluded sequences.
#' @return Repertoire reduced by hub sequence (new network architecture)
#' @examples
#' \dontrun{
#' repertoire <- immuneSIM(number_of_seqs = 50,species = "mm",receptor = "ig", chain = "h")
#' hub_seqs_exclusion(repertoire, top_x = 0.005,report = FALSE)
#' }


hub_seqs_exclusion<-function(repertoire, top_x = 0.005, report = FALSE){
  ###construct similarity network and mutate hub sequences.
  #only based on AA. --> goal: tilt architecture.
  #library(stringdist,quietly=TRUE,warn.conflicts = FALSE)
  #library(igraph,quietly=TRUE,warn.conflicts = FALSE)

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
  if(report==TRUE){
    report_excluded<-repertoire[repertoire$junction_aa %in% seqs_to_be_excluded,]
    write.csv(report_excluded,file=paste(repertoire$name_repertoire[1],'_hub_seqs_excluded.csv',sep=""))
  }

  return(repertoire_hubs_excluded)

}

