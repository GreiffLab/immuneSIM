# Drops out VDJ genes from germline gene list
#
# @param vdj_list List containing germline gene frequencies
# @param dropout_vdj named vector with numeric entry for V,D,J that denotes number of germline genes that should be dropped out
# @return VDJ list with modified germline gene frequencies.
# @examples
# dropout_VDJ(vdj_list, c(V=1,D=1,J=1))

.dropout_VDJ<-function(vdj_list,dropout_vdj){
  #allows dropout of V,D,J genes
  #specific or number
  #check whether V-genes should be dropped out
  if(dropout_vdj[["V"]]!=0){
    #sample V-genes to dropout
    number_of_genes_to_drop <- dropout_vdj[["V"]]
    if(number_of_genes_to_drop < nrow(vdj_list$V)){
      rows_to_drop <- sample(1:nrow(vdj_list$V),number_of_genes_to_drop)
    }else{
      #if number to drop out is larger than germline diversity dropout all but one gene.
      rows_to_drop <- sample(1:nrow(vdj_list$V),nrow(vdj_list$V)-1)
    }
    #set frequency of chosen dropouts to zero and adjust overall frequency (to sum up to 1)
    vdj_list$V[rows_to_drop,"frequency"] <- 0
    vdj_list$V$frequency <- vdj_list$V$frequency*(1/sum(vdj_list$V$frequency))
  }

  #repeat process for D and J gene
  if(dropout_vdj[["D"]]!=0){
    number_of_genes_to_drop <- dropout_vdj[["D"]]
    if(number_of_genes_to_drop < nrow(vdj_list$D)){
      rows_to_drop <- sample(1:nrow(vdj_list$D),number_of_genes_to_drop)
    }else{
      rows_to_drop <- sample(1:nrow(vdj_list$D),nrow(vdj_list$D)-1)
    }
    vdj_list$D[rows_to_drop,"frequency"] <- 0
    vdj_list$D$frequency <- vdj_list$D$frequency*(1/sum(vdj_list$D$frequency))
  }
  if(dropout_vdj[["J"]]!=0){
    number_of_genes_to_drop <- dropout_vdj[["J"]]
    if(number_of_genes_to_drop < nrow(vdj_list$J)){
      rows_to_drop <- sample(1:nrow(vdj_list$J),number_of_genes_to_drop)
    }else{
      rows_to_drop <- sample(1:nrow(vdj_list$J),nrow(vdj_list$J)-1)
    }
    vdj_list$J[rows_to_drop,"frequency"] <- 0
    vdj_list$J$frequency <- vdj_list$J$frequency*(1/sum(vdj_list$J$frequency))
  }

  return(vdj_list)
}
