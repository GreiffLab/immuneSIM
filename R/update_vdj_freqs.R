.update_vdj_freqs<-function(vdj_list,VDJ_used){
  #since not every vdj gets picked at the desired rate (because non fct sims are discarded)
  #we need to adjust vdj freqs to reflect changed probabilities after a set number of sims

  #generate df of VDJs and the frequencies they were used at
  vdj_used_datatable<-data.table::rbindlist(VDJ_used)
  vdj_used_df<-data.table::setDF(vdj_used_datatable)

  v_factored<-factor(vdj_used_df$v_gene,levels=vdj_list$V$gene)
  d_factored<-factor(vdj_used_df$d_gene,levels=vdj_list$D$gene)
  j_factored<-factor(vdj_used_df$j_gene,levels=vdj_list$J$gene)

  #generate dataframe for comparison and set the ones that were at 0 to the minimal achieved frequency (to avoid dividing by zero)
  V_used<-reshape2::melt(table(v_factored)/sum(table(v_factored)))
  V_used[V_used$value==0.00,'value']<-min(V_used[V_used$value!=0.00,'value'])
  D_used<-reshape2::melt(table(d_factored)/sum(table(d_factored)))
  D_used[D_used$value==0.00,'value']<-min(D_used[D_used$value!=0.00,'value'])
  J_used<-reshape2::melt(table(j_factored)/sum(table(j_factored)))
  J_used[J_used$value==0.00,'value']<-min(J_used[J_used$value!=0.00,'value'])

  #get initial frequencies and save into new list
  vdj_list_new<-vdj_list

  #for each of V,D,J scale initial frequency by how far off the freqs currently are
  vdj_list_new$V$frequency<-vdj_list_new$V$frequency*(vdj_list_new$V$frequency/(V_used$value))
  vdj_list_new$V$frequency<-vdj_list_new$V$frequency/sum(vdj_list_new$V$frequency)

  vdj_list_new$D$frequency<-vdj_list_new$D$frequency*(vdj_list_new$D$frequency/(D_used$value))
  vdj_list_new$D$frequency<-vdj_list_new$D$frequency/sum(vdj_list_new$D$frequency)

  vdj_list_new$J$frequency<-(vdj_list_new$J$frequency*(vdj_list_new$J$frequency/(J_used$value)))
  vdj_list_new$J$frequency<-vdj_list_new$J$frequency/sum(vdj_list_new$J$frequency)

  return(vdj_list_new)

}
