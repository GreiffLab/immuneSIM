#plotting helper function for vdj occurrence
.plot_vdj_occ<-function(curr_repertoire,curr_directory){

  list_order_genes<-list()
  #evaluate VDJ occurrence
  occurrence_v <- as.data.frame(table(curr_repertoire$v_call))
  occurrence_v$Freq <- occurrence_v$Freq/sum(occurrence_v$Freq)
  occurrence_v$gene <- "V"
  list_order_genes[["V"]]<-names(sort(-table(curr_repertoire$v_call)))

  occurrence_d <- as.data.frame(table(curr_repertoire$d_call))
  occurrence_d$Freq <- occurrence_d$Freq/sum(occurrence_d$Freq)
  occurrence_d$gene <- "D"
  list_order_genes[["D"]]<-names(sort(-table(curr_repertoire$d_call)))

  occurrence_j <- as.data.frame(table(curr_repertoire$j_call))
  occurrence_j$Freq <- occurrence_j$Freq/sum(occurrence_j$Freq)
  occurrence_j$gene <- "J"
  list_order_genes[["J"]]<-names(sort(-table(curr_repertoire$j_call)))

  vdj_occurrence <- base::rbind(occurrence_v,occurrence_d,occurrence_j)
  vdj_occurrence$repertoire <- as.character(unique(curr_repertoire$name_repertoire))

  name_rep_vdj<-vdj_occurrence$repertoire[[1]]

  #make ready for plot by ordering V,D,J
  vdj_occurrence_plot_df <- vdj_occurrence
  vdj_occurrence_plot_df$gene <- factor(vdj_occurrence_plot_df$gene,levels=c("V","D","J"))

  #for each of V,D,J create plot
  vdj_occurrence_plot_list<-list()
  for(i in 1:length(unique(vdj_occurrence_plot_df$gene))){
    curr_vdj_occurrence_plot_df <- vdj_occurrence_plot_df[vdj_occurrence_plot_df$gene==unique(vdj_occurrence_plot_df$gene)[i],]

    #in case of J allow for 'random' J gene (used for a particular case)
    if(unique(vdj_occurrence_plot_df$gene)[i]=="J"){
      curr_vdj_occurrence_plot_df$Var1 <- factor(curr_vdj_occurrence_plot_df$Var1,levels=c("IGHJ1","IGHJ2","IGHJ3","IGHJ4","IGHJrandom"))
    }

    curr_vdj_occurrence_plot_df$Var1<-factor(curr_vdj_occurrence_plot_df$Var1,levels=list_order_genes[[unique(vdj_occurrence_plot_df$gene)[i]]])

    #add variable definitions pre plotting.
    Var1 <- Freq <- NULL

    #name_plot <- paste(folder,"/vdj_occurrence_",unique(vdj_occurrence_plot_df$gene)[i],"_",vdj_occurrence$repertoire,"_axis.pdf",sep="")
    vdj_occurrence_plot_list[[i]] <- ggplot2::ggplot(data=curr_vdj_occurrence_plot_df, ggplot2::aes(Var1, Freq))+
      ggplot2::geom_bar(stat="identity",position="dodge") +
      ggplot2::theme_bw()+
      ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
      .theme.akbar()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=7, angle = 90))

  }

  name_plot <- paste("vdj_occurrence_",name_rep_vdj,".pdf",sep="")
  #output plots for V,D,J
  grDevices::pdf(file.path(curr_directory, name_plot),  width = 15,  height = 5)
  grid::pushViewport( grid::viewport(layout=grid::grid.layout(1,3, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(vdj_occurrence_plot_list[[1]], vp=vplayout(1,1))
  base::print(vdj_occurrence_plot_list[[2]], vp=vplayout(1,2))
  base::print(vdj_occurrence_plot_list[[3]], vp=vplayout(1,3))
  grDevices::dev.off()

}
