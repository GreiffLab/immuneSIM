#' Comparative plots of main repertoire features of two input repertoires (length distribution, amino acid frequency, VDJ usage, kmer occurrence)
#' @param repertoire_A An annotated AIRR-compliant immuneSIM repertoire.
#'
#' (http://docs.airr-community.org/en/latest/)
#' @param repertoire_B An annotated AIRR-compliant immuneSIM repertoire.
#' @param names_repertoires A vector containing two strings denoting the names of the repertoires / repertoire descriptions.
#' @param length_aa_plot Defines sequence length for which the amino acid frequency plot will be made.
#' @param output_dir String containing full path of desired output folder. If empty, figures will be output in tempdir().
#' @param verbose Determines whether messages on plot locations are output to user. (Default: TRUE)
#' @return TRUE (plots saved as pdfs into subfolder 'figures')
#' @examples
#' repertoire_A <- list_example_repertoires[["example_repertoire_A"]]
#' repertoire_B <- list_example_repertoires[["example_repertoire_B"]]
#' plot_repertoire_A_vs_B(
#' repertoire_A,
#' repertoire_B,
#' c("Sim_repertoire_1","Sim_repertoire_2"),
#' length_aa_plot = 14,
#' output_dir="")


plot_repertoire_A_vs_B<-function(repertoire_A,repertoire_B,names_repertoires=c("Repertoire_A","Repertoire_B"),length_aa_plot=14,output_dir="", verbose=TRUE){

  #set colorpalette
  my_spectral <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,'Spectral'))(14)

  #set output folder
  if(output_dir == ""){
    wd <- tempdir()
  }else{
    wd <- output_dir
  }

  #set names
  name_a<-names_repertoires[1]
  name_b<-names_repertoires[2]

  name<-paste(name_a, "_vs_", name_b,sep="")

  #get sequences for amino acid freq plotting
  list_plot_dfs<-list()

  #make sure all is character type
  repertoire_A_cdr3<-as.character(repertoire_A$junction_aa)
  repertoire_A_cdr3_nt<-as.character(repertoire_A$junction)

  repertoire_B_jct_aa<-as.character(repertoire_B$junction_aa)
  repertoire_B_jct_nt<-as.character(repertoire_B$junction)

  #evaluate lengths
  lengths_rep_A<-nchar(repertoire_A_cdr3)
  lengths_rep_B<-nchar(repertoire_B_jct_aa)

  #choose sequences of length defined by user
  repertoire_A_cdr3_len_14<-repertoire_A_cdr3[nchar(repertoire_A_cdr3)==length_aa_plot]
  repertoire_B_cdr3_len_14<-repertoire_B_jct_aa[nchar(repertoire_B_jct_aa)==length_aa_plot]

  #if there are sequences of this length in both create amino acid occurrence plot
  if((length(repertoire_B_cdr3_len_14)>0) & (length(repertoire_A_cdr3_len_14) > 0)){
    #subsample the larger one so there is an equal number of seqs underlying each plot
    if(length(repertoire_B_cdr3_len_14) > length(repertoire_A_cdr3_len_14)){
      repertoire_B_cdr3_len_14<-sample(repertoire_B_cdr3_len_14,length(repertoire_A_cdr3_len_14))
    }else if(length(repertoire_B_cdr3_len_14) < length(repertoire_A_cdr3_len_14)){
      repertoire_A_cdr3_len_14<-sample(repertoire_A_cdr3_len_14,length(repertoire_B_cdr3_len_14))
    }

    #get amino acid frequencies
    input_aa_freqs<-.get_AA_freq_distribution(as.character(repertoire_B_cdr3_len_14))

    #get aa freq plot list per entry (parameter combo, incl temp)
    plots_aa_freq_list_imgt<-.plot_report_repertoire_AA_freq_comparative(repertoire_seqs=repertoire_A_cdr3_len_14,
                                                           repertoire_name=name_a,freqs_input=input_aa_freqs,reference_name=name_b)

    plots_aa_freq_list_cellreports_np2<-.plot_report_repertoire_AA_freq_comparative(repertoire_seqs=repertoire_B_cdr3_len_14,
                                                                      repertoire_name=name_b,freqs_input=input_aa_freqs,reference_name=name_b)

    #output plot
    grDevices::pdf(file.path(wd, paste("aa_freq_",name,".pdf",sep="")), family="Helvetica", width = 20,  height = 10)
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1,2, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
    vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
    base::print(plots_aa_freq_list_imgt[[1]], vp=vplayout(1,1))   # plot for row 1, column 1
    base::print(plots_aa_freq_list_cellreports_np2[[1]], vp=vplayout(1,2))   # plot for row 1, column 1
    grDevices::dev.off()
  }else{
    #if there are not enough sequences to evaluate notify user
    base::warning("\nNot enough sequences of length ", length_aa_plot ," to make amino acid frequency plot\n")
  }



  #vdj recovery plot
  #merge the two vdj occurrences
  vdj_occurrence_df_rep_A <- .calc_vdj_occ(repertoire_A,name_a)#"Repertoire_A")
  vdj_occurrence_df_rep_A <-vdj_occurrence_df_rep_A[,c("Var1","Freq","gene","repertoire")]
  vdj_occurrence_df_rep_B <- .calc_vdj_occ(repertoire_B,name_b)#"Repertoire_B")

  vdj_occurrence_df_rep_A_B<-base::merge(vdj_occurrence_df_rep_B,vdj_occurrence_df_rep_A,by="Var1",all=TRUE)
  #fill in NAs that were created. first freqs
  vdj_occurrence_df_rep_A_B[is.na(vdj_occurrence_df_rep_A_B$Freq.y),"Freq.y"]<-0
  vdj_occurrence_df_rep_A_B[is.na(vdj_occurrence_df_rep_A_B$Freq.x),"Freq.x"]<-0
  #then repertoire names
  vdj_occurrence_df_rep_A_B[is.na(vdj_occurrence_df_rep_A_B$repertoire.x),"repertoire.x"]<-vdj_occurrence_df_rep_A_B[!is.na(vdj_occurrence_df_rep_A_B$repertoire.x),"repertoire.x"][[1]]
  vdj_occurrence_df_rep_A_B[is.na(vdj_occurrence_df_rep_A_B$repertoire.y),"repertoire.y"]<-vdj_occurrence_df_rep_A_B[!is.na(vdj_occurrence_df_rep_A_B$repertoire.y),"repertoire.y"][[1]]
  #finally make a new gene column that brings together all the info
  #turn into factor and delete helper columns
  vdj_occurrence_df_rep_A_B$gene.x[is.na(vdj_occurrence_df_rep_A_B$gene.x)]<-vdj_occurrence_df_rep_A_B$gene.y[is.na(vdj_occurrence_df_rep_A_B$gene.x)]
  vdj_occurrence_df_rep_A_B$gene<-factor(vdj_occurrence_df_rep_A_B$gene.x,levels=c("V","D","J"))
  vdj_occurrence_df_rep_A_B$gene.y<-NULL
  vdj_occurrence_df_rep_A_B$gene.x<-NULL

  vdj_occurrence_df_rep_A_B$r_pearson<-0
  vdj_occurrence_df_rep_A_B$r_spearman<-0

  #for each of V,D and J evaluate pearson,spearman correlation
  for(i in 1:length(levels(vdj_occurrence_df_rep_A_B$gene))){
    current_level <- levels(vdj_occurrence_df_rep_A_B$gene)[i]
    curr_vdj_occurrence_plot_df <- vdj_occurrence_df_rep_A_B[vdj_occurrence_df_rep_A_B$gene==current_level,]

    r_value_spearman<-stats::cor(curr_vdj_occurrence_plot_df$Freq.x,curr_vdj_occurrence_plot_df$Freq.y,method="spearman")
    r_value_pearson<-stats::cor(curr_vdj_occurrence_plot_df$Freq.x,curr_vdj_occurrence_plot_df$Freq.y,method="pearson")

    vdj_occurrence_df_rep_A_B[vdj_occurrence_df_rep_A_B$gene==current_level,"r_spearman"]<-r_value_spearman
    vdj_occurrence_df_rep_A_B[vdj_occurrence_df_rep_A_B$gene==current_level,"r_pearson"]<-r_value_pearson

  }

  #prepare and then plot
  name_plot <- paste("vdj_occurrence_",unique(vdj_occurrence_df_rep_A_B$gene)[i],"_axis.pdf",sep="")
  scale <- max(1.5*max(vdj_occurrence_df_rep_A_B$Freq.x),1.5*max(vdj_occurrence_df_rep_A_B$Freq.y))
  pos_text_x <- .75*scale
  pos_text_y <- .15*scale

  Freq.y <- Freq.x <- r_pearson <-  r_spearman <- NULL

  vdj_occurrence_plot <- ggplot2::ggplot(data=vdj_occurrence_df_rep_A_B, ggplot2::aes(Freq.y, Freq.x))+
    ggplot2::geom_point(alpha = 0.7, size = 3, color=my_spectral[3]) +
    ggplot2::facet_grid(.~gene,scales='free')+
    ggplot2::theme_bw()+
    ggplot2::geom_abline()+
    ggplot2::geom_density_2d(colour = my_spectral[3])+
    ggplot2::scale_x_continuous(limits=c(0,scale),breaks=seq(0,1,0.1))+
    ggplot2::scale_y_continuous(limits=c(0,scale),breaks=seq(0,1,0.1))+
    ggplot2::geom_text(ggplot2::aes(pos_text_x, pos_text_y, label=paste("r_pearson = ",round(r_pearson,3),sep="")))+
    ggplot2::geom_text(ggplot2::aes(pos_text_x, .7*pos_text_y, label=paste("r_spearman = ",round(r_spearman,3),sep="")))+
    ggplot2::labs(x = name_a, y = name_b, fill = "") +
    .theme.akbar()

  grDevices::pdf(file.path(wd, paste("vdj_recovery_",name,".pdf",sep="")), family="Helvetica", width = 10,  height = 5)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(vdj_occurrence_plot, vp=vplayout(1,1))   # plot for row 1, column 1
  grDevices::dev.off()


  #length distribution plot
  #evaluate length os sequences and the frequencies with which they occur.
  lengths_A_df<-as.data.frame(table(lengths_rep_A)/sum(table(lengths_rep_A)))
  lengths_B_df<-as.data.frame(table(lengths_rep_B)/sum(table(lengths_rep_B)))
  names(lengths_A_df)<-c("length","freq")
  lengths_A_df$repertoire<-name_a
  names(lengths_B_df)<-c("length","freq")
  lengths_B_df$repertoire<-name_b

  lengths_df<-base::rbind(lengths_A_df,lengths_B_df)

  #find min and max lenghts for plot bounding
  maxlen_gen_in<-max(as.integer(as.character(lengths_df$length)))+1
  minlen_gen_in<-min(as.integer(as.character(lengths_df$length)))-1

  #order sequence lenghts numerically along x-axis
  lengths_df$repertoire_fill <- factor(lengths_df$repertoire,levels=names_repertoires)
  lengths_df$length<-factor(lengths_df$length,levels=sort(as.integer(levels(lengths_df$length))))#<-factor(lengths_df$length,levels=as.numeric(lengths_df$length))

  length <- freq <- repertoire_fill <- NULL

  df_seq_length <- ggplot2::ggplot(lengths_df, ggplot2::aes(x=length,y=freq,fill=repertoire_fill))+
    ggplot2::geom_bar(stat="identity",position=ggplot2::position_dodge(0.9),show.legend = TRUE)+
    ggplot2::geom_text(data = lengths_df, ggplot2::aes(label = round(freq, digits=2), vjust=-0.7), color='black', position=ggplot2::position_dodge(0.9),fontface='bold',size = 1.5, show.legend = FALSE) +
    ggplot2::theme_bw()+
    ggplot2::labs(x = "CDR3 sequence length", y = "Frequency [%]", fill = "", colour = "") +
    ggplot2::scale_fill_manual(values = c("black","red"))+
    .theme.akbar()

  grDevices::pdf(file.path(wd, paste("length_distribution_",name,".pdf",sep="")), family="Helvetica", width = 8,  height = 6)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(df_seq_length, vp=vplayout(1,1))   # plot for row 1, column 1
  grDevices::dev.off()


  #kmer plots
  #if larger then subsample to 10000 sequences (to keep runtime in check)
  if(length(repertoire_B_jct_nt)>10000){
    repertoire_B_jct_nt<-sample(repertoire_B_jct_nt,10000)
  }
  if(length(repertoire_A_cdr3_nt)>10000){
    repertoire_A_cdr3_nt<-sample(repertoire_A_cdr3_nt,10000)
  }

  #Get kmer dictionary
  dictionary_counts<-.make_kmer_dictionary(k_nt=3,k_aa=1,gap_size_nt=c(0,1,2,3),gap_size_aa=c(0,1,2))
  k_mers_in_dict <- as.integer(names(dictionary_counts[["nt"]]))
  gap_sizes <- dictionary_counts[["gap_nt"]]

  #evaluate kmer occurrence
  kmer_input<-.kmer_occurrence_opt_k_m_gen(sequences_to_kmerize=repertoire_B_jct_nt,AA_nt="nt",opt_k=k_mers_in_dict,opt_m=gap_sizes,dictionary_ref=dictionary_counts)
  kmer_gen<-.kmer_occurrence_opt_k_m_gen(sequences_to_kmerize=repertoire_A_cdr3_nt,AA_nt="nt",opt_k=k_mers_in_dict,opt_m=gap_sizes,dictionary_ref=dictionary_counts)

  #name kmer distributions
  kmer_rep_1<-kmer_gen
  name_rep_1<-name_a#"Repertoire_A"
  kmer_rep_2<-kmer_input
  name_rep_2<-name_b#"Repertoire_B"

  #prep plot dataframe
  params_df_all <- data.frame(pair=kmer_rep_1$pairs,occurrence_rep_1=as.vector(100*kmer_rep_1[["freq_counts"]]),occurrence_rep_2=as.vector(100*kmer_rep_2[["freq_counts"]]),rep_1=name_rep_1,rep_2=name_rep_2)
  params_df_all$kmer_length <- factor(kmer_rep_1$kmer_length,levels=unique(kmer_rep_1$kmer_length))
  params_df_all$kmer_dist <- factor(kmer_rep_1$kmer_dist,levels=unique(kmer_rep_1$kmer_dist))

  #delete all rows that have 0 input AND 0 generated to clean up plot (make it easier to open)
  params_df_all_no_zeros <- params_df_all[!(params_df_all$occurrence_rep_1==0 & params_df_all$occurrence_rep_2==0),]

  #give new names for naming column
  new_denom_with_accurate_gap<-sapply(1:length(params_df_all_no_zeros$pair),function(x) paste(base::strsplit(as.character(params_df_all_no_zeros$pair[[x]]),",")[[1]][1],paste(rep(".",as.integer(as.character(params_df_all_no_zeros$kmer_dist[[x]]))),collapse=""),base::strsplit(as.character(params_df_all_no_zeros$pair[[x]]),",")[[1]][2],sep=""))
  params_df_all_no_zeros$pair_in_one<-new_denom_with_accurate_gap

  #evaluate spearman an pearson correlation
  r_value <- stats::cor(params_df_all_no_zeros$occurrence_rep_1,params_df_all_no_zeros$occurrence_rep_2,method='pearson')
  r_value_s <- stats::cor(params_df_all_no_zeros$occurrence_rep_1,params_df_all_no_zeros$occurrence_rep_2,method='spearman')

  params_df_all_no_zeros$r_value <- paste("r_pearson = ",round(r_value,2),sep="")
  params_df_all_no_zeros$r_value_spearman <- paste("r_spearman = ",round(r_value_s,2),sep="")

  #set plotting parameters
  scale <- max(1.5*max(params_df_all$occurrence_rep_1),1.5*max(params_df_all$occurrence_rep_2))
  params_df_all_no_zeros$diff<-abs(params_df_all_no_zeros$occurrence_rep_1-params_df_all_no_zeros$occurrence_rep_2)
  top_diffs<-utils::tail(sort(params_df_all_no_zeros$diff),5*length(unique(kmer_rep_1$kmer_dist)))
  params_df_all_no_zeros$top_diff<-ifelse(params_df_all_no_zeros$diff %in% top_diffs,1,0)

  pos_text_x <- .75*scale
  pos_text_y <- .15*scale

  occurrence_rep_1 <- occurrence_rep_2 <- pair_in_one <- NULL

  #plot
  kmer_plot_sim_v_input<-ggplot2::ggplot(params_df_all_no_zeros, ggplot2::aes(x=occurrence_rep_1,y=occurrence_rep_2,label=pair_in_one)) +
    ggplot2::geom_point(alpha = 0.7, size = 3, color=my_spectral[3])+
    ggplot2::theme_bw()+
    ggplot2::geom_density_2d(colour = my_spectral[3])+
    ggplot2::labs(x = unique(params_df_all_no_zeros$rep_1), y =  unique(params_df_all_no_zeros$rep_2), fill = "", colour = "") +
    ggplot2::scale_x_continuous(limits=c(0,scale),breaks=seq(0,100,1))+
    ggplot2::scale_y_continuous(limits=c(0,scale),breaks=seq(0,100,1))+
    ggplot2::geom_abline()+
    ggplot2::geom_text(ggplot2::aes(pos_text_x, pos_text_y, label=r_value))+
    ggplot2::geom_text(ggplot2::aes(pos_text_x, 0.7*pos_text_y, label=r_value_spearman))+
    .theme.akbar()

  grDevices::pdf(file.path(wd, paste("kmer_plot_",name,".pdf",sep="")), family="Helvetica", width = 8,  height = 7)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  #grid.text("", vp=vplayout(1,1:2))
  base::print(kmer_plot_sim_v_input, vp=vplayout(1,1))   # plot for row 1, column 1
  grDevices::dev.off()

  if(verbose==TRUE){
    cat("Plots were saved in PDF format in: ",wd,"\n")
  }

}
