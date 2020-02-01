#' Plots main repertoire features (length distribution,amino acid frequencies and VDJ usage)
#' @param repertoire An annotated AIRR-compliant immuneSIM repertoire.
#'
#' (http://docs.airr-community.org/en/latest/)
#' @param output_dir String containing full path of desired output folder. If empty figures will be output in tempdir().
#' @param verbose Determines whether messages on plot locations are output to user. (Default: TRUE)
#' @return TRUE (plots saved as pdfs into subfolder 'figures')
#' @examples
#' repertoire <- list_example_repertoires[["example_repertoire_A"]]
#' plot_report_repertoire(repertoire,output_dir="")

#plotting function
#plots aafreqs plot for most common length and rough VDJ distribution plots
plot_report_repertoire<-function(repertoire,output_dir = "", verbose = TRUE){

  #set output folder
  if(output_dir == ""){
    wd <- tempdir()
  }else{
    wd <- output_dir
  }

  #prepare coloring scheme
  darken_color <- function(color, factor=1.4){
    col <- grDevices::col2rgb(color)
    col <- col/factor
    col <- grDevices::rgb(t(col), maxColorValue=255)
  }

  colorset<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','black')
  color_list <- list()
  for(i in 1:length(colorset)){
    color_list[[i]] <- darken_color(colorset[i],factor=1.)
  }
  darkened_colors<- unlist(color_list)

  #set name and evaluate lengths
  name_rep<-repertoire$name_repertoire[[1]]
  lengths_rep<-nchar(as.character(repertoire$sequence_aa))
  lengths_rep_df<-as.data.frame(table(lengths_rep)/sum(table(lengths_rep)))
  names(lengths_rep_df)<-c("length","freq")
  lengths_rep_df$repertoire_name<-name_rep

  maxlen_gen_in<-max(as.integer(as.character(lengths_rep_df$length)))+1
  minlen_gen_in<-min(as.integer(as.character(lengths_rep_df$length)))-1

  length <- freq <- repertoire_name <- NULL
  #plot
  df_seq_length <- ggplot2::ggplot(lengths_rep_df, ggplot2::aes(x=length,y=freq,fill=repertoire_name))+
    ggplot2::geom_bar(stat="identity",position=ggplot2::position_dodge(0.9),show.legend = FALSE)+
    ggplot2::geom_text(data = lengths_rep_df, ggplot2::aes(label = round(freq, digits=2), vjust=-0.7), color='black', position=ggplot2::position_dodge(0.9),fontface='bold',size = 2, show.legend = FALSE) +
    ggplot2::theme_bw()+
    ggplot2::labs(x = "Full VDJ length (a.a.)", y = "Frequency [%]", fill = "", colour = "") +
    ggplot2::scale_fill_manual(values = c("black"))+
    .theme.akbar()

  grDevices::pdf(file.path(wd, paste("length_distribution_",name_rep,".pdf",sep="")), family="Helvetica", width = 8,  height = 6)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(df_seq_length, vp=vplayout(1,1))   # plot for row 1, column 1
  grDevices::dev.off()

  #evaluate amino acid frequencies for most common length
  lengths_cdr3_aa<-nchar(as.character(repertoire$junction_aa))
  lengths_rep_df<-table(lengths_cdr3_aa)/sum(table(lengths_cdr3_aa))
  most_common_length<-names(tail(sort(lengths_rep_df),1))

  repertoire_cdr3_len_14<-as.character(repertoire$junction_aa)[nchar(as.character(repertoire$junction_aa))==most_common_length]

  input_aa_freqs<-.get_AA_freq_distribution(repertoire_cdr3_len_14,length_seqs=most_common_length)

  #get aa freq plot list per entry (parameter combo, incl temp)
  plots_aa_freq_list_imgt<-.plot_report_repertoire_AA_freq(input_aa_freqs=input_aa_freqs,name_rep=name_rep)

  #Plot
  grDevices::pdf(file.path(wd, paste("aa_freq_",name_rep,".pdf",sep="")), family="Helvetica", width = 13,  height = 10)
  grid::pushViewport( grid::viewport(layout = grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(plots_aa_freq_list_imgt[[1]], vp=vplayout(1,1))   # plot for row 1, column 1
  grDevices::dev.off()

  #plot vdj occurrence
  .plot_vdj_occ(curr_repertoire=repertoire,curr_directory=wd)

  if(verbose==TRUE){
    cat("Plots saved in PDF format in: ",wd,"\n")
  }

}

