#plotting function
#comparative plots aafreqs plot for all present lengths and rough VDJ distribution plots
.plot_report_repertoire_AA_freq_comparative<-function(repertoire_seqs,repertoire_name,freqs_input,reference_name="exp_mm_igh"){
  #library(reshape2,quietly=TRUE,warn.conflicts = FALSE)
  #library(Metrics,quietly=TRUE,warn.conflicts = FALSE)
  #library(ggplot2,quietly=TRUE,warn.conflicts = FALSE)

  #set color scheme
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

  #get AA frequency distribution
  AA_freq<-.get_AA_freq_distribution(repertoire_seqs=repertoire_seqs)#repertoire[[1]])

  #find which lengths have information in list of frequencies (lengths that did not occur will have all 0 freqs and need not be used for plotting)
  where_is_info<-which(sapply(1:length(AA_freq),function(x) sum(AA_freq[[x]][[1]])!=0))

  list_aa_plots<-list()
  amino_acid_freqs<-list()
  for(z in 1:length(where_is_info)){

    #calculate AA freq for AA length of interest for first repertoire
    AA_freqs <- AA_freq[[where_is_info[[z]]]]
    curr_length<-names(AA_freq)[where_is_info[[z]]]

    #find current aa freqs for second repertoire
    AA_freqs_input<-freqs_input[[where_is_info[[z]]]]

    #melt both AAfreq into dataframe and turn values into %
    #for rep1
    AA_freq_plot <- reshape2::melt(AA_freqs)
    AA_freq_plot$value <- AA_freq_plot$value*100
    #for rep2
    AA_freqs_input_melt<-reshape2::melt(AA_freqs_input)
    AA_freqs_input_melt$value <- AA_freqs_input_melt$value*100

    #calculate mse btw rep 1 and 2
    mse_v_input<-mean(sapply(1:length(AA_freqs_input_melt$value), function(x) Metrics::mse(AA_freq_plot$value[[x]]/100,AA_freqs_input_melt$value[[x]]/100)))

    #prep plot
    AA_freq_plot$name_repertoire <- repertoire_name
    names(AA_freq_plot) <- c('amino_acid','freqs','position','name_repertoire')
    number_of_positions <- length(AA_freq_plot$amino_acid)/21

    #set position for labels
    AA_freq_plot <- base::transform(AA_freq_plot, mid_y = stats::ave(AA_freq_plot$freqs, AA_freq_plot$position, FUN = function(val) base::cumsum(val) - (0.5 * val)))
    AA_freq_plot$mid_y <- 100-AA_freq_plot$mid_y

    amino_acid_freqs[[z]]<-AA_freqs

    #initialize variables
    amino_acid <- freqs <- mid_y <- position <- NULL

    #plot
    df_length <- ggplot2::ggplot(AA_freq_plot, ggplot2::aes(x=position,y=freqs,fill=amino_acid,label=amino_acid))+
      ggplot2::geom_bar(stat="identity",show.legend = FALSE)+
      ggplot2::geom_text(ggplot2::aes(y = mid_y),colour='black',fontface='bold',show.legend=FALSE)+
      ggplot2::geom_hline(yintercept=c(25,50,75),linetype="dotted")+
      ggplot2::ggtitle(paste("mmse ",repertoire_name," vs ",reference_name,": ",round(mse_v_input,6),sep=""))+
      ggplot2::theme_bw()+
      ggplot2::labs(x = "CDR3 sequence position", y = "Occurrence AA [%]", fill = "", colour = "") +
      ggplot2::scale_fill_manual(values = colorset)+
      ggplot2::scale_colour_manual(values = darkened_colors) +
      ggplot2::scale_x_discrete(limit = c(1:number_of_positions))+
      .theme.akbar()

    list_aa_plots[[z]]<-df_length
    names(list_aa_plots)[z]<-repertoire_name#plot_name
  }

  return(list_aa_plots)
}
