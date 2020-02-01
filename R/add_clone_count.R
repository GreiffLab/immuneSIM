# Add clone count.
#
# @param simulated_repertoire An annotated AIRR (Adaptive Immune Receptor Repertoire) compliant immuneSIM repertoire.
# @param user_defined_alpha Parameter that for the simulation of powerlaw distribution (recommended range 2-5)
# @param equal_cc overrides other parameters if TRUE and clone count will be equal for all sequences
# @return immuneSIM repertoire including clone count and frequency column.
# @examples
# add_clone_count(repertoire, 2, FALSE)

.add_clone_count<-function(simulated_repertoire,user_defined_alpha,equal_cc=FALSE){
  #library(poweRlaw,quietly=TRUE,warn.conflicts = FALSE)

  if(equal_cc==FALSE){
    #set parameters
    xmin<-1
    alpha_value<-user_defined_alpha
    number_of_simulated_seqs<-length(simulated_repertoire[,1])#

    #number of entries required
    x<-xmin:number_of_simulated_seqs #CDR3 rank (= 100 unique CDR3s)#

    #get powerlaw distribution
    freq_counts<-poweRlaw::dpldis(x,xmin,alpha_value)
    freq_counts <- freq_counts*(1/sum(freq_counts))
    #randomize order of frequencies
    simulated_repertoire$freqs<-sample(freq_counts)

    #calculate counts (rounded from frequency distribution)
    simulated_repertoire$counts <-  round(freq_counts*(max(freq_counts)/min(freq_counts))/max(freq_counts))
    #recalculate frequencies from rounded counts
    simulated_repertoire$freqs <-  simulated_repertoire$counts/sum(simulated_repertoire$counts)

  }else{
    #if equal_cc = TRUE make distribution uniform
    simulated_repertoire$freqs<-1/length(simulated_repertoire[,1])
  }

  return(simulated_repertoire)
}
