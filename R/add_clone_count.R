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
    xmin <- 1
    number_of_simulated_seqs <- length(simulated_repertoire[,1])
    
    #simulate counts from a powerlaw distribution
    simulated_repertoire$counts <-  poweRlaw::rpldis(number_of_simulated_seqs, xmin, user_defined_alpha)
    #recalculate frequencies from rounded counts
    simulated_repertoire$freqs <-  simulated_repertoire$counts/sum(simulated_repertoire$counts)

  }else{
    #if equal_cc = TRUE make distribution uniform
    simulated_repertoire$freqs<-1/length(simulated_repertoire[,1])
  }

  return(simulated_repertoire)
}
