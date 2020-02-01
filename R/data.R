#' Translation dictionary amino acid <-> nucleotide codon
#'
#' A dataframe containing a mapping from each of 64 codons
#' to amino acids.
#'
#' @format A data frame with 64 rows and  variables:
#' \describe{
#'   \item{aa}{amino acid}
#'   \item{codon}{nucleotide codon}
#' }
#' @source \url{https://www.genscript.com/tools/codon-table}
"gen_code"



#' Hotspot dataframe for SHM
#'
#' A dataframe containing mutation probabilities for every possible 5mer pattern
#'
#' @format A data frame with 1024 rows and  variables:
#' \describe{
#'   \item{pattern}{amino acid}
#'   \item{toA}{probability of mutation to adenine}
#'   \item{toC}{probability of mutation to cytosine}
#'   \item{toG}{probability of mutation to guanine}
#'   \item{toT}{probability of mutation to thymine}
#'   \item{Source}{source of probability}
#' }
#' @source \url{https://cran.r-project.org/package=AbSim}
"hotspot_df"




#' Dataframe containing insertion sequences and deletion lengths
#'
#' A dataframe containing all insertions and deletions
#' observed in experimental data (pooled across all samples, Greiff, 2017)
#' This dataframe is a subset of the dataframe used in the application note.
#' The original dataframe which contains 11363603 rows can be downloaded
#' from:
#'
#' https://github.com/GreiffLab/immuneSIM or using the provided
#' function: load_insdel_data()
#' @format A data frame with 500000 rows and variables:
#' \describe{
#'   \item{n1}{np1 insertions}
#'   \item{n2}{np2 insertions}
#'   \item{del_v}{lengths of V gene deletions}
#'   \item{del_d_5}{lengths of 5' end D gene deletions}
#'   \item{del_d_3}{lengths of 3' end D gene deletions}
#'   \item{del_j}{lengths of J gene deletions}
#' }
#' @source \url{https://doi.org/10.1016/j.celrep.2017.04.054}
"insertions_and_deletion_lengths_df"


#' Vector containing VDJ length distributions
#'
#' A vector containing 10000 VDJ lengths for simulating
#' of fully random sequences (independent of germline genes)
#'
#' @format A vector with 10000 entries:
#' \describe{
#'   \item{length}{VDJ nucleotide lengths sampled from murine naive
#'   follicular B-cell data, Greiff 2017}
#' }
#' @source \url{https://doi.org/10.1016/j.celrep.2017.04.054}
"length_dist_simulation"



#' Collection of germline genes and frequencies
#'
#' A list containing sublists for species ("hs","mm") which in turn
#' contain sublists for receptors ("ig","tr") which are subset in
#' chains ("h", "k", "l" and "b", "a", respectively). Each entry
#' contains a list of three dataframes ("V","D" and "J") with the major IMGT
#' annotated germline genes including name, sequence based on IMGT and frequencies based on
#' experimental data from DeWitt(2017), Emerson (2017), Greiff (2017) and Madi (2017)
#'
#' @format A list of lists containing dataframes with up to 126 entries:
#' \describe{
#'   \item{gene}{name of germline gene}
#'   \item{allele}{allele number (presently restricted to allele 01)}
#'   \item{sequence}{nucleotide sequence of germline gene}
#'   \item{species}{name of species}
#'   \item{frequency}{Frequencies of germline genes based on experimental data}
#' }
#' @source {
#'    \url{http://www.imgt.org/vquest/refseqh.html}
#'
#'    \url{https://doi.org/10.1371/journal.pone.0160853}
#'
#'    \url{https://doi.org/10.1038/ng.3822}
#'
#'    \url{https://doi.org/10.1016/j.celrep.2017.04.054}
#'
#'    \url{https://doi.org/10.7554/eLife.22057}
#'    }
"list_germline_genes_allele_01"


#' One Spot
#'
#' A dataframe containing a mutation probabilities to base per 5mer
#' (inherited from AbSim package)
#'
#' @format A dataframe with 32 entries:
#' \describe{
#'   \item{pattern}{amino acid}
#'   \item{toA}{probability of mutation to adenine}
#'   \item{toC}{probability of mutation to cytosine}
#'   \item{toG}{probability of mutation to guanine}
#'   \item{toT}{probability of mutation to thymine}
#'   \item{Source}{source of probability}
#' }
#' @source {
#'    \url{https://cran.r-project.org/package=AbSim}
#'
#'    \url{https://doi.org/10.1093/bioinformatics/btx533}
#' }
"one_spot_df"



#' Example repertoires
#'
#' A list containing two example repertoires (100 sequences each)
#' simulated with immuneSIM using default parameters.
#' These repertoires are used in the examples.
#'
#' @format A list with 2 entries:
#' \describe{
#'   \item{example_repertoire_A}{Repertoire simulated using standard parameters (A)}
#'   \item{example_repertoire_A}{Repertoire simulated using standard parameters (B)}
#' }
#' @source {
#'    \url{https://immunesim.readthedocs.io}
#' }
"list_example_repertoires"
