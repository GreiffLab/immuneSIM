#' Generates a dataframe from separate heavy and light or beta and alpha chain dataframes
#'
#' @param repertoire_heavy A repertoire containing heavy/beta chain data
#' @param repertoire_light A repertoire containing light/alpha chain data
#' @return immuneSIM repertoire containing heavy/beta and light/alpha chain data.
#' @examples
#' repertoire_heavy <- immuneSIM(number_of_seqs = 5,species = "mm",receptor = "ig", chain = "h")
#' repertoire_light <- immuneSIM(number_of_seqs = 5,species = "mm",receptor = "ig", chain = "kl")
#' paired_repertoire <- combine_into_paired(repertoire_heavy,repertoire_light)

combine_into_paired <- function(repertoire_heavy,repertoire_light){
    repertoire_light <- repertoire_light[,c(1:5,7:10,13,14,16)]
    names(repertoire_light) <- c("sequence_L", "sequence_aa_L", "junction_L", "junction_aa_L", "v_call_L", "j_call_L", "np1_L", "np2_L", "del_v_L","del_j_L", "v_sequence_alignment_L", "j_sequence_alignment_L")

    repertoire_paired <- cbind(repertoire_heavy,repertoire_light)

    return(repertoire_paired)
}
