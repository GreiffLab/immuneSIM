#' Loads full insertion/deletion data from GitHub
#' @return Dataframe containing insertions and deletions (11363603 rows, 6 columns)
#' @examples
#' \donttest{
#' full_insertions_and_deletion_df <- load_insdel_data()
#' }

load_insdel_data <- function(){
  #initialize variable
  insdel_data <- NULL
  #load source data from github
  repmis::source_data("https://github.com/GreiffLab/immuneSIM/raw/master/additional_data/full_insertions_and_deletion_data.Rdata?raw=true")

  return(insdel_data)
}
