#' Replaces codons with synonymous codons
#' @param repertoire An annotated AIRR compliant immuneSIM repertoire.
#'
#' (http://docs.airr-community.org/en/latest/)
#' @param mode Defines whether codons should be replaced in the nt or AA sequence or in both ("nt","AA","both")
#' @param codon_replacement_list List containing instructions for which codons should be replaced and how
#' @param skip_probability Probability with which a sequence gets skipped in the codon replacement process between 0,1
#' @return immuneSIM repertoire with replaced codons
#' @examples
#' repertoire <- list_example_repertoires[["example_repertoire_A"]]
#' rep_codon_repl <- codon_replacement(repertoire, "both",
#' list(tat = "tac", agt = "agc", gtt = "gtg"), 0)


codon_replacement<-function(repertoire,mode="both",codon_replacement_list,skip_probability=0){
  ###replace codons either with or without change in AA
  ### for control of kmer freqs (through synyonym replacements)
  ## and/or pattern injection through a combo of synyonym and nonsynonym replacements
  ## needs replacement list. named list where every entry is name=current codon, entry=replacement codon
  ## skip probabilty allows control of how likely a row gets skipped.
  #if synonym. only change on nt level goal: tilt kmer occurrence.
  #library(stringr,quietly=TRUE,warn.conflicts = FALSE)
  #library(Biostrings,quietly=TRUE,warn.conflicts = FALSE)

  repertoire$sequence<-as.character(repertoire$sequence)
  repertoire$sequence_aa<-as.character(repertoire$sequence_aa)
  repertoire$junction<-as.character(repertoire$junction)
  repertoire$junction_aa<-as.character(repertoire$junction_aa)

  #prepare dataframe for replacement event recording
  repertoire$codon_replacement<-""
  replacements_df<-reshape2::melt(codon_replacement_list)
  names(replacements_df)<-c("replacement","initial")
  replacements_df$count<-0
  replacements_df$replacement<-as.character(replacements_df$replacement)
  replacements_df$initial<-as.character(replacements_df$initial)

  #string to add if no replacement occurs
  replacements_df_no_cr <- paste(sapply(1:nrow(replacements_df),function(x) paste(paste(replacements_df[x,"initial"],replacements_df[x,"replacement"],sep=","),replacements_df[x,"count"],sep=":")),collapse="|")

  if(mode %in% c("both","AA","aa","Aa")){
    for(i in 1:nrow(repertoire)){

      #evaluate whether current sequence should undergo codon replacement or be skipped
      skip_yes_no<-sample(c("yes","no"), 1, prob=c(skip_probability,1-skip_probability))

      if(skip_yes_no=="no"){

        #prep replacement in vdj and cdr3 for sequence (split seqs in codons)
        vdj_nt<-as.character(repertoire$sequence[[i]])
        start = seq(1, nchar(vdj_nt), 3)
        stop  = pmin(start + 2, nchar(vdj_nt))
        vdj_split<-stringr::str_sub(vdj_nt, start, stop)

        cdr3_nt<-as.character(repertoire$junction[[i]])
        start_cdr3 = seq(1, nchar(cdr3_nt), 3)
        stop_cdr3  = pmin(start_cdr3 + 2, nchar(cdr3_nt))
        cdr3_split<-stringr::str_sub(cdr3_nt, start_cdr3, stop_cdr3)

        #replacement in vdj and cdr3
        for(j in 1:length(codon_replacement_list)){
          #get replacement codon
          vdj_replacement<-names(codon_replacement_list)[j]

          #record number of occurrences of current replacement
          nb_of_occurrences<-length(vdj_split[vdj_split==vdj_replacement])
          replacements_df$count[j]<-nb_of_occurrences

          #replace codons in full VDJ
          vdj_split[vdj_split==vdj_replacement]<-codon_replacement_list[[vdj_replacement]]#"agc"

          #do the same for CDR3 (l.73 for naming clarity)
          cdr3_replacement<-names(codon_replacement_list)[j]
          cdr3_split[cdr3_split==cdr3_replacement]<-codon_replacement_list[[cdr3_replacement]]
        }

        #record replacement event
        repertoire$codon_replacement[i]<-paste(sapply(1:nrow(replacements_df),function(x) paste(paste(replacements_df[x,"initial"],replacements_df[x,"replacement"],sep=","),replacements_df[x,"count"],sep=":")),collapse="|")

        #collapse seqs again and translate into a.a.
        vdj_new_nt<-paste(vdj_split,collapse="")
        vdj_new_aa<-as.character(Biostrings::translate(Biostrings::DNAStringSet(vdj_new_nt)))

        cdr3_new_nt<-paste(cdr3_split,collapse="")
        cdr3_new_aa<-as.character(Biostrings::translate(Biostrings::DNAStringSet(cdr3_new_nt)))

        #replace existing sequences depending on mode (either only a.a., only nt or both sequences)
        if(mode %in% c("both","AA","Aa","aa")){
          repertoire$sequence_aa[[i]]<-vdj_new_aa
          repertoire$junction_aa[[i]]<-cdr3_new_aa
        }

        if(mode %in% c("both","nt","Nt","NT")){
          repertoire$junction[[i]]<-cdr3_new_nt
          repertoire$sequence[[i]]<-vdj_new_nt
        }
      }else{
        #if sequence is skipped record event as as zero replacements occured
        repertoire$codon_replacement[i]<-replacements_df_no_cr
      }

    }
  }

  return(repertoire)
}
