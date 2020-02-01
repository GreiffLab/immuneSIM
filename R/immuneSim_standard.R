
.immuneSim_standard<-function(number_of_seqs,vdj_list,species="mm",receptor="ig",chain="h",insertions_and_deletion_lengths,user_defined_alpha=2,name_repertoire,freq_update_time=0.5*number_of_seqs,max_cdr3_length=100,min_cdr3_length=6,shm.mode="none",shm.prob=15/350,verbose=TRUE){
  #based on input VDJ list, insertion and deletion frequencies and clonecounts a repertoire
  #of size number_of_seqs is simulated
  #name repertoire: additional column that identifies repertoire

  #prepare insertion/deletions df by creating mod1,2,3 subsets for easier access.
  if(verbose==TRUE){
    cat("initializing sim")
  }
  insertion_pairs_mod1 <- insertions_and_deletion_lengths[(nchar(insertions_and_deletion_lengths$n1)+nchar(insertions_and_deletion_lengths$n2))%%3==1,]
  if(verbose==TRUE){
    cat(".")
  }
  insertion_pairs_mod2 <- insertions_and_deletion_lengths[(nchar(insertions_and_deletion_lengths$n1)+nchar(insertions_and_deletion_lengths$n2))%%3==2,]
  if(verbose==TRUE){
    cat(".")
  }
  insertion_pairs_mod3 <- insertions_and_deletion_lengths[(nchar(insertions_and_deletion_lengths$n1)+nchar(insertions_and_deletion_lengths$n2))%%3==0,]
  if(verbose==TRUE){
    cat("\n")
  }
  #set number of sequences to be simulated
  number_of_simulated_seqs<-number_of_seqs

  #set threshold for sim progress output. also allow toggle off of printouts
  if(verbose==TRUE){
    if(number_of_simulated_seqs<=100){
      cat_threshold <- 10
    }else if(number_of_simulated_seqs<=1000){
      cat_threshold <- 100
    }else if(number_of_simulated_seqs<=10000){
      cat_threshold <- 500
    }else if(number_of_simulated_seqs<=100000){
      cat_threshold <- 5000
    }else{
      cat_threshold <- 10000
    }
  }else{
    cat_threshold <- Inf
  }

  simulated_seqs<-list()
  vdj_used<-list()
  vdj_list_ref<-vdj_list

  #common frameshifts and premature cdr3 endings can occur in immuneSIM
  #due to this a lookaround is required for the generated sequences/cdr3s
  #depending on the germline genes the patterns that signify a shift can can look different
  #below we assign checks for these cases.
  #identified common frameshifts (imgt). if these occur avoid them --> extend cdr3 search
  if(species=="hs" & receptor=="ig"){
    common_frameshifts<-c("WTS","LTT","LMS","STP","RSL","**")#,"SST","SIS","GTT","VTT","TTT","RTT","MTT","PTT","GCL","STT")
    check_within_cdr3<-c("LTT","LMS") #GCL  MTT  WTT  TTT  RTT  PTT  SIS  VTT  SST  GTT RSL  STP
    jct_correction<-c("DYW","DSW","FDSW","DVW","DLW","QHW","FDYW","YFDYW","YFDLW","AFDVW","DFDYW","YYFDYW","YYYYYGMDVW","YYYYGMDVW")
  }else if(species=="mm" & receptor=="ig"){
    common_frameshifts<-c("WTT","LLT","WLT","GLT","**")#,"CLL","VLT","TLT","RLT","SLT","PLT")
    check_within_cdr3<-c("LLT","GLT")
    jct_correction<-c("DFDYW","WYFDVW","DYYAMDYW","YYAMDYW","YWYFDVW","WFAYW","YFDYW","FDYW","YAMDYW","FAYW","YFDVW","DVW","AYW","DYW")
  }else if(species=="hs" & receptor=="tr"){
    common_frameshifts<-c("AVF","LHL","HIF","PST","FST","SFL","**") # GFL  AST  TVF PPL  PDF  NCF
    check_within_cdr3<-c("AVF", "LHL","PST","HIF", "SFL","FST")#  PDF  SFL  HIF  PST  LHL   )
    jct_correction<-c("YNEQFF","EQYF","NEQFF","YEQYF","SYNEQFF","SYEQYF","F")
  }else if(species=="mm" & receptor=="tr"){
    common_frameshifts<-c("AVF","LVL","PVL","ALF","IIF","**")
    check_within_cdr3<-"\\*\\*"
    #common_frameshifts<-c("")
    #check_within_cdr3<-c("AVF","TVL","ALF","IIF")
    jct_correction<-c("QNTLYF","GSYEQYF","TEVFF","YAEQFF","DTQYF","EQYF","NYAEQFF","SYEQYF","YEQYF","F")
  }else{
    common_frameshifts<-c("")
    check_within_cdr3<-"\\*\\*"
    jct_correction <- c("")
  }

  count<-number_of_simulated_seqs

  while(count>0){
    #in case we reach a freq update point: update frequencies to improve sampling
    if((count!=number_of_seqs) & ((number_of_seqs-count)%%freq_update_time==0)){

      VDJ_used <- vdj_used
      vdj_list<-.update_vdj_freqs(vdj_list=vdj_list_ref,VDJ_used=vdj_used)
    }

    #if(insertions_deletions_corr==TRUE){
    sample_insertions_and_deletion_lengths<-insertions_and_deletion_lengths[sample(nrow(insertions_and_deletion_lengths),1),]
    sample_insertions<-sample_insertions_and_deletion_lengths[c("n1","n2")]
    sample_deletions<-sample_insertions_and_deletion_lengths[c("del_v","del_d_5","del_d_3","del_j")]
    #}else{
    #  sample_insertions<-insertions_and_deletion_lengths[sample(nrow(insertions_and_deletion_lengths),1),c("n1","n2")]
    #  sample_deletions<-insertions_and_deletion_lengths[sample(nrow(insertions_and_deletion_lengths),1),c("del_v","del_d_5","del_d_3","del_j")]
    #}

    #sample V gene name and lookup V-gene sequence
    sampled_v_name<-as.character(sample(vdj_list$V$gene,1,replace=TRUE,prob=vdj_list$V$frequency))
    sampled_v_seq_raw<-as.character(vdj_list$V[vdj_list$V$gene==sampled_v_name,"sequence"])

    #if seq is longer than sampled deletion apply deletion (check in case very short synthetic Vgenes )
    if(nchar(sampled_v_seq_raw)<=sample_deletions$del_v){
      sample_deletions$del_v<-0
    }
    sampled_v_seq<-base::substr(sampled_v_seq_raw,1,nchar(sampled_v_seq_raw)-sample_deletions$del_v)

    #sample insertion n1
    sampled_n1<-as.character(sample(sample_insertions$n1,1,replace=TRUE))

    #sample D gene name and lookup V-gene sequence
    sampled_d_name<-as.character(sample(vdj_list$D$gene,1,replace=TRUE,prob=vdj_list$D$frequency))
    sampled_d_seq_raw<-as.character(vdj_list$D[vdj_list$D$gene==sampled_d_name,"sequence"])
    #if sampled d-seq is longer than sampled deletion apply deletion
    if(nchar(sampled_d_seq_raw)<=(sample_deletions$del_d_5+sample_deletions$del_d_3)){
      sample_deletions$del_d_5<-0
      sample_deletions$del_d_3<-0
    }
    sampled_d_seq<-base::substr(sampled_d_seq_raw,sample_deletions$del_d_5+1,nchar(sampled_d_seq_raw)-sample_deletions$del_d_3)

    #sample insertion n2
    sampled_n2<-as.character(sample(sample_insertions$n2,1,replace=TRUE))

    #sample J gene name and lookup V-gene sequence
    sampled_j_name<-as.character(sample(vdj_list$J$gene,1,replace=TRUE,prob=vdj_list$J$frequency))
    sampled_j_seq_raw<-as.character(vdj_list$J[vdj_list$J$gene==sampled_j_name,"sequence"])

    anchor_end_name <- sample(c("tgg", "ttt", "ttc"))
    anchor_end <- sapply(1:length(anchor_end_name),function(x) paste("(?=",anchor_end_name[x],")",sep=""))

    #find valid anchor j
    #locate pattern
    locations_j_anchor<-stringr::str_locate_all(sampled_j_seq_raw, anchor_end)
    names(locations_j_anchor)<-anchor_end_name
    #check whether post pattern is in frame.
    #extract anchor and position.
    pos_start_end_j<-NA
    chosen_anchor_j<-""
    for(i in 1:3){
      #check that the option contains a valid start,end point
      j_options<-locations_j_anchor[[i]]
      if(nrow(j_options)>0){
        chosen_anchor_j<-names(locations_j_anchor)[i]
        chosen_row<-sample(nrow(j_options),1)
        pos_start_end_j<-j_options[chosen_row,]
      }
    }

    #update deletion sampling
    sampled_j_seq_raw_pre_anchor<-base::substr(sampled_j_seq_raw,1,pos_start_end_j["start"]-1)

    #check post anchor: make sure that it is not a common frameshift case
    sampled_j_seq_post_anchor<-base::substr(sampled_j_seq_raw,pos_start_end_j[["start"]],nchar(sampled_j_seq_raw))
    if(nchar(sampled_j_seq_post_anchor)>=9){
      check_sequence<-tryCatch(
        {
          as.character(Biostrings::translate(Biostrings::DNAStringSet(substr(sampled_j_seq_post_anchor,1,9))))
        }, error = function(e){
          p<-"**"
          return(p)
        }
      )
      if(check_sequence %in% common_frameshifts){
        chosen_anchor_j<-""
      }
    }

    if(sample_deletions$del_j>nchar(sampled_j_seq_raw_pre_anchor)){
      sample_deletions$del_j<-sample(seq(0,nchar(sampled_j_seq_raw_pre_anchor),1),1)
    }

    sampled_j_seq_pre_anchor<-base::substr(sampled_j_seq_raw_pre_anchor,sample_deletions$del_j+1,nchar(sampled_j_seq_raw_pre_anchor))

    #if sampled j-seq is longer than sampled deletion apply deletion
    sampled_j_seq<-base::substr(sampled_j_seq_raw,(sample_deletions$del_j+1),nchar(sampled_j_seq_raw))

    #optimize insertion to not get out of frame seqs
    #paste together without insertions and find anchors then choose insertions from limited length pool so that it will be in frame
    sequencesim_no_ins_noj<-paste(sampled_v_seq,sampled_d_seq,sampled_j_seq_pre_anchor,sep="")

    #anchor points (same as IGoR) (randomize order to since first one the one to be picked down the line if it fulfills conditions)
    anchor_beg_name <- sample(c("tgt","tgc"))
    anchor_beg <- sapply(1:length(anchor_beg_name),function(x) paste("(?=",anchor_beg_name[x],")",sep=""))

    #take subsequence of relevant length
    jct_estimate <- base::substr(sequencesim_no_ins_noj,nchar(sampled_v_seq)-15,nchar(sequencesim_no_ins_noj))

    #find valid anchor v (locate pattern)
    locations_v_anchor<-stringr::str_locate_all(sampled_v_seq, anchor_beg)
    names(locations_v_anchor)<-anchor_beg_name
    #check whether pre pattern is in frame and far enough down the sequence
    eval_validity_v<-lapply(1:2, function(x) ((locations_v_anchor[[x]][,"start"]-1)%%3==0) & locations_v_anchor[[x]][,"start"]>260)

    #extract anchor and position.
    pos_start_end_v<-NA
    for(i in 1:2){
      if(TRUE %in% eval_validity_v[[i]]){
        pos_start_end_v<-locations_v_anchor[[i]][which(eval_validity_v[[i]]==TRUE),]
      }
    }

    ##check whetherthe position is unique or not (TRBV1 in mm trb can have nonunique positions)
    if(length(pos_start_end_v)>2){
      pos_start_end_v<-pos_start_end_v[sample(1:nrow(pos_start_end_v),1),]
    }

    ##check whether we have a sequence if so choose n1,n2 of lengths that in sum make it in frame
    if(is.na(pos_start_end_v[1]) | chosen_anchor_j==""){
      sequencesim <- "**"
      cdr3_found_tmp<-""
    }else{
      cdr3_estimate_no_ins_noj<-paste(substr(sampled_v_seq,pos_start_end_v["start"],nchar(sampled_v_seq)),sampled_d_seq,sampled_j_seq_pre_anchor,sep="")
      if(nchar(cdr3_estimate_no_ins_noj)%%3==0){
        #choose inserts of sum length 0,3,6,9
        n1_n2<-insertion_pairs_mod3[sample(nrow(insertion_pairs_mod3),1),c("n1","n2")]
        sampled_n1 <- n1_n2[[1,"n1"]]
        sampled_n2 <- n1_n2[[1,"n2"]]
      }else if(nchar(cdr3_estimate_no_ins_noj)%%3==1){
        #choose inserts of sum length 2,5,8,11
        n1_n2<-insertion_pairs_mod2[sample(nrow(insertion_pairs_mod2),1),c("n1","n2")]
        sampled_n1 <- n1_n2[[1,"n1"]]
        sampled_n2 <- n1_n2[[1,"n2"]]
      }else if(nchar(cdr3_estimate_no_ins_noj)%%3==2){
        #choose inserts of sum length 1,4,7,10
        n1_n2<-insertion_pairs_mod1[sample(nrow(insertion_pairs_mod1),1),c("n1","n2")]
        sampled_n1 <- n1_n2[[1,"n1"]]
        sampled_n2 <- n1_n2[[1,"n2"]]
      }
      cdr3_found_tmp<-paste(substr(sampled_v_seq,pos_start_end_v["start"],nchar(sampled_v_seq)),sampled_n1,sampled_d_seq,sampled_n2,sampled_j_seq_pre_anchor,chosen_anchor_j,sep="")
      sequencesim<-paste(sampled_v_seq,sampled_n1,sampled_d_seq,sampled_n2,sampled_j_seq,sep="")

      cdr3_found_tmp_pos <- stringr::str_locate_all(pattern = cdr3_found_tmp, sequencesim)[[1]]

    }


    if(shm.mode!="none"){

      lengths_subregions <- c(v=nchar(sampled_v_seq),n1=nchar(sampled_n1),d=nchar(sampled_d_seq),n2=nchar(sampled_n2),j=nchar(sampled_j_seq))

      sequencesim_pre_shm <- sequencesim
      sequencesim <- .somatic_hypermutation_sim(vdj_seq=sequencesim, mut_param=shm.mode,v_seq=sampled_v_seq,d_seq=sampled_d_seq,j_seq=sampled_j_seq, SHM.nuc.prob=shm.prob)
      sampled_v_seq <- base::substr(sequencesim,1,lengths_subregions["v"])
      running_start <- lengths_subregions["v"]
      sampled_n1 <- base::substr(sequencesim,running_start+1,running_start+lengths_subregions[["n1"]])
      running_start <- running_start+lengths_subregions["n1"]
      sampled_d_seq <- base::substr(sequencesim,running_start+1,running_start+lengths_subregions["d"])
      running_start <- running_start+lengths_subregions["d"]
      sampled_n2 <- base::substr(sequencesim,running_start+1,running_start+lengths_subregions["n2"])
      running_start <- running_start+lengths_subregions["n2"]
      sampled_j_seq <- base::substr(sequencesim,running_start+1,running_start+lengths_subregions["v"])

      if(cdr3_found_tmp!=""){
        cdr3_found_tmp<-base::substr(sequencesim,cdr3_found_tmp_pos[1,"start"],cdr3_found_tmp_pos[1,"end"])
      }

      #identify shm events
      shm_record_out <- .record_shm_events(sequencesim_pre_shm, sequencesim)

    }else{
      shm_record_out <- ""
    }

    #if possible translate to AA
    #if successful keep if no stop codon otherwise discard
    if(cdr3_found_tmp!=""){

      if(nchar(sequencesim)%%3==1){
        sequencesim<-base::substr(sequencesim,1,nchar(sequencesim)-1)
      }else if(nchar(sequencesim)%%3==2){
        sequencesim<-base::substr(sequencesim,1,nchar(sequencesim)-2)
      }

      sequencesim_aa<-tryCatch(
        {
          as.character(Biostrings::translate(Biostrings::DNAStringSet(sequencesim)))
        }, error = function(e){
          p<-"**"
          return(p)
        }
      )

      curr_seq<-sequencesim_aa
      curr_seq_nt<-sequencesim

      curr_jct<-tryCatch(
        {
          as.character(Biostrings::translate(Biostrings::DNAStringSet(cdr3_found_tmp)))
        }, error = function(e){
          p<-"**"
          return(p)
        }
      )

      anchor_end <- sample(c("W", "F"))

      #check if there is a stop codon
      if(!("*" %in% names(table(base::strsplit(sequencesim_aa,split=""))))){
        #identify cdr3 in full sequence (nt and aa)
        locations_beg<-stringr::str_locate_all(curr_seq_nt, cdr3_found_tmp)
        locations_beg_aa<-stringr::str_locate_all(curr_seq, curr_jct)

        #check whether it was found and identify it
        if(length(locations_beg_aa[[1]]) > 0){
          cdr3_to_end<-base::substr(curr_seq_nt,locations_beg[[1]][1,"start"],nchar(curr_seq_nt))
          cdr3_to_end_aa<-base::substr(curr_seq,(locations_beg_aa[[1]][1,"start"]),nchar(curr_seq))

          res <- sapply(1:length(anchor_end),function(x) stringr::str_match(cdr3_to_end_aa, paste(base::substr(cdr3_to_end_aa,1,3),"(.*?)",anchor_end[x],sep="")))[1,]
          cdr3_found_aa<-res[!is.na(res)][1]

          #NEW: correct cdr3 if it ends to early
          #get seq from end of curr cdr3 to end of VDJ
          end_cdr3_to_end_aa<-base::strsplit(cdr3_to_end_aa,cdr3_found_aa)[[1]][2]
          #evaluate whether/which jct correction is needed
          find_jct_correction<-sapply(1:length(jct_correction),function(x) base::substr(end_cdr3_to_end_aa,1,nchar(jct_correction[x]))==jct_correction[x])
          add_jct_corr<-jct_correction[find_jct_correction]
          #if several possibilities choose the longest match (best correction)
          if(length(add_jct_corr)>1){
            corr_ordered_by_length<-add_jct_corr[order(nchar(add_jct_corr), add_jct_corr)]
            add_jct_corr<-corr_ordered_by_length[length(corr_ordered_by_length)]
          }
          #paste together previous cdr3 with correction
          cdr3_found_aa<-paste(cdr3_found_aa,add_jct_corr,sep="")

          #based on found aa get nt cdr3 seq
          cdr3_found_nt<-base::substr(cdr3_to_end,1,nchar(cdr3_found_aa)*3)
          post_cdr3_j <- base::strsplit(curr_seq,cdr3_found_aa)[[1]][2]
          #check for frameshifts
          frameshifted_cdr3 <- FALSE
          for(y in 1:length(check_within_cdr3)){
            if(grepl(check_within_cdr3[y],cdr3_found_aa) | grepl(check_within_cdr3[y],post_cdr3_j)){
              frameshifted_cdr3 <- TRUE
            }
          }

          #if all checks passed addd to simulated
          if(!is.na(cdr3_found_aa) & (nchar(cdr3_found_aa)>=min_cdr3_length) & (nchar(cdr3_found_aa) <= max_cdr3_length) & !frameshifted_cdr3){
            #update user with simulation progress depending (depends on cat thresholds based on number of sequences )
            if(cat_threshold!=Inf){
              if((number_of_seqs-count+1)%%cat_threshold==0){
                if(verbose==TRUE){
                  cat("simulated sequences:", (number_of_seqs-count+1),"\n")
                }
              }
            }
            #save with all essential info
            simulated_seqs[[count]]<-data.frame(sequence=sequencesim,
                                                sequence_aa=sequencesim_aa,
                                                jct_nt=cdr3_found_nt,#ndn_jct, #placeholder for orientation of jct analysis
                                                jct_aa=cdr3_found_aa,#j_jct, #placeholder for orientation of jct analysis #substr(sequencesim_aa,nchar(sequencesim_aa)-(j_jct+ndn_jct+3),nchar(sequencesim_aa)-(j_jct-6)),
                                                v_gene=sampled_v_name,
                                                d_gene=sampled_d_name,
                                                j_gene=sampled_j_name,
                                                n1=sampled_n1,
                                                n2=sampled_n2,
                                                del_v=sample_deletions$del_v,
                                                del_d_5=sample_deletions$del_d_5,
                                                del_d_3=sample_deletions$del_d_3,
                                                del_j=sample_deletions$del_j,
                                                seq_v=substr(sampled_v_seq,pos_start_end_v["start"],nchar(sampled_v_seq)),
                                                seq_d=sampled_d_seq,
                                                seq_j=paste(sampled_j_seq_pre_anchor,chosen_anchor_j,sep=""),
                                                freqs=NA,
                                                counts=NA,
                                                shm_events = shm_record_out,
                                                name_repertoire=name_repertoire
                                                )

            vdj_used[[count]]<-data.frame(v_gene=sampled_v_name,
                                          d_gene=sampled_d_name,
                                          j_gene=sampled_j_name)

            count<-count-1

          }
        }
      }
    }
  }

  simulated_seqs_datatable<-data.table::rbindlist(simulated_seqs)
  simulated_seqs_df <- data.table::setDF(simulated_seqs_datatable)

  simulated_seqs_df
}
