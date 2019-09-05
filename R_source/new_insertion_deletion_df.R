#new instertions&deletions
#dropout insertions or deletions or any combination thereof
#(no resampling options as sampling in sim generates enough diversity regarding ins,del, dropout is needed for special repertoires however.

.new_insertion_deletion_df<-function(new_insertion_deletion,options){
  #go through the options chosen by the user and
  #make a new_insertion_deletion dataframe that conforms with
  #chosen parameters (if insertion dropout ="", if deletion dropout = 0)
  for(i in 1:length(options)){
    option<-options[i]

    if(option=="no_insertions"){
      new_insertion_deletion$n1<-""
      new_insertion_deletion$n2<-""
    } else if(option=="no_insertions_n1"){
      new_insertion_deletion$n1<-""
    } else if(option=="no_insertions_n2"){
      new_insertion_deletion$n2<-""
    } else if(option=="no_deletions_v"){
      new_insertion_deletion$del_v<-0
    } else if(option=="no_deletions_d"){
      new_insertion_deletion$del_d_5<-0
      new_insertion_deletion$del_d_3<-0
    } else if(option=="no_deletions_d_5"){
      new_insertion_deletion$del_d_5<-0
    } else if(option=="no_deletions_d_3"){
      new_insertion_deletion$del_d_3<-0
    } else if(option=="no_deletions_j"){
      new_insertion_deletion$del_j<-0
    } else if(option=="no_deletions_vd"){
      new_insertion_deletion$del_v<-0
      new_insertion_deletion$del_d_5<-0
      new_insertion_deletion$del_d_3<-0
     # new_insertion_deletion$del_j<-0
    } else if(option=="no_deletions"){
      new_insertion_deletion$del_v<-0
      new_insertion_deletion$del_d_5<-0
      new_insertion_deletion$del_d_3<-0
      new_insertion_deletion$del_j<-0
    } else if(option=="no_insertions_and_deletions"){
      new_insertion_deletion$n1<-""
      new_insertion_deletion$n2<-""
      new_insertion_deletion$del_v<-0
      new_insertion_deletion$del_d_5<-0
      new_insertion_deletion$del_d_3<-0
      new_insertion_deletion$del_j<-0
    }
  }

  return(new_insertion_deletion)

}
