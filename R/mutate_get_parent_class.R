mutate_get_parent_class <- 
  function(
           classes,
           class_cutoff = 4,
           meta = .MCn.class_tree_data,
           this_class = F
           ){
    ## ------------------------------------- 
    cat("## get_parent_class\n")
    db <- dplyr::filter(meta, hierarchy >= class_cutoff)
    ## -------------------------------------
    ## for name to get id
    db_id <- lapply(db$id, c)
    names(db_id) <- db$name
    ## for id to get parentId
    db_parent <- lapply(db$parentId, c)
    names(db_parent) <- db$id
    ## for id to get name
    db_name <- lapply(db$name, c)
    names(db_name) <- db$id
    ## -------------------------------------
    set_list <- pbapply::pblapply(classes, get_parent_class,
                                  db_id = db_id,
                                  db_parent = db_parent,
                                  db_name = db_name,
                                  this_class = this_class)
    names(set_list) <- classes
    return(set_list)
  }
get_parent_class <- 
  function(
           class,
           db_id,
           db_parent,
           db_name,
           this_class = F
           ){
    set <- c()
    parent <- 0
    id <- db_id[[class]]
    test <- try(db_parent[[id]], silent = T)
    if(class(test) == "try-error")
      return()
    while(is.null(parent) == F){
      if(parent != 0){
        set <- c(set, db_name[[parent]])
        id <- parent
      }
      parent <- db_parent[[id]]
    }
    if(length(set) == 0){
      if(this_class == T)
        return(class)
    }
    return(set)
  }
