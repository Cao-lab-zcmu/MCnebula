## set some empty var to local to storage var
call_fun_mc.space <- 
  function(
           FUN,
           args,
           clear_start = T,
           clear_end = T
           ){
    local <- environment()
    if(clear_start){
      rm_mc.set(envir = parent.env(local))
    }
    ## ---------------------------------------------------------------------- 
    overall_set <- get_mc.global_meta()
    set <- overall_set[names(overall_set) == FUN]
    set <- unlist(set, use.names = F) 
    ## ---------------------------------------------------------------------- 
    lapply(set, function(var){
             assign(var, 0, envir = parent.env(local))
    })
    ## ----------------------------------------------------------------------
    res <- do.call(match.fun(FUN), args)
    ## ------------------------------------- 
    res <- list(envir = parent.env(local), results = res)
    ## ------------------------------------- 
    if(clear_end){
      rm_mc.set(envir = parent.env(local))
    }
    return(res)
  }
get_mc.global_meta <- 
  function(){
    overall_set <- list(build_classes_tree_list = c(".MCn.class_tree_list"),
                        collate_ppcp = c(".MCn.ppcp_dataset",
                                         ".MCn.class_tree_data",
                                         ".MCn.nebula_class",
                                         ".MCn.nebula_index"),
                        collate_structure = c(".MCn.formula_set",
                                              ".MCn.structure_set"),
                        generate_child_nebulae = c(".MCn.child_graph_list"),
                        generate_parent_nebula = c(".MCn.parent_graph",
                                                   ".MCn.parent_nodes",
                                                   ".MCn.parent_edges"),
                        initialize_mcnebula = c(".MCn.sirius",
                                                ".MCn.output",
                                                ".MCn.results",
                                                ".MCn.palette",
                                                ".MCn.palette_stat",
                                                ".MCn.palette_ppcp",
                                                ".MCn.palette_label"))
    return(overall_set)
  }
