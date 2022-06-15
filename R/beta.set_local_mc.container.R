## set some empty var to local to storage var
set_local_mc.container <- 
  function(
           FUN
           ){
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
                                                ".MCn.palette_label")
    )
    if("all" %in% FUN){
      set <- overall_set
    }else{
      set <- overall_set[names(overall_set) %in% FUN]
    }
    set <- unlist(set, use.names = F) 
    ## set var in parent frame
    lapply(set, function(var){
             assign(var, 0, envir = parent.frame(3))
             print(var)
             print(parent.frame(3))
    })
    print(ls())
    print("---")
    print(parent.env(environment()))
  }
