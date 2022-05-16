annotate_node <- 
  function(
           node_id,
           compound_class_list = .MCn.nebula_class,
           output = paste0(.MCn.output, "/", .MCn.results),
           plot_nodes_id = T,
           plot_structure = T,
           plot_ppcp = T,
           ratio_df = NA,
           nodes_mark = NA,
           annotate_ppcp.class.id = T,
           ...
           ){
    nodes <- node_id
    ## get top compound class (nodes_color data)
    ## as well as, collate metadata
    metadata <- lapply(compound_class_list, head, n = 1) %>%
      data.table::rbindlist(idcol = T) %>% # as data.frame
      dplyr::filter(.id %in% nodes) %>% # filter via nodes
      dplyr::select(.id, name) %>%
      dplyr::rename(vis_class = name)
    ## ---------------------------------------------------------------------- 
    ## mark nodes in color
    if(is.data.frame(nodes_mark)){
      ## the secound col as mark col
      colnames(nodes_mark) <- c(".id", "mark")
      ## merge with metadata
      metadata <- merge(metadata, nodes_mark, by = ".id", all.x = T) %>% 
        dplyr::mutate(vis_class = ifelse(is.na(mark), "Others", mark))
    }
    ## ------------------------------------- 
    if(is.data.frame(ratio_df)){
      plot_ratio <- T
    }else{
      plot_ratio <- F
    }
    ## ------------------------------------- 
    tmp_dir <- paste0(output, "/", "tmp")
    ## ------------------------------------- 
    ## plot 2D structure, require ChemmineOB and ChemmineR
    with_structure <- 0
    if(requireNamespace("ChemmineOB", quietly = T)){
      ## structure
      tmp_stru <- paste0(tmp_dir, "/", "structure")
      if(!file.exists(tmp_stru)){
        dir.create(tmp_stru)
      }
      if(plot_structure){
        with_structure <- 1
        batch_mode_structure(metadata = metadata, tmp_stru = tmp_stru)
      }
    }
    ## ------------------------------------- 
    tmp_ppcp <- paste0(tmp_dir, "/", "ppcp")
    ## draw nodes with class id number
    if(plot_ppcp | plot_structure | plot_ratio ){
      do.call(batch_mode_nodes, list(metadata = metadata,
                                     tmp_ppcp = tmp_ppcp,
                                     with_structure = with_structure,
                                     plot_ppcp = plot_ppcp,
                                     plot_ratio = plot_ratio,
                                     ratio_df = ratio_df,
                                     annotate_ppcp.class.id = annotate_ppcp.class.id,
                                     get_ppcp_legend = T,
                                     ...))
      filepath <- paste0(tmp_ppcp, "/", node_id, ".svg")
      ## mv file
      file.copy(filepath, output)
      ## ------------------------------------- 
    }
  }
