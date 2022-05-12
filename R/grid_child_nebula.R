grid_child_nebula <-
  function(
           graph,
           anno = c(nebula_index = "classification"),
           class,
           layout = "fr",
           title_palette = .MCn.palette_label,
           palette = .MCn.palette,
           print_into = T,
           save_layout_df = NA,
           remove_nodes = F,
           legend_fill = F,
           legend_size = F,
           remove_legend_lab = F,
           theme_args = NA,
           ...
           ){
    ## reformat graph, add with class
    graph <- tidygraph::as_tbl_graph(graph)
    ## nodes ------------------------------------------------------
    nodes <- merge(graph, class, by.x = "name" , by.y = ".id", all.x = TRUE, sort = F)
    nodes <- dplyr::as_tibble(nodes)
    ## edges ------------------------------------------------------
    edges <- dplyr::as_tibble(tidygraph::activate(graph, edges))
    if(nrow(edges) >= 1){
      ## "dotproduct" or other attributes of compare spectra method.
      edges <- dplyr::rename(edges, similarity = 3)
    }else{
      edges <- dplyr::mutate(edges, similarity = NA) 
    }
    ## ------------------------------------------------------------
    ## gather nodes and edges
    graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges)
    ## create network layout
    if(layout == "fr" & nrow(nodes) >= 500)
      layout = "kk"
    layout_n <- create_layout(graph, layout = layout)
    if(is.environment(save_layout_df)){
      assign("layout_n", layout_n, envir = save_layout_df)
    }
    ## plot
    p <- base_vis_c_nebula(layout_n, palette,
                           title = anno[["nebula_index"]],
                           title_fill = title_palette[as.numeric(anno[["hierarchy"]])],
                           remove_nodes = remove_nodes,
                           theme_args = theme_args,
                           ...)
    ## if print into grid panel
    if(!print_into){
      return(p)
    }
    ## remove legend or not
    p <- p + guides(size = ifelse(legend_size, "legend", "none"),
                             fill = ifelse(legend_fill, "legend", "none"))
    ## color bar
    if(class(class$vis_class) == "numeric" & legend_fill){
      p <- p + guides(fill = guide_colorbar(direction = "horizontal", barheight = 0.3))
    }
    ## rm legend labal
    if(remove_legend_lab){
      p <- p + labs(size = "", fill = "")
    }
    print(p, vp = grid::viewport(layout.pos.row = anno[["row"]],
                                 layout.pos.col = anno[["col"]]))
  }
base_vis_c_nebula <-
  function(
           nebula,
           title = NULL,
           palette = .MCn.palette,
           title_fill = "grey",
           nodes_size_range = c(3, 7),
           nodes_stroke = 0.2,
           edges_width_range = c(0.1, 0.7),
           edges_color = "black",
           title_size = 20,
           remove_nodes = F,
           legend_position = "right",
           scale_fill_expression = "scale_fill_manual(values = palette)",
           theme_args = NA,
           ...
           ){
    if(is.vector(attr(palette, "name")))
      palette <- palette[which(names(palette) %in% nebula$vis_class)]
    ## ------------------------------------- 
    if(remove_nodes){
      nodes_size_range = 0
    }
    ## ------------------------------------- 
    theme_args.default <- list(
        text = element_text(family = "Times"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        ## cunstom defined legend position
        legend.position = legend_position,
        plot.title = ggtext::element_textbox(
          ## nebula name textbox
          size = title_size,
          color = "white", fill = title_fill, box.color = "white",
          halign = 0.5, linetype = 1, r = unit(5, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)
        ))
    ## ------------------------------------- 
    if(is.list(theme_args)){
      theme_args.default <- theme_args.default[!names(theme_args.default) %in% names(theme_args)]
      theme_args <- c(theme_args.default, theme_args)
    }else{
      theme_args <- theme_args.default
    }
    ## ------------------------------------- 
    p <- ggraph(nebula) + 
      ## compound MS2 similarity mapping as edge width
      geom_edge_fan(aes(edge_width = similarity), color = edges_color, show.legend = F) + 
      ## nodes size mapping as compound idenfication similarity
      geom_node_point(aes(size = ifelse(!is.na(tanimotoSimilarity),
                                                ## if the tanimotoSimilarity is NA, set to 0.2
                                                tanimotoSimilarity, 0.2),
                                  ## nodes fill. mapping as classification or custom mark
                                  fill = vis_class),
                              shape = 21,
                              stroke = nodes_stroke) + 
      ## for custum commound expression
      eval(parse(text = scale_fill_expression)) +
      ## set range for edge width
      scale_edge_width(range = edges_width_range) + 
      ## for this setting, if range set to 0, remove the nodes
      scale_size(range = nodes_size_range) +
      ## if the nodes is removed, the override.aes setting will retain nodes shape and color in legend
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ## the title is the annotation of classification
      ggtitle(stringr::str_wrap(title, width = 30)) +
      labs(size = "Tanimoto\nsimilarity", fill = "Compound mark") +
      theme_grey() +
      do.call(theme, theme_args)
    return(p)
  }
