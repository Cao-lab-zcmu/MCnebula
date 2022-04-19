#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param graph_list PARAM_DESCRIPTION, Default: .MCn.child_graph_list
#' @param compound_class_list PARAM_DESCRIPTION, Default: .MCn.nebula_class
#' @param output PARAM_DESCRIPTION, Default: paste0(.MCn.output, "/", .MCn.results)
#' @param layout PARAM_DESCRIPTION, Default: 'fr'
#' @param width PARAM_DESCRIPTION, Default: 23
#' @param height PARAM_DESCRIPTION, Default: 30
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{rename}}, \code{\link[dplyr]{reexports}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{mutate}}
#'  \code{\link[svglite]{svglite}}
#'  \code{\link[grid]{grid.newpage}}, \code{\link[grid]{Working with Viewports}}
#'  \code{\link[pbapply]{pbapply}}
#' @rdname visualize_child_nebulae
#' @export 
#' @importFrom data.table rbindlist
#' @importFrom dplyr select rename as_tibble arrange mutate
#' @importFrom svglite svglite
#' @importFrom grid grid.newpage pushViewport
#' @importFrom pbapply pbmapply
visualize_child_nebulae <-
  function(
           graph_list = .MCn.child_graph_list,
           compound_class_list = .MCn.nebula_class,
           output = paste0(.MCn.output, "/", .MCn.results),
           layout = "fr",
           width = 23,
           height = 30,
           nodes_mark = NA,
           ...
           ){
    cat("[INFO] MCnebula run: visualize_child_nebulae\n")
    ## get top compound class (nodes_color data)
    metadata <- lapply(compound_class_list, head, n = 1) %>%
      data.table::rbindlist(idcol = T) %>%
      dplyr::select(.id, name) %>%
      dplyr::rename(vis_class = name)
    ## ------------------------------------- 
    if(is.data.frame(nodes_mark)){
      ## the secound col as mark col
      colnames(nodes_mark) <- c(".id", "mark")
      ## merge with metadata
      metadata <- merge(metadata, nodes_mark, by = ".id", all.x = T) %>% 
        dplyr::mutate(vis_class = ifelse(is.na(mark), "Others", mark))
    }
    ## ---------------------------------------------------------------------- 
    ## draw network via ggplot, and print into grid palette
    ## number of child_nebulae
    n = length(graph_list)
    ## specification of grid (cols * rows)
    cols = n^(1/2)
    if(round(cols) != cols){
      cols = round(cols)
      rows = cols + 1
    }else{
      rows = cols
    }
    ## ------------------------------------- 
    ## grid position of all child_nebulae
    graph_anno <- names(graph_list) %>% # names
      dplyr::as_tibble() %>%
      dplyr::rename(nebula_index = value) %>%
      merge(.MCn.class_tree_data[,c("name", "hierarchy")], by.x = "nebula_index", by.y = "name", all.x = T) %>%
      dplyr::arrange(desc(hierarchy)) %>%
      ## calculate position
      dplyr::mutate(seq = 1:n, 
                    col = ifelse(seq %% cols != 0, seq %% cols, cols),
                    row = (seq - col)/cols + 1)
    ## ------------------------------------- 
    ## re-set rows
    rows <- max(graph_anno$row)
    ## as list
    nebula_index <- graph_anno$nebula_index
    graph_anno <- by_group_as_list(graph_anno, "nebula_index")
    ## re-order the graph list according to annotation
    graph_list <- lapply(nebula_index, function(x){
                           graph_list[[x]]
                         })
    ## ------------------------------------- 
    ## prepare grid panel
    svglite::svglite(paste0(output, "/", "child_nebulae.svg"), width = width, height = height)
    grid::grid.newpage()
    grid::pushViewport(viewport(layout = grid.layout(rows, cols)))
    ## draw child_nebulae in grid
    pbapply::pbmapply(grid_child_nebula, ## function
                      graph_list, ## graph list
                      graph_anno, ## graph annotation
                      MoreArgs = list( ## args
                                      class = metadata,
                                      layout = layout,
                                      ...
                                      ))
    dev.off()
    cat("[INFO] MCnebula Job Done: visualize_child_nebulae\n")
  }
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
    if(!is.na(save_layout_df)){
      assign("layout_n", layout_n, envir = save_layout_df)
    }
    ## plot
    p <- base_vis_c_nebula(layout_n, palette,
                           title = anno[["nebula_index"]],
                           title_fill = title_palette[as.numeric(anno[["hierarchy"]])],
                           remove_nodes = remove_nodes,
                           ...)
    if(!print_into){
      return(p)
    }
    p <- p + guides(size = ifelse(legend_size, "legend", "none"),
                             fill = ifelse(legend_fill, "legend", "none"))
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
           edges_width_range = c(0.1, 0.7),
           title_size = 20,
           remove_nodes = F,
           legend_position = "right",
           scale_fill_expression = "scale_fill_manual(values = palette)",
           ...
           ){
    if(is.vector(attr(palette, "name")))
      palette <- palette[which(names(palette) %in% nebula$vis_class)]
    ## ------------------------------------- 
    if(remove_nodes){
      nodes_size_range = 0
    }
    ## ------------------------------------- 
    p <- ggraph(nebula) + 
      ## compound MS2 similarity mapping as edge width
      geom_edge_fan(aes(edge_width = similarity), color = "black", show.legend = F) + 
      ## nodes size mapping as compound idenfication similarity
      geom_node_point(aes(size = ifelse(is.na(tanimotoSimilarity) == F,
                                                ## if the tanimotoSimilarity is NA, set to 0.2
                                                tanimotoSimilarity, 0.2),
                                  ## nodes fill. mapping as classification or custom mark
                                  fill = vis_class),
                              shape = 21) + 
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
      labs(size = "Tanimoto\nsimilarity", fill = "Compound class") +
      theme_grey() +
      theme(
        text = element_text(family = "Times"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        ## cunstom defined legend position
        legend.position = legend_position,
        panel.grid = element_blank(),
        plot.title = ggtext::element_textbox(
          ## nebula name textbox
          size = title_size,
          color = "white", fill = title_fill, box.color = "white",
          halign = 0.5, linetype = 1, r = unit(5, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)
        )
      )
    return(p)
  }
