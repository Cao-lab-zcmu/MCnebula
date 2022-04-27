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
        dplyr::mutate(vis_class = ifelse(!is.na(mark), mark,
                                         ifelse(is.numeric(mark), NA, "Others")))
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
