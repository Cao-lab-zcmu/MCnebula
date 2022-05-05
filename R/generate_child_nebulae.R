#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nodes PARAM_DESCRIPTION, Default: .MCn.parent_nodes
#' @param edges PARAM_DESCRIPTION, Default: .MCn.parent_edges
#' @param max_edges PARAM_DESCRIPTION, Default: 5
#' @param nebula_index PARAM_DESCRIPTION, Default: .MCn.nebula_index
#' @param output PARAM_DESCRIPTION, Default: paste0(.MCn.output, "/", .MCn.results)
#' @param output_format PARAM_DESCRIPTION, Default: 'graphml'
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
#'  \code{\link[pbapply]{pbapply}}
#' @rdname generate_child_nebulae
#' @export 
#' @importFrom pbapply pblapply
generate_child_nebulae <-
  function(
           nodes = .MCn.parent_nodes,
           edges = .MCn.parent_edges,
           max_edges = 5,
           nebula_index = .MCn.nebula_index,
           output = paste0(.MCn.output, "/", .MCn.results),
           output_format = "graphml",
           ...
           ){
    ##
    cat("[INFO] MCnebula run: generate_child_nebulae\n")
    assign("envir_nebula", environment(), envir = parent.env(environment()))
    ## get names of all classes
    names <- unique(nebula_index$name)
    ## for using lapply, first, trans the data.frame into list
    nebula_index <- by_group_as_list(nebula_index, "relativeIndex")
    ## push names
    names(nebula_index) <- names
    ## facet parent nebula
    ## create dir for placing graph file
    dir = paste0(output, "/", "child_nebula")
    if(file.exists(dir) == F){
      dir.create(dir)
    }
    .MCn.child_graph_list <<- pbapply::pblapply(nebula_index, separate_nebula,
                      output = dir,
                      max_edges = max_edges,
                      output_format = output_format,
                      ...)
    cat("[INFO] MCnebula Job Done: generate_child_nebulae\n")
  }
separate_nebula <-
  function(
           df,
           nodes = get("nodes", envir = get("envir_nebula")),
           edges = get("edges", envir = get("envir_nebula")),
           write_output = T,
           output = paste0(.MCn.output, "/", .MCn.results, "/", "child_nebula"),
           output_format = "graphml",
           max_edges = 5,
           write_extra = F
           ){
    id <- unique(df$".id")
    ## get the child nebula name
    ## note that some character in name caused fail to write as file into dir
    name <- gsub("/", "#", df[1,]$"name")
    nodes <- nodes[nodes$".id" %in% id, ]
    edges <- edges[edges$".id_1" %in% id & edges$".id_2" %in% id, ] 
    ## an edges number cut-off
    edges <- better_vis_nebula(edges, max_edges = max_edges)
    child_nebula <- igraph::graph_from_data_frame(edges, directed = T, vertices = nodes)
    if(write_output == T){
      igraph::write_graph(child_nebula,
                  file = paste0(output, "/", name, ".", output_format),
                  format = output_format)
      if(write_extra == T){
        write_tsv(edges, paste0(output, "/", name, "_edges.tsv"))
        write_tsv(nodes, paste0(output, "/", name, "_nodes.tsv"))
      }
    }
    return(child_nebula)
  }
better_vis_nebula <-
  function(
           edges,
           max_edges = 5
           ){
    ## order
    edges <- dplyr::arrange(edges, desc(edges[,3]))
    ta <- table(c(edges[[1]], edges[[2]]))
    ## at least loop number
    n <- length(which(ta > max_edges))
    if(n == 0){
      return(edges)
    }
    ## ----------------------------------------------------------------------
    ## copy data for override
    df <- dplyr::mutate(edges[, 1:2], SEQ = 1:nrow(edges))
    assign("envir_meta", environment(), envir = parent.env(environment()))
    ## use sapply instead of while loop
    continue = 1
    sapply(1:n, edges_cut_off, max_edges = max_edges)
    ## ---------------------------------------------------------------------- 
    edges <- edges[df$SEQ, ]
    return(edges)
  }
edges_cut_off <-
  function(
           i,
           max_edges = 5
           ){
    continue = get("continue", envir = get("envir_meta"))
    if(continue == 1){
      edges = get("df", envir = get("envir_meta"))
      ## stat edges number of an id
      ta <- table(c(edges[[1]], edges[[2]]))
      ## ---------------------------------------------------------------------- 
      if(max(ta) > max_edges){ ## greater than threshold, hence perform exclude
        ## select an id to exclude its excess edges
        key_id <- names(ta[which(ta == max(ta))])[1]
        ## get SEQ of the edges which need to be excluded
        incude_id_edges <- edges[which(edges[[1]] == key_id | edges[[2]] == key_id),]
        exclude_edges_seq <- incude_id_edges[-(1:max_edges), ]$SEQ
        ## exclude edges
        edges <- edges[which(!edges$SEQ %in% exclude_edges_seq), ]
        assign("df", edges, envir = get("envir_meta"))
      }else{
      ## ---------------------------------------------------------------------- 
        ## signature of stop exclude
        assign("continue", 0, envir = get("envir_meta"))
      }
    }else{
      return()
    }
  }
