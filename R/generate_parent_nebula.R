#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param write_output PARAM_DESCRIPTION, Default: T
#' @param output_format PARAM_DESCRIPTION, Default: 'graphml'
#' @param output PARAM_DESCRIPTION, Default: paste0(.MCn.output, "/", .MCn.results)
#' @param edges_file PARAM_DESCRIPTION, Default: paste0(output, "/parent_nebula/parent_nebula_edges.tsv")
#' @param edges_method PARAM_DESCRIPTION, Default: 'method_formula_based_spec_compare'
#' @param nodes_attributes PARAM_DESCRIPTION, Default: .MCn.formula_set
#' @param nodes_other_attributes PARAM_DESCRIPTION, Default: .MCn.structure_set
#' @param edge_filter PARAM_DESCRIPTION, Default: 0.5
#' @param cpu_cores PARAM_DESCRIPTION, Default: 8
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
#'  \code{\link[dplyr]{reexports}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{mutate_all}}, \code{\link[dplyr]{rename}}
#'  \code{\link[igraph]{as_data_frame}}
#' @rdname generate_parent_nebula
#' @export 
#' @importFrom dplyr as_tibble mutate mutate_at rename
#' @importFrom igraph graph_from_data_frame
generate_parent_nebula <-
  function(
           write_output = T,
           output_format = "graphml",
           output = paste0(.MCn.output, "/", .MCn.results),
           # exists edges file
           edges_file = paste0(output, "/parent_nebula/parent_nebula_edges.tsv"), 
           # or NULL
           edges_method = "method_formula_based_spec_compare", 
           nodes_attributes = .MCn.formula_set,
           nodes_other_attributes = .MCn.structure_set,
           rm_parent_isolate_nodes = F,
           edge_filter = 0.5,
           cpu_cores = 8, 
           ...
           ){
    cat("[INFO] MCnebula run: generate_parent_nebula\n")
    ## main body
    ## ---------------------------------------------------------------------- 
    ## generate edges data
    if(is.null(edges_method)){
      ## no edges_method
      cat("## generate_parent_nebula: no edges_uethod used\n")
      edges <- dplyr::as_tibble(cbind(".id_1" = nodes_other_attributes$".id",
                     ".id_2" = nodes_other_attributes$".id")) %>%
        dplyr::mutate(dotproduct = 1, mass_diff = 0)
    }else if(edges_method == "method_formula_based_spec_compare"){
      ## with edges_method
      if(!is.null(edges_file) & file.exists(edges_file)){
        cat("## generate_parent_nebula: file.exists(edges_file) == T. Escape from time-consuming computation\n")
        edges <- read_tsv(edges_file) %>%
          dplyr::mutate_at(c(".id_1", ".id_2"), as.character) %>%
          dplyr::mutate_at(c(colnames(.)[3:4]), as.numeric)
      }else{
        cat("## generate_parent_nebula: method_formula_based_spec_compare\n")
        edges = method_formula_based_spec_compare(edge_filter = edge_filter, cpu_cores = cpu_cores, ...)
      }
    }
    ## ---------------------------------------------------------------------- 
    ## generate nodes data
    nodes <- nodes_attributes
    ## additional nodes attributes
    if(is.data.frame(nodes_other_attributes)){
      nodes <- merge(nodes, nodes_other_attributes, by = ".id", all.x = T, sort = T) %>%
        ## rename the column name, otherwise the column will be choosed as key column in igraph
        dplyr::rename(compound_name = name)
    }
    if(rm_parent_isolate_nodes){
      non_iso_nodes <- unique(c(edges$.id_1, edges$.id_2))
      nodes.parent <- dplyr::filter(nodes, .id %in% all_of(non_iso_nodes))
    }else{
      nodes.parent <- nodes
    }
    ## ---------------------------------------------------------------------- 
    ## graph
    parent_nebula <- igraph::graph_from_data_frame(edges, directed = T, vertices = nodes.parent)
    if(write_output){
      dir = paste0(output, "/", "parent_nebula")
      if(!file.exists(dir)){
        dir.create(dir)
      }
      igraph::write_graph(parent_nebula,
                  file = paste0(dir, "/", "parent_nebula.", output_format),
                  format = output_format)
      write_tsv(edges, paste0(dir, "/", "parent_nebula_edges.tsv"))
      write_tsv(nodes, paste0(dir, "/", "parent_nebula_nodes.tsv"))
    }
    ## ---------------------------------------------------------------------- 
    ## set as global var for next stage
    .MCn.parent_graph <<- parent_nebula
    .MCn.parent_nodes <<- as_tibble(nodes)
    .MCn.parent_edges <<- as_tibble(edges)
    cat("[INFO] MCnebula Job Done: generate_parent_nebula\n")
  }
