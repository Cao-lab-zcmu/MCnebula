#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param apset PARAM_DESCRIPTION
#' @param cluster_cutoff PARAM_DESCRIPTION, Default: 0.7
#' @param generate_neighbors PARAM_DESCRIPTION, Default: 20
#' @param share_numbers PARAM_DESCRIPTION, Default: 10
#' @param share_limite PARAM_DESCRIPTION, Default: 'average'
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
#'  \code{\link[ChemmineR]{jarvisPatrick}}
#'  \code{\link[dplyr]{rename}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{distinct}}
#'  \code{\link[tidyr]{separate}}
#' @rdname method_rerank_jarvis_patrick_cluster
#' @export 
#' @importFrom ChemmineR jarvisPatrick
#' @importFrom dplyr rename mutate arrange distinct
#' @importFrom tidyr separate
method_rerank_jarvis_patrick_cluster <-
  function(
           apset,
           cluster_cutoff = 0.7,
           generate_neighbors = 20,
           share_numbers = 10,
           share_limite = "average",
           ...
           ){
    ## cluster via function of ChemmineR
    cat("## method_rerank_jarvis_patrick_cluster: ChemmineR::jarvisPatrick\n")
    meta_rank <- ChemmineR::jarvisPatrick(nearestNeighbors(apset,
                                                numNbrs = generate_neighbors,
                                                cutoff = cluster_cutoff),
                               k = share_numbers,
                               linkage = share_limite,
                               mode = "a1b") %>%
      vector_as_df() %>%
      merge(., vector_as_df(table(.$expr)), by.x = "expr", by.y = "ids", all.x = T) %>%
      dplyr::rename(cluster_id = expr, size = expr.y) %>%
      tidyr::separate(col = "ids", into = c(".id", "structure_rank"), sep = "_", remove = T) %>%
      dplyr::mutate(structure_rank = as.numeric(structure_rank)) %>%
      dplyr::arrange(desc(size), structure_rank) %>%
      ## for each .id, only the top 1 (according to 'size', 'structure_rank', sequentialy) retain
      dplyr::distinct(.id, .keep_all = T)
    return(meta_rank)
  }
vector_as_df <-
  function(
           vector
           ){
    df <- data.table(cbind(ids = names(vector), expr = unname(vector)))
    return(df)
  }
