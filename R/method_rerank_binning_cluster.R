#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param apset PARAM_DESCRIPTION
#' @param cluster_cutoff PARAM_DESCRIPTION, Default: seq(0.9, 0.1, by = -0.1)
#' @param least_size PARAM_DESCRIPTION, Default: 3
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[ChemmineR]{cmp.cluster}}
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{filter}}, \code{\link[dplyr]{distinct}}
#'  \code{\link[reshape2]{melt}}
#'  \code{\link[tidyr]{separate}}
#' @rdname method_rerank_binning_cluster
#' @export 
#' @importFrom ChemmineR cmp.cluster
#' @importFrom dplyr select mutate arrange filter distinct
#' @importFrom reshape2 melt
#' @importFrom tidyr separate
method_rerank_binning_cluster <-
  function(
           apset,
           cluster_cutoff = seq(0.9, 0.1, by = -0.1),
           least_size = 3
           ){
    ## cluster via function of ChemmineR
    cat("## method_rerank_binning_cluster: ChemmineR::cmp.cluster\n")
    meta_rank <- ChemmineR::cmp.cluster(db = apset, cutoff = cluster_cutoff) %>%
      dplyr::select(ids, starts_with("CLSZ_")) %>%
      ## convert into long table
      reshape2::melt(id.var = "ids", variable.name = "cutoff", value.name = "size") %>%
      ## get 'cutoff' and as.numeric
      dplyr::mutate(cutoff = as.numeric(gsub("CLSZ_", "", cutoff))) %>%
      ## get '.id' and 'structure_rank'
      tidyr::separate(col = "ids", into = c(".id", "structure_rank"), sep = "_", remove = T) %>%
      dplyr::mutate(structure_rank = as.numeric(structure_rank)) %>%
      dplyr::arrange(desc(cutoff), desc(size), structure_rank) %>%
      ## at least, the size of cluster reach 'least_size', contribute to re-rank
      dplyr::filter(!(cutoff >= min(cluster_cutoff) & size <= least_size)) %>%
      ## for each .id, only the top 1 (according to 'cutoff', 'size', 'structure_rank', sequentialy) retain
      dplyr::distinct(.id, .keep_all = T)
    return(meta_rank)
  }
