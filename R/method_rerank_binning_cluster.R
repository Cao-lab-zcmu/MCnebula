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
           structure_set,
           reference_compound,
           cluster_cutoff = seq(0.9, 0.1, by = -0.01),
           csi_score_weight,
           class_similarity_weight
           ){
    cat("## convert data: SMILES_set -> SDF_set -> AP_set\n")
    sdfset <- smiles_to_sdfset(structure_set)
    apset <- sdf2ap(sdfset)
    ## ---------------------------------------------------------------------- 
    ## cluster via function of ChemmineR
    cat("## method_rerank_binning_cluster: ChemmineR::cmp.cluster\n")
    meta_rank <- ChemmineR::cmp.cluster(db = apset, cutoff = cluster_cutoff)
    ## ---------------------------------------------------------------------- 
    ## preparation for calculating
    meta_rank <- binning_reshape_meta1(meta_rank, reference_compound)
    ## ---------------------------------------------------------------------- 
    cluster_score <- binning_get_cluster_score(meta_rank)
    ## ---------------------------------------------------------------------- 
    ## merge to get cluster_score
    meta_rank <- binning_reshape_meta2(meta_rank, cluster_score)
    ## ---------------------------------------------------------------------- 
    structure_set <- binning_reshape_structure_set(structure_set, meta_rank,
                                                   csi_score_weight = csi_score_weight,
                                                   class_similarity_weight = class_similarity_weight)
    ## ---------------------------------------------------------------------- 
    return(structure_set)
  }
binning_reshape_meta1 <- 
  function(
           meta_rank,
           reference_compound
           ){
    meta_rank <- dplyr::select(meta_rank, ids, starts_with("CLID"))
    ## reshape the data as long data frame
    meta_rank <- reshape2::melt(meta_rank, id.var = "ids", variable.name = "name", value.name = "number")
    ## get the cutoff of cluster results
    meta_rank <- dplyr::mutate(meta_rank, cutoff = stringr::str_extract(name, "(?<=_).{1,}$"),
                               ## get origin .id
                               .id = stringr::str_extract(ids, ".*(?=_)"),
                               ## get origin structure_rank
                               structure_rank = stringr::str_extract(ids, "(?<=_)[0-9]{1,}$"),
                               ## whether reference compound
                               reference = ifelse(.id %in% reference_compound$.id & structure_rank == 1, T, F))
    ## merge to get tanimotoSimilarity of reference compound
    meta_rank <- merge(meta_rank,
                       reference_compound[, c(".id", "tanimotoSimilarity")],
                       by = ".id", all.x = T)
    ## if the compound is reference compound, assign tanimotoSimilarity as reference_score
    meta_rank <- dplyr::mutate(meta_rank,
                               reference_score = ifelse(reference, tanimotoSimilarity, 0),
                               ## set group, according to cutoff and cluster ID
                               group = paste0(name, "_", number))
    return(meta_rank)
  }
binning_get_cluster_score <- 
  function(
           meta_rank
           ){
    ## to calculate score of each cluster, set group
    cluster_score <- dplyr::group_by(meta_rank, group)
    ## calculate
    cluster_score <- dplyr::summarise_at(cluster_score, "reference_score", sum)
    ## rename
    cluster_score <- dplyr::rename(cluster_score, cluster_score = reference_score)
    return(cluster_score)
  }
binning_reshape_meta2 <- 
  function(
           meta_rank,
           cluster_score
           ){
    meta_rank <- merge(meta_rank, cluster_score, by = "group", all.x = T)
    ## according to cutoff, re-size the score
    meta_rank <- dplyr::mutate(meta_rank, cluster_score = cluster_score * as.numeric(cutoff)^3)
    ## group via id and structure_rank
    meta_rank <- dplyr::group_by(meta_rank, ids)
    ## for all scores of one id, get sum score
    meta_rank <- dplyr::summarise_at(meta_rank, "cluster_score", sum)
    ## ---------------------------------------------------------------------- 
    ## reformat and output
    ## ------------------------------------- 
    meta_rank <- dplyr::mutate(meta_rank,
                               ## get origin .id
                               .id = stringr::str_extract(ids, ".*(?=_)"),
                               ## get origin structure_rank
                               structure_rank = stringr::str_extract(ids, "(?<=_)[0-9]{1,}$"),
                               ## as.numeric
                               structure_rank = as.numeric(structure_rank))
    ## by id as list
    meta_rank <- by_group_as_list(meta_rank, ".id")
    ## normalize the score, via dividing by top score
    meta_rank <- lapply(meta_rank, function(df){
                          dplyr::mutate(df, norm_score = cluster_score / max(cluster_score))
                               })
    meta_rank <- data.table::rbindlist(meta_rank)
    return(meta_rank)
  }
binning_reshape_structure_set <- 
  function(
           structure_set,
           meta_rank,
           csi_score_weight,
           class_similarity_weight
           ){
    structure_set <- merge(structure_set,
                           meta_rank[, c(".id", "structure_rank", "norm_score")],
                           by = c(".id", "structure_rank"), all.x = T)
    ## ------------------------------------- 
    ## get comprehensive score
    structure_set <- by_group_as_list(structure_set, ".id") %>% 
      lapply(function(df){
               dplyr::mutate(df, norm_csi_score = (score / max(score))^-1)
                           }) %>% 
      data.table::rbindlist()
    ## ------------------------------------- 
    ## gather score
    structure_set <- structure_set %>% 
      dplyr::mutate(re_rank_score = csi_score_weight * norm_csi_score +
                    class_similarity_weight * norm_score) %>% 
      dplyr::arrange(.id, desc(re_rank_score)) %>% 
      dplyr::distinct(.id, .keep_all = T) %>% 
      dplyr::as_tibble()
    return(structure_set)
  }
