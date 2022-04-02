#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nebula_name PARAM_DESCRIPTION
#' @param top_n PARAM_DESCRIPTION, Default: 10
#' @param match_pattern PARAM_DESCRIPTION, Default: c("precursorFormula")
#' @param collate_factor PARAM_DESCRIPTION, Default: 0.85
#' @param revise_MCn_formula_set PARAM_DESCRIPTION, Default: T
#' @param revise_MCn_structure_set PARAM_DESCRIPTION, Default: T
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
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{reexports}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{distinct}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[tidyr]{separate}}
#' @rdname nebula_re_rank
#' @export 
#' @importFrom dplyr filter as_tibble arrange mutate select distinct
#' @importFrom pbapply pblapply
#' @importFrom data.table rbindlist
#' @importFrom tidyr separate
nebula_re_rank <-
  function(
           nebula_name,
           top_n = 50,
           match_pattern = NULL, # c("precursorFormula"), or c("precursorFormula", "adduct")
           collate_factor = NA,
           only_gather_structure = F,
           ## ------------------------------------- 
           reference_compound = NA,
           reference_ratio = 0.5,
           ## ------------------------------------- 
           cluster_method = NA,
           csi_score_weight = 0.6,
           class_similarity_weight = 0.3,
           ## ------------------------------------- 
           filter_via_classification = F,
           ## ------------------------------------- 
           rt_set = NA,
           rt_weight = 0.1,
           rt_window = 1.5,
           ## ------------------------------------- 
           revise_MCn_formula_set = F,
           revise_MCn_structure_set = F,
           ...
           ){
    cat("[INFO] MCnebula run: nebula_re_rank\n")
    ## ---------------------------------------------------------------------- 
    structure_set <- get_nebula_structure_candidates(nebula_name, top_n, match_pattern, collate_factor,
                                                     ...)
    if(only_gather_structure == T){
      return(dplyr::as_tibble(structure_set))
    }
    ## ---------------------------------------------------------------------- 
    if(is.data.frame(reference_compound) == F){
      reference_compound <- set_reference_compound(structure_set, reference_ratio)
    }else{
      reference_compound <- reference_compound
    }
    ## ---------------------------------------------------------------------- 
    ## cluster method
    if(is.na(cluster_method) == F){
      method_fun <- match.fun(cluster_method)
      cat("## netbula_re_rank:", paste0(cluster_method), "\n")
      structure_set <- method_fun(structure_set, reference_compound, csi_score_weight = csi_score_weight,
                                  class_similarity_weight = class_similarity_weight,
                                  ...)
    }
    ## ---------------------------------------------------------------------- 
    ## retrive class of candidates via classyfire
    if(filter_via_classification == T){
      cat("## netbula_re_rank: method_filter_candidates_upon_classyfire\n")
      structure_set <- method_filter_candidates_upon_classyfire(structure_set, nebula_name, ...)
    }
    ## ---------------------------------------------------------------------- 
    ## rt prediction
    if(is.data.frame(rt_set)){
      cat("## netbula_re_rank: method_predict_candidates_rt\n")
      structure_set <- method_predict_candidates_rt(structure_set, reference_compound, rt_set,
                                                    rt_weight = rt_weight, rt_window = rt_window, ...)
    }
    ## ---------------------------------------------------------------------- 
    ## revise .GlobalVar .MCn.formula_set
    if(revise_MCn_formula_set == T){
      revise_MCn_formula_set(structure_set)
    }
    ## ---------------------------------------------------------------------- 
    ## revise .GlobalVar .MCn.structure_set -------
    if(revise_MCn_structure_set == T){
      revise_MCn_structure_set(structure_set)
    }
    ## ---------------------------------------------------------------------- 
    cat("[INFO] MCnebula Job Done: nebula_re_rank\n")
    return(dplyr::as_tibble(structure_set))
  }
## ---------------------------------------------------------------------- 
smiles_to_sdfset <-
  function(
           structure_set
           ){
    ##
    smiles_set <- structure_set$smiles
    names(smiles_set) <- paste0(structure_set$".id", "_", structure_set$structure_rank)
    ## this function automaticly set the vector name as name of each subset
    sdf_set <- ChemmineR::smiles2sdf(smiles_set)
    return(sdf_set)
  }
df_get_structure <-
  function(
           x,
           top_n = 10,
           collate_factor = 0.85,
           ...
           ){
    df <- get_structure(
                        x[[".id"]],
                        x[["precursorFormula"]],
                        x[["adduct"]],
                        return_row = 1:top_n,
                        ...)
    if(nrow(df) == 0){
      return(df)
    }
    df <- dplyr::mutate(df, .id = x[[".id"]]) ## add key_id
    if(is.na(collate_factor) == F){
      top_simi <- df[1, "tanimotoSimilarity"]
      df <- dplyr::filter(df, tanimotoSimilarity >= top_simi * collate_factor)
    }
    return(df)
  }
rename_file <-
  function(
           file,
           suffix = "prefix"
           ){
    if(file.exists(file) == T){
      file.rename(file, paste0(file, ".", suffix))
    }
  }
revise_MCn_formula_set <- 
  function(
           structure_set
           ){
    ## prepare replace data
    rp <- dplyr::arrange(structure_set, .id) %>%
      tidyr::separate(col = "file_name", sep = "_", into = c("precursorFormula", "adduct")) %>%
      dplyr::mutate(adduct = gsub("\\+(?!$)", " \\+ ", adduct, perl = T),
                    adduct = gsub("\\-(?!$)", " \\- ", adduct, perl = T)) %>%
      dplyr::select(.id, precursorFormula, adduct, molecularFormula)
    ## replace
    fset <- dplyr::arrange(.MCn.formula_set, .id)
    fset[fset$".id" %in% rp$".id", c(".id", "precursorFormula", "adduct", "molecularFormula")] <- rp
    .MCn.formula_set <<- fset
    return()
  }
revise_MCn_structure_set <- 
  function(
           structure_set
           ){
    sset <- dplyr::arrange(.MCn.structure_set, .id)
    ## prepare replace data
    rp <- dplyr::arrange(structure_set, .id) %>%
      dplyr::select(colnames(sset))
    ## replace
    sset <- dplyr::distinct(rbind(rp, sset), .id, .keep_all = T)
    .MCn.structure_set <<- sset
    ## rename exist structure picture -------
    tmp_stru <- paste0(.MCn.output, "/", .MCn.results, "/tmp/structure")
    if(file.exists(tmp_stru) == T){
      lapply(paste0(tmp_stru, "/", rp$".id", ".svg"), rename_file)
    }
  }
get_nebula_structure_candidates <- 
  function(
           nebula_name,
           top_n = 50,
           match_pattern = NULL,
           collate_factor = NA,
           ... 
           ){
    ## get formula
    id_set <- dplyr::filter(.MCn.nebula_index, name == nebula_name)
    formula_adduct <- dplyr::filter(.MCn.formula_set, .id %in% id_set$".id")
    ## ---------------------------------------------------------------------- 
    ## match patern
    if("precursorFormula" %in% match_pattern == F){
      formula_adduct$precursorFormula = NULL
    }
    if("adduct" %in% match_pattern == F){
      formula_adduct$adduct = NULL
    }
    ## ---------------------------------------------------------------------- 
    ## catch file
    formula_adduct <- by_group_as_list(formula_adduct, ".id")
    ## then, use lapply match file
    cat("## netbula_re_rank: get_structure\n")
    structure_set <- pbapply::pblapply(formula_adduct, df_get_structure,
                                       top_n = top_n,
                                       collate_factor = collate_factor,
                                       ...)
    structure_set <- data.table::rbindlist(structure_set, fill = T)
    cat("## STAT of structure_set:",
        paste0(nrow(structure_set), " (structure sum)/", length(unique(structure_set$".id")), "(.id sum)"), "\n")
    return(structure_set)
  }
set_reference_compound <- 
  function(
           structure_set,
           reference_ratio = 0.5
           ){
    reference_compound <- dplyr::filter(structure_set, structure_rank == 1) %>% 
      dplyr::select(.id, structure_rank, tanimotoSimilarity) %>% 
      dplyr::arrange(tanimotoSimilarity) %>% 
      dplyr::slice(1:(round(reference_ratio * nrow(.))))
    return(reference_compound)
  }
