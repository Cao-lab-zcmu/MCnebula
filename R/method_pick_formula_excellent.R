#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dir PARAM_DESCRIPTION, Default: NULL
#' @param key_id PARAM_DESCRIPTION, Default: NULL
#' @param exclude_element PARAM_DESCRIPTION, Default: NULL
#' @param ppm_error PARAM_DESCRIPTION, Default: 20
#' @param fc PARAM_DESCRIPTION, Default: 1.5
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}
#' @rdname method_pick_formula_excellent
#' @export 
#' @importFrom dplyr mutate
method_pick_formula_excellent <- 
  function(
           key_id = NULL, 
           dir = NULL,
           exclude_element = NULL,
           ppm_error = 20,
           fc = 1.5
           ){
    ## ----------------------------------------------------------------------
    if(length(dir) >= 1){
      meta <- data.table::data.table(dir = dir)
      meta <- dplyr::mutate(meta, key_id = lapply(dir, grep_id))
    }else{
      meta <- data.table::data.table(key_id = key_id, dir = get_dir(key_id))
    }
    ## ---------------------------------------------------------------------- 
    meta <- dplyr::mutate(meta, formula_dir = paste0(.MCn.sirius, "/", dir, "/", "formula_candidates.tsv"),
                          structure_dir = paste0(.MCn.sirius, "/", dir, "/", "structure_candidates.tsv"),
                          fingerid = paste0(.MCn.sirius, "/", dir, "/", "fingerid"))
    ## ------------------------------------- 
    cat("## Method part: batch_get_formula\n")
    formula_list <- pbapply::pblapply(meta$key_id, mutate_get_formula,
                                      ppm_error = ppm_error, exclude_element = exclude_element)
    formula_df <- data.table::rbindlist(formula_list)
    ## ------------------------------------- 
    ## get top ZodiacScore within a key_id formula set
    df_fz <- lapply(formula_list, head, n = 1)
    df_fz <- data.table::rbindlist(df_fz)
    ## ------------------------------------- 
    cat("## Method part: batch_get_structure\n")
    structure_list <- pbapply::pbmapply(mutate_get_structure, meta$structure_dir, meta$key_id,
                                        SIMPLIFY = F)
    structure_df <- data.table::rbindlist(structure_list)
    ## -------------------------------------
    ## merge structure_df with formula_df, to get ZodiacScore of top structure
    structure_df <- merge(structure_df, formula_df, by = c(".id", "molecularFormula", "adduct"))
    ## ---------------------------------------------------------------------- 
    ## ---------------------------------------------------------------------- 
    if(is.na(fc) == F){
      df_sz <- dplyr::rename(structure_df, sz_score = ZodiacScore)
      ## ------------------------------------- 
      ## compare fz_score with sz_score
      compare <- merge(df_fz[, c(".id", "ZodiacScore"), with = F],
                       df_sz[, c(".id", "sz_score"), with = F],
                       all.x = T, by = ".id")
      compare <- dplyr::mutate(compare,
                               use_zodiac = ifelse(is.na(sz_score), T,
                                                   ifelse(ZodiacScore >= sz_score * fc, T, F)
                                                   ))
      ## -------------------------------------
      ## the .id which formula of use_zodiac or not
      fz <- dplyr::filter(compare, use_zodiac == T)$".id"
      sz <- dplyr::filter(compare, use_zodiac == F)$".id"
    }else{
      sz <- structure_df$".id"
      fz <- df_fz$".id"[which(!df_fz$".id" %in% sz)]
    }
    ## ----------------------------------------------------------------------
    ## ----------------------------------------------------------------------
    df_fz <- mutate(df_fz[.id %in% fz, ], use_zodiac = T)
    df_sz <- mutate(structure_df[.id %in% sz, ], use_zodiac = F)
    ## ------------------ 
    formula_adduct <- dplyr::bind_rows(df_fz, df_sz)
    formula_adduct <- dplyr::as_tibble(formula_adduct) %>%
      dplyr::select(.id, colnames(.))
    ## ----------------------------------------------------------------------
    return(formula_adduct)
  }
mutate_get_formula <- 
  function(
           key_id,
           ppm_error,
           exclude_element
           ){
    formula_df <- try(silent = T, get_formula(key_id, rank = "all", ppm_error = ppm_error,
                                              exclude_element = exclude_element))
    if(class(formula_df)[1] == "try-error"){
      return()
    }else{
      formula_df <- dplyr::mutate(formula_df,
                                  ZodiacScore = ifelse(grepl("[0-9]", ZodiacScore),
                                                       ZodiacScore, 0),
                                  ZodiacScore = as.numeric(ZodiacScore))
      formula_df$".id" <- key_id
      return(formula_df)
    }
  }
mutate_get_structure <-
  function(
           structure_dir,
           key_id
           ){
    structure_df <- try(silent = T, read_tsv(structure_dir))
    if(class(structure_df)[1] == "try-error"){
      return()
    }
    if(nrow(structure_df) == 0){
      return()
    }
    max <- max(structure_df$"CSI:FingerIDScore")
    structure_df <- structure_df[`CSI:FingerIDScore` == max, c("molecularFormula", "adduct"), with = F]
    if(nrow(structure_df) > 1){
      structure_df <- head(structure_df, n = 1)
    }
    structure_df$".id" <- key_id
    return(structure_df)
  }
get_dir <- function(
                    key_id,
                    path = .MCn.sirius
                    ){
  dir <- list.files(path = path,
                    pattern=paste0("^[0-9](.*)_(.*)_", key_id, "$"),
                    full.names = F)
  check <- check_dir(dir)
  if(check == T){
    return(dir)
  }
}
