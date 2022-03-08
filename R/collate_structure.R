#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dirs PARAM_DESCRIPTION, Default: 'all'
#' @param path PARAM_DESCRIPTION, Default: .MCn.sirius
#' @param output PARAM_DESCRIPTION, Default: paste0(.MCn.output, "/", .MCn.results)
#' @param write_output PARAM_DESCRIPTION, Default: T
#' @param write_picked_formula_adduct PARAM_DESCRIPTION, Default: T
#' @param collate_method PARAM_DESCRIPTION, Default: 'method_pick_formula_excellent'
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
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{reexports}}, \code{\link[dplyr]{filter}}
#' @rdname collate_structure
#' @export 
#' @importFrom pbapply pbsapply pblapply
#' @importFrom data.table rbindlist
#' @importFrom dplyr mutate as_tibble filter
collate_structure <- 
  function(
           dirs = "all",
           path = .MCn.sirius,
           output = paste0(.MCn.output, "/", .MCn.results),
           write_output = T,
           write_picked_formula_adduct = T,
           collate_method = "method_pick_formula_excellent", # "top_score", "top_similarity", "top_zodiac_score"
           ...
           ){
  cat( paste0("[INFO] MCnebula run: collate_structure\n") )
  ## -----------------------------------------------------------------
  ## check dirs
  cat("## collate_structure: check_dir\n")
  if(dirs == "all"){
    dirs <- list.files(path = path, pattern="^[0-9](.*)_(.*)_(.*)$", full.names = F)
    check <- pbapply::pbsapply(dirs, check_dir) %>% unname
  }else{
    check <- pbapply::pbsapply(dirs, check_dir) %>% unname
  }
  dirs <- dirs[which(check == T)]
  ## build a new envir to place data
  formula_cache <- new.env()
  structure_cache <- new.env()
  ## -----------------------------------------------------------------
  cat("## collate_structure:", paste0(collate_method), "\n")
  method_fun <- match.fun(collate_method)
  ## method
  pbapply::pblapply(dirs, method_fun,
                    return_formula = F,
                    ## the data are placed into cache envir
                    formula_cache = formula_cache,
                    structure_cache = structure_cache,
                    ...)
  ## -----------------------------------------------------------------
  ## structure collate
  structure_dataset <- eapply(structure_cache, data.table) 
  structure_dataset <- data.table::rbindlist(structure_dataset, idcol = T)
  .MCn.structure_set <<- dplyr::mutate(structure_dataset,
                                       tanimotoSimilarity = as.numeric(tanimotoSimilarity)) %>%
    dplyr::as_tibble()
  ## write output
  if(write_output == T){
    write_tsv( structure_dataset, paste0(output, "/", collate_method, ".structure.tsv"))
  }
  ## ------------------------------------- 
  ## formula_adduct collate
  formula_adduct_set <- eapply(formula_cache, data.table) %>%
    data.table::rbindlist(idcol = T)
  .MCn.formula_set <<- dplyr::filter(formula_adduct_set, is.na(precursorFormula) == F)
  ## write output
  if(write_picked_formula_adduct == T){
    write_tsv(formula_adduct_set, paste0(output, "/", collate_method, ".tsv"))
  }
  ## -----------------------------------------------------------------
  cat( paste0("[INFO] MCnebula Job Done: collate_structure.\n") )
}
grep_id <- function(x, sep = "_"){
  v <- unlist(strsplit(x, split="_"))
  id <- v[length(v)]
  return(id)
}
## ----
check_dir <- function(dir, path = .MCn.sirius, file = "compound.info"){
  if(file.exists(paste0(path, "/", dir, "/", file)) == T){
    check = T
  }else{
    check = F
  }
  return(check)
}
