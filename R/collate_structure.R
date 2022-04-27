#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dirs PARAM_DESCRIPTION, Default: 'all'
#' @param write_output PARAM_DESCRIPTION, Default: T
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
           write_output = T,
           ...
           ){
  cat( paste0("[INFO] MCnebula run: collate_structure\n") )
  ## ---------------------------------------------------------------------- 
  ## check dirs
  cat("## collate_structure: check_dir\n")
  if(dirs == "all"){
    dirs <- list.files(path = .MCn.sirius, pattern="^[0-9](.*)_(.*)_(.*)$", full.names = F)
    check <- pbapply::pbsapply(dirs, check_dir) %>%
      unname()
  }else{
    check <- pbapply::pbsapply(dirs, check_dir) %>%
      unname()
  }
  dirs <- dirs[which(check == T)]
  ## ---------------------------------------------------------------------- 
  cat("## collate_structure: method_pick_formula_excellent\n")
  formula_set <- method_pick_formula_excellent(dir = dirs, ...)
  ## ------------------------------------- 
  ## set as global var
  .MCn.formula_set <<- formula_set
  ## ---------------------------------------------------------------------- 
  ## structure collate
  cat("## collate_structure: re-collate structure\n")
  structure_dataset <- pbapply::pbmapply(mutate2_get_structure,
                                         formula_set$.id,
                                         formula_set$precursorFormula,
                                         formula_set$adduct,
                                         SIMPLIFY = F)
  structure_dataset <- data.table::rbindlist(structure_dataset, fill = T) %>%
    ## debug
    dplyr::mutate(tanimotoSimilarity = as.numeric(tanimotoSimilarity)) %>% 
    dplyr::relocate(.id) %>% 
    dplyr::as_tibble()
  ## ------------------------------------- 
  ## set as global var
  .MCn.structure_set <<- structure_dataset
  ## ---------------------------------------------------------------------- 
  ## write output
  if(write_output == T){
    output = paste0(.MCn.output, "/", .MCn.results)
    write_tsv(structure_dataset, paste0(output, "/", "method_pick_formula_excellent", ".structure.tsv"))
    write_tsv(formula_set, paste0(output, "/", "method_pick_formula_excellent", ".tsv"))
  }
  ## ---------------------------------------------------------------------- 
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
mutate2_get_structure <- 
  function(
           key_id,
           formula,
           adduct
           ){
    df <- try(silent = T, get_structure(key_id,
                                        precursor_formula = formula,
                                        adduct = adduct,
                                        return_row = 1))
    if(class(df)[1] == "try-error"){
      return()
    }else{
      df$".id" <- key_id
      return(df)
    }
  }
