#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param key_id PARAM_DESCRIPTION, Default: NULL
#' @param dir PARAM_DESCRIPTION, Default: NULL
#' @param precursor_formula PARAM_DESCRIPTION, Default: 'method_pick_formula_excellent'
#' @param adduct PARAM_DESCRIPTION, Default: NULL
#' @param reformat PARAM_DESCRIPTION, Default: T
#' @param filter PARAM_DESCRIPTION, Default: T
#' @param filter_threshold PARAM_DESCRIPTION, Default: 0.1
#' @param class_index PARAM_DESCRIPTION, Default: 'canopus.tsv'
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname get_ppcp
#' @export 
get_ppcp <- 
  function(
           key_id = NULL,
           dir = NULL,
           precursor_formula = "method_pick_formula_excellent",
           adduct = NULL,
           reformat = T,
           filter = T,
           filter_threshold = 0.1,
           class_index = "canopus.tsv",
           ...
           ){
    ## get dir path
    if(is.null(dir) == T & is.null(key_id) == T){
      return()
    }else if(is.null(dir) == T){
      dir <- get_dir(key_id)
    }
    ## ---------------------------------------------------
    ## aquire formula via the method
    if( precursor_formula == "method_pick_formula_excellent" ){
      meta <- method_pick_formula_excellent(dir = dir)
      precursor_formula <- meta$precursorFormula
      adduct <- meta$adduct
    }
    ## ---------------------------------------------------
    ## read ppcp data
    file <- list.files(path = paste0(.MCn.sirius, "/", dir, "/", "canopus"),
                       pattern = paste0("^", precursor_formula, "(.*)", escape_ch(adduct), "(.*)", ".fpt$"),
                       full.names = T)
    ppcp <- read_fpt(file)
    ## ---------------------------------------------------
    ## reformat section
    if(reformat == F){
      return(ppcp)
    }
    ## check meta list
    if(exists(".MCn.class_tree_list") == F){
      build_classes_tree_list(class_index = class_index)
    }
    ## merge with meta table, and filter
    ppcp <- lapply(.MCn.class_tree_list, merge_class_ppcp,
                   ## parameter
                   key_id = key_id,
                   values = ppcp,
                   filter = filter,
                   filter_threshold = filter_threshold)
    return(ppcp)
  }
## a small function to get data of ppcp
read_fpt <- function(file){
  fpt = data.table::fread(input = file, header = F, quote = "")
  fpt$relativeIndex = seq(0, nrow(fpt) - 1)
  return(fpt)
}
## specific character in adduct description need to be revise, for pattern matching
escape_ch <- function(x){
  x <- gsub("\\[", "\\\\\\[", x)
  x <- gsub("\\]", "\\\\\\]", x)
  x <- gsub("\\+", "\\\\\\+", x)
  x <- gsub("\\-", "\\\\\\-", x)
  x <- gsub(" ", "", x)
  return(x)
}
## the function to merge raw ppcp with meta list
merge_class_ppcp <-
  function(
           class,
           values,
           filter = T,
           filter_threshold = 0.1,
           key_id = NULL,
           filter_col = "V1"
           ){
    df <- merge(class, values, all.x = T, by = "relativeIndex", sort = F)
    df <- df[which(df[[filter_col]] > ifelse(filter == T, filter_threshold, 0)),] %>%
      dplyr::as_tibble()
    return(df)
  }
