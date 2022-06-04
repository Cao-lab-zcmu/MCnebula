#' @title get_formula
#' @description A function to read table of formula data of feature in SIRIUS project.
#' @param key_id Character.
#' @param exclude_element Vector, Default: NULL
#' @param formula_method Character, Default: 'top_zodiac'
#' @param rank Vector of number, Default: 1:5
#' @param ppm_error A bumber, Default: 20
#' @param return_col Vector of character, Default: c("rank", "precursorFormula",
# '   "molecularFormula", "adduct", "ZodiacScore", 
#'    "massErrorPrecursor(ppm)")
#' @param ... ...
#' @return A data.frame
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname get_formula
#' @export 
get_formula <-
  function(
           key_id,
           exclude_element = NULL, ## e.g., c("S", "B", "P", "Si")
           formula_method = "top_zodiac",
           rank = 1:5, # or "all"
           ppm_error = 20,
           return_col = c("rank", "precursorFormula", "molecularFormula",
                        "adduct", "ZodiacScore", "massErrorPrecursor(ppm)"),
           ...
           ){
    path <- list.files(path = .MCn.sirius, pattern=paste0("*_", key_id, "$"), full.names=T)
    file <- read_tsv(paste0(path, "/", "formula_candidates.tsv"))
    file$rank <- as.numeric(file$rank)
    ## ---------------------------------------------------------------------- 
    if("ZodiacScore" %in% colnames(file) == F){
      file$ZodiacScore = 0
    }
    ## ---------------------------------------------------------------------- 
    if(is.null(exclude_element) == F){
      file <- file[!unname(sapply(file$precursorFormula, grep_element,
                                  exclude_element = exclude_element)), ]
    }
    ## ---------------------------------------------------------------------- 
    if(formula_method == "top_zodiac"){
      if(rank[1] == "all"){
        rank <- unique(file$rank)
      }
      ## to escape the key in data.table
      filter = rank
      file <- file[abs(file$"massErrorPrecursor(ppm)") <= ppm_error &
                   file$"rank" %in% filter,
                   return_col, with = F]
    }
    return(file)
}
## ---------------------------------------------------------------------- 
grep_element <-
  function(
           formula,
           exclude_element = c("S", "P", "B")
           ){
    check <- grepl(paste(exclude_element, collapse = "|"), formula)
    return(check)
  }
