#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param key_id PARAM_DESCRIPTION
#' @param precursor_formula PARAM_DESCRIPTION, Default: NULL
#' @param adduct PARAM_DESCRIPTION, Default: NULL
#' @param structure_method PARAM_DESCRIPTION, Default: 'top_score'
#' @param order PARAM_DESCRIPTION, Default: T
#' @param return_row PARAM_DESCRIPTION, Default: 1:10
#' @param path PARAM_DESCRIPTION, Default: .MCn.sirius
#' @param as_tibble PARAM_DESCRIPTION, Default: F
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
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[dplyr]{reexports}}
#' @rdname get_structure
#' @export 
#' @importFrom data.table rbindlist
#' @importFrom dplyr as_tibble
get_structure <-
  function(
           key_id,
           precursor_formula = NULL,
           adduct = NULL,
           structure_method = "top_score", # or top_similarity
           order = T,
           return_row = 1:10, # or "all"
           path = .MCn.sirius,
           as_tibble = F,
           ...
           ){
    path <- list.files(path = path, pattern = paste0("*_", key_id, "$"), full.names = T)
    files <- list.files(paste0(path, "/fingerid"),
                        pattern = paste0(precursor_formula, "(.*)", escape_ch(adduct), "(.*).tsv$"),
                        full.names = F)
    ## ---------------------------------------------------
    ## read file
    list <- lapply(paste0(path, "/fingerid/", files), read_tsv)
    names(list) <- gsub(".tsv", "", files)
    ## ---------------------------------------------------
    ## reformat the data
    if(order == T){
      ## bind row as data frame
      df <- data.table::rbindlist(list, idcol = T)
      colnames(df)[which(colnames(df) == ".id")] <- "file_name"
      if(nrow(df) == 0){
        return(df)
      }
      ## ---------------------------------------------------
      ## order upon CSI:fingerID score
      if(structure_method == "top_score"){
        df <- df[order(-df$score),]
        df$structure_rank = as.numeric(1:nrow(df))
      ## order upon tanimoto similarity
      }else if(structure_method == "top_similarity"){
        df <- df[order(-df$similarity),]
        df$structure_rank = as.numeric(1:nrow(df))
      }
      ## ---------------------------------------------------
      if(as_tibble == T){
        df <- dplyr::as_tibble(df)
      }
      # return with top n
      if(return_row[1] != "all"){
        return(df[which(1:nrow(df) %in% return_row),])
      }else{
        return(df)
      }
    }
    return(list)
}
