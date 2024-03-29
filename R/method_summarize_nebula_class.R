#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PPCP dataset
#' @param ppcp_threshold PARAM_DESCRIPTION, Default: 0.5
#' @param max_classes PARAM_DESCRIPTION, Default: 5
#' @param hierarchy_priority PARAM_DESCRIPTION, Default: c(6, 5, 4, 3)
#' @param class_data_type PARAM_DESCRIPTION, Default: 'classes_tree_list'
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
#'  \code{\link[dplyr]{rename}}, \code{\link[dplyr]{filter}}
#' @rdname method_summarize_nebula_class
#' @export 
#' @importFrom data.table rbindlist
#' @importFrom dplyr rename filter
method_summarize_nebula_class <-
  function(
           data,
           ppcp_threshold = 0.5,
           max_classes = NA,
           hierarchy_priority = c(6, 5, 4, 3), ## level 5, subclass, class, superclass
           class_data_type = "classes_tree_list", ## or "classes_tree_data"
           ...
           ){
    ## input data 
    if(class_data_type == "classes_tree_list"){
      class_data = .MCn.class_tree_list
      metadata <- data.table::rbindlist(class_data, idcol = T)
      metadata <- dplyr::rename(metadata, hierarchy = .id)
    }else if(class_data_type == "classes_tree_data"){
      metadata <- class_data <- get("metadata", envir = get("envir_meta"))
    }
    ## main body
    df <- dplyr::filter(data, V1 >= ppcp_threshold)
    df <- merge(df, metadata[, 1:5], all.x = T, by = "relativeIndex", sort = F)
    df <- dplyr::filter(df, hierarchy %in% hierarchy_priority)
    df <- df[order(factor(df$hierarchy, levels = hierarchy_priority), -df$V1), ]
    if(is.na(max_classes) == F)
      df <- head(df, n = max_classes)
    return(df)
  }
