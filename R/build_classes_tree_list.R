#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param class_index PARAM_DESCRIPTION, Default: 'canopus.tsv'
#' @param path PARAM_DESCRIPTION, Default: .MCn.sirius
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{reexports}}
#' @rdname build_classes_tree_list
#' @export 
#' @importFrom dplyr as_tibble
build_classes_tree_list <-
  function(
           class_index="canopus.tsv", path=.MCn.sirius
           ){
    data <- read_tsv(paste0(path, "/", class_index))
    ## ---------------------------------------------------------------------- 
    ## separate each levels of classes into sub-list
    root <- data[which(data$parentId==""), ]
    list <- list()
    n = 1
    list[[n]] <- root %>% dplyr::as_tibble()
    df <- data[data$parentId %in% root$id, ]
    ## ---------------------------------------------------------------------- 
    while(nrow(df) > 0){
      n = n + 1
      list[[n]] <- df %>% dplyr::as_tibble()
      df <- data[data$parentId %in% df$id, ]
    }
    .MCn.class_tree_list <<- list
    cat("INFO: Classification Index in.MCn.sirius project --->", class_index, "\nA total of 11 levels. These classes (upon ClassyFire and CANOPUS) were separated into sub-lists.
        Use following arguments to get some specific classes:
        .MCn.class_tree_list[[3]] >>> superclass
        .MCn.class_tree_list[[4]] >>> class
        .MCn.class_tree_list[[5]] >>> subclass
        .MCn.class_tree_list[[6]] >>> level 5 \n")
}