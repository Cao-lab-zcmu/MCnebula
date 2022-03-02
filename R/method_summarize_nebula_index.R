#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ppcp_dataset PARAM_DESCRIPTION
#' @param nebula_class PARAM_DESCRIPTION, Default: .MCn.nebula_class
#' @param ppcp_threshold PARAM_DESCRIPTION, Default: 0.5
#' @param min_possess PARAM_DESCRIPTION, Default: 10
#' @param max_possess_pct PARAM_DESCRIPTION, Default: 0.2
#' @param identical_merge PARAM_DESCRIPTION, Default: T
#' @param identical_factor PARAM_DESCRIPTION, Default: 0.8
#' @param merge_allowed_hierarchy PARAM_DESCRIPTION, Default: c(top_level = 4)
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
#'  \code{\link[dplyr]{distinct}}, \code{\link[dplyr]{filter}}, \code{\link[dplyr]{rename}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{group_by}}
#'  \code{\link[pbapply]{pbapply}}
#' @rdname method_summarize_nebula_index
#' @export 
#' @importFrom data.table rbindlist
#' @importFrom dplyr distinct filter rename select group_by
#' @importFrom pbapply pblapply
method_summarize_nebula_index <-
  function(
           ppcp_dataset,
           nebula_class = .MCn.nebula_class,
           ppcp_threshold = 0.5,
           # min number of compounds allowed exist in a child-nebula. if less than, filter the nebula
           min_possess = 10, 
           ## max percentage of compounds allowed exist in a child-nebula
           max_possess_pct = 0.2, 
           identical_merge = T,
           identical_factor = 0.8,
           ## if set to 4 (class), the level of or below this hierarchy (e.g., subclass) will perform merge
           merge_allowed_hierarchy = c("top_level" = 4), 
           ...
           ){
    classes <- data.table::rbindlist(nebula_class) %>%
      dplyr::distinct(relativeIndex)
    ## get classes
    classes <- classes$relativeIndex
    ## environment for lapply function
    assign("envir_classes", environment(), envir = parent.env(environment()))
    cat("## Method part: class_retrieve\n")
    index_list <- pbapply::pblapply(ppcp_dataset, class_retrieve,
                         ...)
    index_df <- data.table::rbindlist(index_list, idcol = T)
    ## ---------------------------------------------------------------------- 
    ## filter via max_possess and min_possess
    stat <- table(index_df$relativeIndex)
    stat <- stat[which(stat >= min_possess & stat <= max_possess_pct * length(unique(index_df$".id")))]
    index_df <- index_df %>%
      dplyr::filter(relativeIndex %in% names(stat))
    ## gather with classes annotation
    index_df <- data.table::rbindlist(.MCn.class_tree_list, idcol = T) %>%
      dplyr::rename(hierarchy = .id) %>%
      dplyr::select(relativeIndex, name, hierarchy) %>%
      merge(index_df, by = "relativeIndex", all.y = T, sort = F)
    ## ---------------------------------------------------------------------- 
    if(identical_merge == T){
      ## filter identical or similar classes
      ## enumerate combination
      class_for_merge <- index_df %>% 
        dplyr::filter(hierarchy >= merge_allowed_hierarchy[["top_level"]]) %>%
        dplyr::distinct(relativeIndex) %>%
        unlist() %>%
        combn(m = 2) %>%
        t() %>%
        data.frame()
      cat("## Method part: identical_filter\n")
      discard = pbapply(class_for_merge, 1, identical_filter,
                      identical_factor = identical_factor,
                      ...) %>%
        unlist() %>%
        unique()
    }
    index_df <- index_df[!index_df$relativeIndex %in% discard, ]
    ## cluster id in each classes
    nebula_index <- dplyr::group_by(index_df, relativeIndex)
    return(nebula_index)
  }
class_retrieve <-
  function(
           data,
           the_relativeIndex = get("classes", envir = get("envir_classes")),
           ppcp_threshold = 0.5
           ){
    ##
    classes <- the_relativeIndex
    data <- dplyr::filter(data, relativeIndex %in% classes, V1 >= ppcp_threshold)
    return(data)
  }
identical_filter <- 
  function(
           couple,
           index_df = get("index_df", envir = get("envir_classes")),
           identical_factor = 0.7
           ){
    ##
    x = unique(index_df[which(index_df$relativeIndex %in% couple[1]), ]$".id")
    y = unique(index_df[which(index_df$relativeIndex %in% couple[2]), ]$".id")
    p_x = table(x %in% y)
    p_y = table(y %in% x)
    if("TRUE" %in% names(p_x) == F | "TRUE" %in% names(p_y) == F){
      return()
    }
    p_x = prop.table(p_x)[["TRUE"]]
    p_y = prop.table(p_y)[["TRUE"]]
    if(p_x >= identical_factor & p_y >= identical_factor){
      idn = ifelse(length(x) >= length(y), couple[2], couple[1])
      return(idn)
    }else{
      return()
    }
  }
