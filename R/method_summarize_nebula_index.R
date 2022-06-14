#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ppcp_dataset PARAM_DESCRIPTION
#' @param nebula_class PARAM_DESCRIPTION, Default: .MCn.nebula_class
#' @param ppcp_threshold PARAM_DESCRIPTION, Default: 0.5
#' @param min_possess PARAM_DESCRIPTION, Default: 10
#' @param max_possess_pct PARAM_DESCRIPTION, Default: 0.2
#' @param identical_factor PARAM_DESCRIPTION, Default: 0.8
#' @param filter_identical PARAM_DESCRIPTION, Default: c(top_hierarchy = 4)
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
           ppcp_dataset = .MCn.ppcp_dataset,
           nebula_class = .MCn.nebula_class,
           ppcp_threshold = 0.5,
           # min number of compounds allowed exist in a child-nebula. if less than, filter the nebula
           min_possess = 10, 
           ## max percentage of compounds allowed exist in a child-nebula
           max_possess_pct = 0.2, 
           ## identical filter
           filter_identical = c("top_hierarchy" = 4), 
           identical_factor = 0.7,
           rm_position_describe_class = T,
           ## in the nebula, if too many structure score is too low, filter the nebula.
           ## or NA
           filter_via_struc_score = "tanimotoSimilarity", 
           struc_score_cutoff = 0.3,
           min_reached_pct = 0.6,
           target_classes = NA,
           ...
           ){
    if(is.list(nebula_class)){
      classes <- data.table::rbindlist(nebula_class)
      if(rm_position_describe_class){
        classes <- dplyr::filter(classes, !grepl("[0-9]", name))
      }
      ## classes is merely a classes index number set
      classes <- dplyr::distinct(classes, relativeIndex)$relativeIndex
    }else{
      classes <- c()
    }
    ## user defined classes
    if(is.vector(target_classes)){
      target_classes <- dplyr::filter(.MCn.class_tree_data, name %in% target_classes)
      ## merge
      classes <- c(classes, target_classes$relativeIndex)
    }
    ## get classes
    cat("## Method part: class_retrieve\n")
    index_list <- pbapply::pblapply(ppcp_dataset, class_retrieve,
                                    the_relativeIndex = classes,
                                    ...)
    index_df <- data.table::rbindlist(index_list, idcol = T)
    ## ---------------------------------------------------------------------- 
    ## filter via max_possess and min_possess
    stat <- table(index_df$relativeIndex)
    stat <- stat[which(stat >= min_possess & stat <= max_possess_pct * length(unique(index_df$".id")))]
    index_df <- index_df %>%
      dplyr::filter(relativeIndex %in% names(stat))
    ## ------------------------------------- 
    ## gather with classes annotation
    index_df <- data.table::rbindlist(.MCn.class_tree_list, idcol = T) %>%
      dplyr::rename(hierarchy = .id) %>%
      dplyr::select(relativeIndex, name, hierarchy) %>%
      merge(index_df, by = "relativeIndex", all.y = T, sort = F) %>%
      data.table::data.table()
    ## ---------------------------------------------------------------------- 
    if(!is.na(identical_factor)){
      ## filter identical or similar classes
      ## enumerate combination
      class_for_merge <- index_df %>% 
        dplyr::filter(hierarchy >= filter_identical[["top_hierarchy"]]) %>%
        dplyr::distinct(relativeIndex) %>%
        unlist() %>% combn(m = 2) %>%
        t() %>% data.frame()
      ## ------------------------------------- 
      cat("## Method part: identical_filter\n")
      discard = pbapply::pbapply(class_for_merge, 1, identical_filter,
                                 index_df = index_df,
                                 identical_factor = identical_factor,
                                 ...) %>%
        unlist() %>% unique()
      index_df <- dplyr::filter(index_df, !relativeIndex %in% discard)
    }
    ## ---------------------------------------------------------------------- 
    if(!is.na(filter_via_struc_score)){
      df <- merge(index_df, .MCn.structure_set[, c(".id", filter_via_struc_score)],
                  by = ".id", all.x = T)
      list <- by_group_as_list(df, "relativeIndex")
      ## ------------------------------------- 
      cat("## Method part: fun_filter_via_struc_score\n")
      select_index <- pbapply::pblapply(list, fun_filter_via_struc_score,
                                        filter_via_struc_score,
                                        struc_score_cutoff,
                                        min_reached_pct)
      select_index <- unlist(select_index)
      index_df <- dplyr::filter(index_df, relativeIndex %in% all_of(select_index))
    }
    ## ---------------------------------------------------------------------- 
    ## cluster id in each classes
    nebula_index <- dplyr::group_by(index_df, relativeIndex)
    return(nebula_index)
  }
class_retrieve <-
  function(
           data,
           the_relativeIndex,
           ppcp_threshold = 0.5,
           ...
           ){
    ##
    classes <- the_relativeIndex
    data <- dplyr::filter(data, relativeIndex %in% classes, V1 >= ppcp_threshold)
    return(data)
  }
identical_filter <- 
  function(
           couple,
           index_df,
           identical_factor = 0.7,
           ...
           ){
    ## index_df is a data.table project
    x = unique(index_df[relativeIndex %in% couple[1], ]$".id")
    y = unique(index_df[relativeIndex %in% couple[2], ]$".id")
    p_x = table(x %in% y)
    p_y = table(y %in% x)
    ## ------------------------------------- 
    if("TRUE" %in% names(p_x) == F | "TRUE" %in% names(p_y) == F){
      return()
    }
    p_x = prop.table(p_x)[["TRUE"]]
    p_y = prop.table(p_y)[["TRUE"]]
    ## ------------------------------------- 
    if(p_x >= identical_factor & p_y >= identical_factor){
      idn = ifelse(length(x) >= length(y), couple[2], couple[1])
      return(idn)
    }else{
      return()
    }
  }
fun_filter_via_struc_score <- 
  function(
           df,
           score = "tanimotoSimilarity",
           cutoff = 0.4,
           min_reached_pct = 0.5
           ){
    x <- df[[score]]
    df <- dplyr::mutate(df, reach = ifelse(x >= cutoff &
                                           is.na(x) == F,
                                         T, F))
    check <- prop.table(table(df[["reach"]]))
    if("TRUE" %in% names(check) == F)
      return()
    if(check[["TRUE"]] >= min_reached_pct){
      return(df[1, ]$relativeIndex)
    }else{
      return()
    }
  }
## ---------------------------------------------------------------------- 
method_summarize_target_index <- 
  function(
           target_classes
           ){
    target_index <- method_summarize_nebula_index(nebula_class = NA,
                                                  target_classes = target_classes,
                                                  identical_factor = NA,
                                                  filter_via_struc_score = NA,
                                                  max_possess_pct = 1,
                                                  min_possess = 1
    )
    return(target_index)
  }
