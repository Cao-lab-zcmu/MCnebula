#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nebula_name PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname visualize_nebula
#' @export 
visualize_nebula <-
  function(
           nebula_name,
           ...
           ){
    p <- annotate_child_nebulae(
                                nebula_name = nebula_name,
                                write_output = F,
                                plot_structure = F,
                                plot_ppcp = F,
                                merge_image = F,
                                return_plot = T,
                                ...)
    return(p)
  }
