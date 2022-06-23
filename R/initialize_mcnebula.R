#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sirius_path PARAM_DESCRIPTION
#' @param output_path PARAM_DESCRIPTION, Default: sirius_path
#' @param output_file PARAM_DESCRIPTION, Default: 'mcnebula_results'
#' @param palette PARAM_DESCRIPTION, Default: unique(c((ggsci::pal_simpsons())(16), (ggsci::pal_igv("default"))(51)))
#' @param palette_stat PARAM_DESCRIPTION, Default: palette
#' @param palette_label PARAM_DESCRIPTION, Default: colorRampPalette(c("#C6DBEFFF", "#3182BDFF", "red"))(10)
#' @param rm_var PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[ggsci]{pal_simpsons}}, \code{\link[ggsci]{pal_igv}}
#' @rdname initialize_mcnebula
#' @export 
#' @importFrom ggsci pal_simpsons pal_igv

initialize_mcnebula <-
  function(
           sirius_path,
           output_path = sirius_path,
           output_file = "mcnebula_results",
           palette = unique(c(ggsci::pal_simpsons()(16),
                              ggsci::pal_igv("default")(51),
                              ggsci::pal_ucscgb()(26),
                              ggsci::pal_d3("category20")(20)
                              )),
           palette_stat = palette,
           palette_ppcp = palette,
           palette_label = colorRampPalette(c("#C6DBEFFF", "#3182BDFF", "red"))(10),
           rm_mc.set = F
           ){
    if(rm_mc.set){
      rm_mc.set(envir = .GlobalEnv)
    }
    if(!file.exists(sirius_path) | !file.exists(output_path)){
      cat("File path not find.\n")
      return()
    }
    if(!file.exists(paste0(sirius_path, "/", ".format"))){
      cat("SIRIUS project not find.\n")
      return()
    }
    .MCn.sirius <<- sirius_path
    .MCn.output <<- output_path
    .MCn.results <<- output_file
    .MCn.palette <<- palette
    .MCn.palette_stat <<- palette_stat
    .MCn.palette_ppcp <<- palette_ppcp
    .MCn.palette_label <<- palette_label
    dir.create(paste0(.MCn.output, "/", .MCn.results))
    cat("MCnebula project has initialized at ->", .MCn.output, "\n")
  }
## ---------------------------------------------------------------------- 
rm_mc.set <- 
  function(
           envir,
           pattern = "^\\.MCn\\..*"
           ){
    mc.set <- ls(pattern = pattern, envir = envir, all.names = T)
    rm(list = mc.set, envir = envir)
  }
