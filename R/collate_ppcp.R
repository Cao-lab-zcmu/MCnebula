#' @title collate_ppcp
#' @description Collate posterior probability of classification prediction (PPCP) from SIRIUS project
#' and conduct integration to get nebula_class and nebula-index.
#' @param dirs Vector, Default: 'all'
#' @param write_output Logic, Default: T
#' @param nebula_class Logic, Default: T
#' @param nebula_index Logic, Default: T
#' @param ... ...
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[dplyr]{rename}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{filter}}, \code{\link[dplyr]{reexports}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname collate_ppcp
#' @export 
#' @importFrom pbapply pbsapply pblapply
#' @importFrom dplyr rename mutate filter as_tibble
#' @importFrom data.table rbindlist
collate_ppcp <- 
  function(
           dirs = "all",
           write_output = T,
           nebula_class = T,
           nebula_index = T,
           ...
           ){
    cat("[INFO] MCnebula run: collate_ppcp\n")
    ## ---------------------------------------------------------------------- 
    ## check dirs ---- canopus
    cat("## collate_ppcp: check_dir\n")
    if(length(dirs) == 1 & dirs[1] == "all"){
      dirs <- list.files(path = .MCn.sirius, pattern = "^[0-9](.*)_(.*)_(.*)$", full.names = F)
      check <- pbapply::pbsapply(dirs, check_dir, file = "canopus") %>% unname
    }else{
      check <- pbapply::pbsapply(dirs, check_dir, file = "canopus") %>% unname
    }
    ## ---------------------------------------------------------------------- 
    ## lock on file location
    meta_dir <- dirs[which(check == T)] %>%
      data.frame() %>%
      dplyr::rename(dir = ".") %>%
      dplyr::mutate(.id = sapply(dir, grep_id)) %>%
      merge(.MCn.formula_set, by = ".id", all.x = T, sort = F) %>%
      dplyr::mutate(adduct_trans = gsub(" ", "", adduct),
             target = paste0(precursorFormula, "_", adduct_trans, ".fpt"), 
             full.name = paste0(.MCn.sirius, "/", dir, "/", "canopus", "/", target), 
             ## these files need to be check and filter (whether exist)
             ## note that some formula is no fingerprint computed
             ppcp = file.exists(full.name))
      ## ------------------------------------- 
    meta_dir_filter <- dplyr::filter(meta_dir, ppcp == T)
    cat("## STAT of PPCP dataset:",
        paste0(nrow(meta_dir_filter), "(formula with PPCP)", "/", nrow(meta_dir), "(all formula)"), 
        "\n")
    ## ---------------------------------------------------------------------- 
    ## load all ppcp dataset
    if(!exists(".MCn.ppcp_dataset")){
      ppcp_dataset <- pbapply::pblapply(meta_dir_filter$full.name, read_fpt)
    }else{
      ppcp_dataset <- .MCn.ppcp_dataset
    }
    names(ppcp_dataset) <- meta_dir_filter$".id"
    .MCn.ppcp_dataset <<- ppcp_dataset
    ## ---------------------------------------------------------------------- 
    ## summarize nebula_class
    if(nebula_class){
      cat("## collate_ppcp: method_summarize_nebula_class\n")
      metadata <- data.table::rbindlist(.MCn.class_tree_list, idcol = T) %>%
        dplyr::rename(hierarchy = .id)
      ## ------------------------------------- 
      ## set as global var
      .MCn.class_tree_data <<- dplyr::as_tibble(metadata)
      ## ------------------------------------- 
      assign("envir_meta", environment(), envir = parent.env(environment()))
      ## get nebula classes
      nebula_class <- pbapply::pblapply(ppcp_dataset, method_summarize_nebula_class, 
                             class_data_type = "classes_tree_data",
                             ...)
      ## ------------------------------------- 
      .MCn.nebula_class <<- nebula_class
    }
    ## ---------------------------------------------------------------------- 
    if(nebula_index){
      cat("## collate_ppcp: method_summarize_nebula_index.\n")
      ## gather all nebula classes
      nebula_index <- method_summarize_nebula_index(ppcp_dataset,
                                                    ...)
      .MCn.nebula_index <<- nebula_index
    ## ---------------------------------------------------------------------- 
      if(write_output){
        output = paste0(.MCn.output, "/", .MCn.results)
        write_tsv(nebula_index, file = paste0(output, "/", "nebula_index.tsv"))
      }
    }
    ## ---------------------------------------------------------------------- 
    rm("envir_meta", envir = parent.env(environment()))
    ## ---------------------------------------------------------------------- 
    cat("[INFO] MCnebula Job Done: collate_ppcp.\n")
    return(nebula_index)
  }
