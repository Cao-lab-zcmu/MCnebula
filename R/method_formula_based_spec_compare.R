#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param path PARAM_DESCRIPTION, Default: .MCn.sirius
#' @param dirs PARAM_DESCRIPTION, Default: 'all'
#' @param cpu_cores PARAM_DESCRIPTION, Default: 8
#' @param compare_fun PARAM_DESCRIPTION, Default: 'dotproduct'
#' @param precursor_mass_diff PARAM_DESCRIPTION, Default: T
#' @param edge_filter PARAM_DESCRIPTION, Default: 0.3
#' @param only_identical_class PARAM_DESCRIPTION, Default: T
#' @param min_hierarchy PARAM_DESCRIPTION, Default: 5
#' @param filter_only_max PARAM_DESCRIPTION, Default: 2000
#' @param min_zodiac PARAM_DESCRIPTION, Default: 0.9
#' @param min_tanimoto PARAM_DESCRIPTION, Default: 0.4
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
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[dplyr]{rename}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{filter}}, \code{\link[dplyr]{group_by}}, \code{\link[dplyr]{summarise_all}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{distinct}}, \code{\link[dplyr]{reexports}}
#' @rdname method_formula_based_spec_compare
#' @export 
#' @importFrom pbapply pbsapply pbmapply pbapply pblapply
#' @importFrom dplyr rename mutate filter group_by summarise_at ungroup select distinct as_tibble

method_formula_based_spec_compare <- 
  function(
           path = .MCn.sirius,
           dirs = "all",
           cpu_cores = 8, 
           compare_fun = "dotproduct",
           precursor_mass_diff = T,
           edge_filter = 0.3,
           ## only identical classes will be compared (according to results of function:collate_ppcp
           only_identical_class = T, 
           # only hierarchy >= 5 (class, subclass...) will be considered)
           min_hierarchy = 5, 
           # only nodes number >= 2000, do zoidacScore and tanimotoSimilarity filter
           filter_only_max = 2000, 
           # only ZodiacScore >= 0.5 ...
           min_zodiac = 0.9, 
           # only the top structure tanimotoSimilarity >= 0.4 ...
           min_tanimoto = 0.4, 
           ...
           ){
    ## ------------------------------------------------------------------------------------------
    ## check dirs ---- spectra
    if(dirs == "all"){
      dirs <- list.files(path = path, pattern="^[0-9](.*)_(.*)_(.*)$", full.names = F)
      cat("## method_formula_based_spec_compare: check_dir\n")
      check <- pbapply::pbsapply(dirs, check_dir, file = "spectra") %>% unname
    }else{
      check <- pbapply::pbsapply(dirs, check_dir, file = "spectra") %>% unname
    }
    ## ------------------------------------------------------------------------------------------
    ## lock on file location
    meta_dir <- dirs[which(check == T)] %>%
      data.frame() %>%
      dplyr::rename(dir = ".") %>%
      dplyr::mutate(.id = sapply(dir, grep_id)) %>%
      merge(.MCn.formula_set, by = ".id", all.x = T) %>%
      merge(.MCn.structure_set[, c(".id", "tanimotoSimilarity")], by = ".id", all.x = T)
    ## some .id were Avoid time-consuming calculation
    if(nrow(meta_dir) >= filter_only_max){
      meta_dir <- dplyr::filter(meta_dir, ZodiacScore >= min_zodiac, tanimotoSimilarity >= min_tanimoto)
    }
    meta_dir <- dplyr::mutate(meta_dir, adduct_trans = gsub(" ", "", adduct),
                              target = paste0(precursorFormula, "_", adduct_trans, ".tsv"), 
                              full.name = paste0(path, "/", dir, "/", "spectra", "/", target), 
                              ## these files need to be check and filter (whether exist)
                              spectra = file.exists(full.name))
    meta_dir_filter <- dplyr::filter(meta_dir, spectra == T)
    cat("## STAT of spectra dataset:",
        paste0(nrow(meta_dir_filter), "(formula with match spectra)", "/", nrow(meta_dir), "(all formula)"), 
        "\n")
    ## ------------------------------------------------------------------------------------------
    ## load all spectra dataset
    spectra_cache <- new.env()
    pbapply::pbmapply(read_as_spectrum2, # function
                      meta_dir_filter$full.name,
                      meta_dir_filter$".id",
                      MoreArgs = list(
                                      cache = spectra_cache
                                      ))
    ## ------------------------------------------------------------------------------------------
    ## enumeration combination
    if(only_identical_class == T){
      ## enumeration combination in each hierarchy
      combn <- dplyr::filter(.MCn.nebula_index,
                             hierarchy >= min_hierarchy,
                             .id %in% meta_dir_filter$".id") %>%
        dplyr::group_by(hierarchy) %>%
        ## dispose in each group
        dplyr::summarise_at(c(".id"), unique) %>%
        dplyr::summarise_at(c(".id"), sort) %>%
        ## enumerate possible
        dplyr::summarise_at(c(".id"), combn_edge) %>%
        dplyr::ungroup() %>%
        dplyr::select(.id) %>%
        dplyr::distinct()
      ## this column has two sub-colunm
      combn <- dplyr::as_tibble(combn$".id")
    ## ------------------------------------------------------------------------------------------
    ## this cost too much time !!!
    }else{
      combn <- combn_edge(meta_dir_filter$".id")
    }
    ## ------------------------------------------------------------------------------------------
    ## compareSpectra (ms2) (via MSnbase)
    cat("## Method part: compare_spectra: sum:", nrow(combn), "\n")
    combn[[compare_fun]] <- pbapply::pbapply(combn, 1, couple_ms2_compare,
                                             fun = compare_fun,
                                             cl = cpu_cores,
                                             cache = spectra_cache,
                                             ...)
    ## filter via spectra similarity (edge_filter)
    ## ------------------------------------------------------------------------------------------
    combn <- combn[which(combn[[compare_fun]] >= edge_filter), ]
    if(precursor_mass_diff == T){
      ## dir
      info_path = paste0(path, "/", meta_dir_filter$dir, "/", "compound.info")
      ## load file
      cat("## Method part: load compound info file\n")
      meta_dir_filter$compound_mass <- pbapply::pblapply(info_path, get_precursor_mass) %>%
        unlist()
      ## compute mass difference
      cat("## Method part: diff_precursor_mass: sum:", nrow(combn), "\n")
      combn[["mass_diff"]] <- pbapply::pbapply(combn[,1:2], 1, precursor_mass_diff,
                                               df = data.table::data.table(meta_dir_filter))
    }
    combn <- dplyr::as_tibble(combn)
    return(combn)
    ## ------------------------------------------------------------------------------------------
  }
read_as_spectrum2 <-
  function(filename,
           key_id,
           cache = spectra_cache){
    file <- read_tsv(filename)
    file <- new("Spectrum2", mz = file$mz, intensity = file$rel.intensity)
    assign(paste0(key_id), file, envir = cache)
    return()
  }
combn_edge <-
  function(x){
    combn <- combn(x, 2)
      combn <- t(combn)
      combn <- data.frame(combn)
      colnames(combn) <- c(".id_1", ".id_2")
      return(combn)
  }
couple_ms2_compare <-
  function(x,
           fun = "dotproduct",
           cache = spectra_cache,
           ...){
    simi <- MSnbase::compareSpectra(get(x[1], envir = cache),
                                    get(x[2], envir = cache),
                                    fun = fun,
                                    ...)
    return(simi)
  }
get_precursor_mass <-
  function(
           path
           ){
    df <- read_tsv(path)
    mass <- df[index == "ionMass", 2]
    mass <- as.numeric(mass)
    return(mass)
  }
precursor_mass_diff <-
  function(
           x,
           df
           ){
    ##
    x1 = df[.id == x[1], "compound_mass", with = F]
    x2 = df[.id == x[2], "compound_mass", with = F]
    x = x2 - x1
    return(x)
  }
