extract_rdata_list <- 
  function(
           rdata,
           names = NA
           ){
    if(file.exists(rdata) == F)
      return()
    load(rdata)
    if(is.na(names[1]) == F){
      list <- list[names(list) %in% names]
    }
    return(list)
  }
packing_as_rdata_list <- 
  function(
           path,
           pattern,
           rdata,
           extra = NULL,
           rm_files = T
           ){
    file_set <- list.files(path, pattern = pattern)
    if(length(file_set) == 0)
      return()
    ## read as list
    list <- pbapply::pblapply(paste0(path, "/", file_set), read_tsv)
    names(list) <- file_set
    ## merge
    list <- c(extra, list)
    ## according to name, unique
    df <- data.table::data.table(name = names(list), n = 1:length(list))
    df <- dplyr::distinct(df, name, .keep_all = T)
    list <- list[df$n]
    ## rm origin file sets
    if(rm_files == T){
      lapply(paste0(path, "/", file_set), file.remove)
    }
    ## save as rdata
    save(list, file = paste0(path, "/", rdata))
  }
