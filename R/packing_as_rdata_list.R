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
           extra = NULL
           ){
    file_set <- list.files(path, pattern = pattern)
    if(length(file_set) == 0)
      return()
    ## read as list
    list <- pbapply::pblapply(paste0(path, "/", file_set), read_tsv)
    names(list) <- file_set
    ## merge
    list <- c(extra, list)
    ## save as rdata
    save(list, file = paste0(path, "/", rdata))
  }
