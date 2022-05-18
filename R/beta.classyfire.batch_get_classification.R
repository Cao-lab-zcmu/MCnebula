## ---------------------------------------------------------------------- 
## ---------------------------------------------------------------------- 
## ---------------------------------------------------------------------- 
## ---------------------------------------------------------------------- 
## classyfire curl classification
batch_get_classification <- 
  function(
           inchikey2d,
           dir_pubchem,
           dir_classyfire,
           ...
           ){
    rdata <- paste0(dir_classyfire, "/", "class.rdata")
    classes <- extract_rdata_list(rdata)
    ## ------------------------------------- 
    inchikey2d <- inchikey2d[!inchikey2d %in% names(classes)]
    ## ------------------------------------- 
    if(length(inchikey2d) == 0)
      return()
    ## ------------------------------------- 
    inchikey_set <- extract_rdata_list(paste0(dir_pubchem, "/", "inchikey.rdata"), inchikey2d)
    ## ------------------------------------- 
    list <- lapply(inchikey_set, function(df){
                     if("InChIKey" %in% colnames(df))
                       return(df)
                     return()
           })
    ## ------------------------------------- 
    ## get classyfire classification
    df <- data.table::rbindlist(list)
    df <- dplyr::mutate(df, inchikey2d = stringr::str_extract(InChIKey, "^[A-Z]{1,}"))
    ## use the function which writed based on classyfireR::get_classification
    auto_classyfire(df, dir_classyfire, ...)
    ## ------------------------------------- 
    ## gather classes
    packing_as_rdata_list(dir_classyfire, pattern = "^[A-Z]{14}$", rdata = "class.rdata", extra = classes)
  }
## ---------------------------------------------------------------------- 
auto_classyfire <-
  function(
           df,
           dir_classyfire,
           classyfire_cl = NULL,
           ...
           ){
    ## create access log
    log_file <- paste0(dir_classyfire, "/log")
    if(file.exists(log_file)){
      log_df <- data.table::fread(log_file)
      df <- dplyr::filter(df, !InChIKey %in% log_df$log)
      if(nrow(df) == 0)
        return()
    }
    ## ------------------------------------- 
    list <- by_group_as_list(df, "inchikey2d")
    ## this part can be multi-threads
    log <- pbapply::pblapply(list, base_auto_classyfire,
                             dir_classyfire = dir_classyfire,
                             cl = classyfire_cl) %>% 
      unlist(use.names = F) %>% 
      data.table::data.table(log = .)
    ## ------------------------------------- 
    if(exists("log_df"))
      log <- dplyr::bind_rows(log_df, log)
    write_tsv(log, file = log_file)
  }
base_auto_classyfire <-
  function(
           df,
           dir_classyfire
           ){
    inchikey2d <- df[1,][["inchikey2d"]]
    log <- lapply(df[["InChIKey"]], base2_classyfire,
                  inchikey2d = inchikey2d,
                  dir_classyfire = dir_classyfire)
    return(unlist(log, use.names = F))
  }
base2_classyfire <-
  function(
           inchikey,
           inchikey2d,
           dir_classyfire
           ){
    file = paste0(dir_classyfire, "/", inchikey2d)
    if(file.exists(file) == F){
      ## if not exists
      ch <- classyfireR::get_classification(inchikey)
    }else{
      return()
    }
    if(is.null(ch)){
      return(inchikey)
    }else{
      ch <- classyfireR::classification(ch)
      write_tsv(ch, file)
      return()
    }
  }
