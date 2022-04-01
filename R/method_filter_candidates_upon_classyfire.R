method_filter_candidates_upon_classyfire <- 
  function(
           structure_set, 
           nebula_name,
           ...
           ){
    cat("## netbula_re_rank: method_filter_candidates_upon_classyfire\n")
    ## ---------------------------------------------------------------------- 
    tmp_dir <- paste0(.MCn.results, "/tmp")
    ## dir for pubchem data
    tmp_pub <- paste0(tmp_dir, "/pubchem")
    ## dir for classyfire data
    tmp_classyfire <- paste0(tmp_dir, "/classyfire")
    ## create dir
    lapply(c(tmp_dir, tmp_pub, tmp_classyfire), function(dir){
             if(file.exists(dir) == F)
               dir.create(dir)
           })
    ## ---------------------------------------------------------------------- 
    inchikey2d <- unique(structure_set$inchikey2D)
    pubchem_curl_inchikey(inchikey2d, tmp_pub, ...)
    ## ---------------------------------------------------------------------- 
    cat("## use classyfireR to get classification\n")
    auto_get_classification(inchikey2d, tmp_pub, tmp_classyfire, ...)
    ## ---------------------------------------------------------------------- 
  }
## ---------------------------------------------------------------------- 
## ---------------------------------------------------------------------- 
## ---------------------------------------------------------------------- 
## ---------------------------------------------------------------------- 
## pubchem curl inchikey
pubchem_curl_inchikey <- 
  function(
           inchikey2d,
           dir,
           curl_cl = NULL,
           gather_as_rdata = T,
           ...
           ){
    rdata <- paste0(dir, "/", "inchikey.rdata")
    inchikey_set <- extract_rdata_list(rdata)
    ## -------------------------------------
    inchikey2d <- inchikey2d[!inchikey2d %in% names(inchikey_set)]
    ## ------------------------------------- 
    if(length(inchikey2d) == 0)
      return()
    ## ------------------------------------- 
    cat("## curl pubchem to get possibly InChIKey\n")
    pbapply::pblapply(inchikey2d, base_pubchem_curl_inchikey,
                      dir = dir, cl = curl_cl, ...)
    ## ------------------------------------- 
    cat("## gather InChIKey\n")
    packing_as_rdata_list(dir, pattern = "^[A-Z]{14}$",
                          rdata = "inchikey.rdata", extra = inchikey_set)
  }
base_pubchem_curl_inchikey <- 
  function(
           inchikey2d,
           dir,
           type = "inchikey",
           get = "InChIkey",
           ...
           ){
    file <- paste0(dir, "/", inchikey2d)
    ## ------------------------------------- 
    ## if exists and valid, return
    if(file.exists(file)){
      csv <- read_tsv(file)
      if("CID" %in% colnames(csv))
        return()
    }
    ## ------------------------------------- 
    ## curl via inchikey2d, which the same as InChIKey 
    url_start = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/", type, "/")
    ## to get InChIKey, and get data as csv format
    url_end = paste0("/property/", paste(get, collapse = ","), "/CSV")
    ## gather url
    url = paste0(url_start, "/", inchikey2d, "/", url_end)
    ## ------------------------------------- 
    csv <- RCurl::getURL(url)
    ## ------------------------------------- 
    ## check results, if 404, return()
    if(grepl("Status: 404", csv)){
      write_tsv(csv, file = file)
      return()
    }
    ## ------------------
    ## if 503, system busy, try again
    while(grepl("Status:	503", csv)){
      csv <- RCurl::getURL(url)
    }
    ## ------------------------------------- 
    csv <- data.table::fread(text = csv)
    write_tsv(csv, file = file)
  }
## ---------------------------------------------------------------------- 
## ---------------------------------------------------------------------- 
## ---------------------------------------------------------------------- 
## ---------------------------------------------------------------------- 
## classyfire curl classification
auto_get_classification <- 
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
    df <- data.table::rbindlist(list)
    df <- dplyr::mutate(df, inchikey2d = stringr::str_extract(InChIKey, "^[A-Z]{1,}"))
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
    ## classyfireR
    list <- by_group_as_list(df, "inchikey2d")
    ## this part can be multi-threads
    pbapply::pblapply(list, base_auto_classyfire,
                      dir_classyfire = dir_classyfire,
                      cl = classyfire_cl)
  }
base_auto_classyfire <-
  function(
           df,
           dir_classyfire
           ){
    inchikey2d <- df[1,][["inchikey2d"]]
    lapply(df[["InChIKey"]], base2_classyfire,
           inchikey2d = inchikey2d,
           dir_classyfire = dir_classyfire)
  }
base2_classyfire <-
  function(
           inchikey,
           inchikey2d,
           dir_classyfire
           ){
    file = paste0(dir_classyfire, "/", inchikey2d)
    ch <- try(read_tsv(file), silent = T)
    if(class(ch)[1] == "try-error"){
      ## if not exists
      ch <- classyfireR::get_classification(inchikey)
    }else{
      return()
    }
    if(is.null(ch)){
      return()
    }else{
      ch <- classyfireR::classification(ch)
      write_tsv(ch, file)
    }
  }
