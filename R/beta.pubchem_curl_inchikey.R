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

