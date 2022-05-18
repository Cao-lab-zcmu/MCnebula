pubchem_get_synonyms <- 
  function(
           cid,
           dir,
           curl_cl = NULL,
           gather_as_rdata = T,
           group_number = 50,
           ...
           ){
    rdata <- paste0(dir, "/", "cid.rdata")
    ## extract as list
    cid_set <- extract_rdata_list(rdata)
    ## as data.table
    cid_set <- data.table::rbindlist(cid_set)
    ## !duplicated
    if("cid" %in% colnames(cid_set)){
      cid_set <- cid_set %>% 
        dplyr::distinct(cid, syno)
    }
    ## -------------------------------------
    ## exclude existing
    cid <- cid[!cid %in% cid_set$cid]
    ## ------------------------------------- 
    if(length(cid) == 0)
      return()
    ## ------------------------------------- 
    ## if grouped, the rest number
    rest <- length(cid) %% group_number
    ## assign group
    group <- matrix(cid[1:(length(cid) - rest)], ncol = group_number) %>% 
      ## use apply to multiple list
      apply(1, c, simplify = F) %>% 
      ## gather the rest cid
      c(., list(tail(cid, n = rest))) %>% 
      ## add group name
      mapply(FUN = function(vec, name){
               attr(vec, "name") <- name
               return(vec)
           }, ., paste0("G", 1:length(.)),
           SIMPLIFY = F)
    ## ------------------------------------- 
    if(rest == 0)
      group <- group[1:(length(group) - 1)]
    ## ------------------------------------- 
    pbapply::pblapply(group, base_pubchem_get_synonyms,
                      dir = dir, cl = curl_cl, ...)
    ## ------------------------------------- 
    cat("## gather synonyms\n")
    packing_as_rdata_list(dir, pattern = "^G[0-9]{1,}$",
                          dedup = F,
                          rdata = "cid.rdata",
                          ## data.table as list
                          extra = list(cid_set))
  }
base_pubchem_get_synonyms <- 
  function(
           cid,
           dir,
           ...
           ){
    savename <- attr(cid, "name")
    file <- paste0(dir, "/", savename)
    ## gather cid and sep by ,
    cid <- paste(cid, collapse = ",")
    ## use cid
    url_start <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    ## get as XML
    url_end <- "/synonyms/XML"
    ## paste as url
    url <- paste0(url_start, cid, url_end)
    ## ------------------------------------- 
    check <- 0
    while(check == 0 | class(check)[1] == "try-error"){
      check <- try(text <- RCurl::getURL(url), silent = T)
    }
    ## ------------------------------------- 
    while(grepl("Status:	503", text)){
      text <- RCurl::getURL(url)
    }
    ## "PUGREST.BadRequest"
    ## ------------------------------------- 
    ## convert to list
    text <- XML::xmlToList(text)
    ## only 'information'
    text <- text[names(text) == "Information"]
    ## in list to separate data
    text <- lapply(text, function(list){
                     syno <- list[names(list) == "Synonym"]
                     syno <- lapply(syno,
                                    function(char){
                                      if(is.null(char)){
                                        return(NA)
                                      }else{
                                        return(char)
                                      }
                                    })
                     data.table::data.table(cid = list$CID, syno = unlist(syno))
           })
    text <- data.table::rbindlist(text, fill = T)
    ## ------------------------------------- 
    ## save data
    write_tsv(text, filename = file)
  }

