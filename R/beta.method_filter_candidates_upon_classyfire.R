method_filter_candidates_upon_classyfire <- 
  function(
           structure_set, 
           nebula_name,
           ...
           ){
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
    ## ------------------------------------- 
    cat("## Method part: base_pubchem_curl_inchikey\n")
    pubchem_curl_inchikey(inchikey2d, tmp_pub, ...)
    ## ------------------------------------- 
    cat("## Method part: batch_get_classification (classyfireR::get_classification)\n")
    batch_get_classification(inchikey2d, tmp_pub, tmp_classyfire, ...)
    ## ---------------------------------------------------------------------- 
    structure_set <- calculate_class_score(structure_set, nebula_name, tmp_classyfire, ...)
    return(structure_set)
  }
## ---------------------------------------------------------------------- 
calculate_class_score <- 
  function(
           structure_set,
           nebula_name,
           dir,
           ...
           ){
    ## get parent class of the nebula
    class_hierarchy <- mutate_get_parent_class(nebula_name, class_cutoff = 1) %>% 
      unlist(use.names = F) %>% 
      c(nebula_name, .) %>% 
      data.table::data.table(Classification = ., hierarchy = length(.):1) %>% 
      dplyr::filter(hierarchy >= 4)
    ## ------------------------------------- 
    inchikey2d <- unique(structure_set$inchikey2D)
    ## ------------------ 
    ## load existed class data
    class_set <- extract_rdata_list(paste0(dir, "/class.rdata"), inchikey2d) %>%
      data.table::rbindlist(idcol = T) %>% 
      dplyr::rename(inchikey2d = .id) %>% 
      ## merge with class hierarchy information
      merge(class_hierarchy, ., by = "Classification") %>% 
      dplyr::mutate(class_score = hierarchy * 2) %>% 
      dplyr::group_by(inchikey2d)
    ## calculate score
    class_score <- dplyr::summarise_at(class_set, "class_score", sum)
    ## merge with origin structure_sete
    structure_set <- merge(structure_set, class_score,
                           by.x = "inchikey2D", by.y = "inchikey2d", all.x = T) %>% 
      dplyr::mutate(class_score = ifelse(is.na(class_score), 0, class_score)) %>% 
      ## re-rank according to class_score
      dplyr::arrange(.id, desc(class_score), desc(score)) %>% 
      dplyr::distinct(.id, .keep_all = T) %>% 
      dplyr::relocate(.id) %>% 
      dplyr::as_tibble()
  }

