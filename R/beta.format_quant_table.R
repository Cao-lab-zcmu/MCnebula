format_quant_table <- 
  function(
           file,
           meta.group = c(blank = "blank", raw = "raw", pro = "pro"),
           from = "mzmine",
           get_metadata = F
           ){
    df <- data.table::fread(file) %>% 
      dplyr::rename(.id = 1, mz = 2, rt = 3) %>% 
      dplyr::select(1, 2, 3, contains("Peak area"))
    ## ------------------------------------- 
    colnames(df) <- gsub("\\.mz.{0,1}ML Peak area", "", colnames(df))
    metadata <- meta.group %>% 
      lapply(function(vec){
               str <- .meta_find_and_sort(colnames(df), vec)
           })
    metadata <- mapply(metadata, names(metadata), SIMPLIFY = F,
                       FUN = function(vec, name){
                         df <- data.table::data.table(group = name, sample = vec)
                         return(df)
                       })
    metadata <- data.table::rbindlist(metadata)
    ## ------------------------------------- 
    if(get_metadata)
      return(metadata)
    ## ------------------------------------- 
    stat <- df %>%
      dplyr::select(-mz, -rt) %>% 
      ## as long table
      reshape2::melt(id.vars = ".id", variable.name = "sample", value.name = "value") %>%
      merge(metadata, by = "sample", all.y = T) %>% 
      ## as data.table
      data.table::data.table() %>% 
      dplyr::mutate(value = as.numeric(value)) %>% 
      ## calculate average
      .[, list(mean = mean(value, na.rm = T)), by = list(.id, group)] %>%
      ## NAN as 0
      dplyr::mutate(mean = ifelse(is.nan(mean), 0, mean)) %>%
      ## as wide data
      data.table::dcast(.id ~ group, value.var = "mean") %>%
      ## .id is character
      dplyr::mutate(.id = as.character(.id)) %>% 
      dplyr::as_tibble() 
    return(stat)
  }
## ---------------------------------------------------------------------- 
.meta_find_and_sort <-
  function(
           name_set,
           pattern_set
           ){
    name_set <- lapply(pattern_set, .meta_mutate_grep_get,
                       string_set = name_set) %>%
      unlist()
    return(name_set)
  }
.meta_mutate_grep_get <-
  function(
           pattern,
           string_set
           ){
    string <- string_set %>%
      .[grepl(pattern, .)]
    return(string)
  }
## ---------------------------------------------------------------------- 
