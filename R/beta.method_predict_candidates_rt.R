method_predict_candidates_rt <- 
  function(
           structure_set,
           reference_compound,
           rt_set,
           rt_weight,
           rt_window,
           ...
           ){
    reference_set <- merge(reference_compound, structure_set, by = c(".id", "structure_rank")) %>% 
      merge(rt_set[, c(".id", "RT")], by = ".id", all.x = T) %>% 
      dplyr::select(.id, inchikey2D, smiles, RT) %>% 
      dplyr::rename(NAME = .id, InChIKey = inchikey2D, SMILES = smiles) %>% 
      dplyr::mutate(InChIKey = paste0(InChIKey, "-UHFFFAOYSA-N"))
    ## ------------------------------------- 
    descs <- Retip::getCD(reference_set)
    return(descs)
    db_rt <- proc.data(descs)
  }
