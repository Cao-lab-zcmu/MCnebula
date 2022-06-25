nebula_get_candidate <- 
  function(
           ...,
           path = "mcnebula_results"
           ){
    path <- paste0(path, "/candidates")
    if(file.exists(path) == F)
      dir.create(path)
    args <- list(
                 ...,
                 top_n = 50,
                 match_pattern = NULL, ## or c("precursorFormula", "adduct") or NULL
                 collate_factor = NA,
                 revise_MCn_formula_set = F,
                 revise_MCn_structure_set = F,
                 only_gather_structure = T
    )
    candidates <- do.call(nebula_re_rank, args)
    write_tsv(candidates, paste0(path, "/", args[[1]], ".tsv"))
  }
