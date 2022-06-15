## target compare specture in classes
method_target_spec_compare <- 
  function(
           nebula_name,
           nebula_index = .MCn.nebula_index,
           output = paste0(.MCn.output, "/", .MCn.results, "/", nebula_name, ".spec.tsv"),
           edge_filter = 0.5,
           ...
           ){
    idset <- dplyr::filter(nebula_index, name %in% nebula_name)$.id
    edges <- method_formula_based_spec_compare(
      target_ids = idset,
      only_identical_class = F,
      min_hierarchy = 1,
      ...
    )
    write_tsv(edges, output)
    return(edges)
  }
