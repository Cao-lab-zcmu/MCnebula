vis_via_molconvert_nebulae <-
  function(
           nebula_name
           ){
    df <- dplyr::filter(.MCn.nebula_index, name == nebula_name)
    stru <- dplyr::filter(.MCn.structure_set, .id %in% df$".id")
    vis_via_molconvert(stru$smiles, stru$".id")
    return("Done")
  }
vis_via_molconvert <-
  function(
           smiles_set,
           id_set,
           output = paste0(.MCn.output, "/", .MCn.results, "/tmp/structure")
           ){
    if(file.exists(output) == F){
      dir.create(paste0(.MCn.output, "/", .MCn.results))
      dir.create(output)
    }
    pbapply::pbmapply(molconvert_structure,
           smiles_set,
           id_set,
           MoreArgs = list(output = output)
    )
    return("Done")
  }
## ---------------------------------------------------------------------- 
molconvert_structure <-
  function(
           smiles,
           id,
           output = paste0(.MCn.output, "/", .MCn.results, "/tmp/structure")
           ){
    file = tempfile()
    system(paste0("molconvert mol \"", smiles, "\" -o ", file))
    system(paste0("obabel ", file, " -imol -osvg -O ", file, " > /dev/null 2>&1"))
    system(paste0("sed -i \'s/white/transparent/g\' ", file))
    system(paste0("cairosvg -f svg ", file, " -o ", output, "/", id, ".svg"))
    return()
  }
