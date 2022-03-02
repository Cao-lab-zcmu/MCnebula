#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param metadata PARAM_DESCRIPTION
#' @param tmp_stru PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{filter}}
#'  \code{\link[pbapply]{pbapply}}
#' @rdname batch_mode_structure
#' @export 
#' @importFrom dplyr mutate filter
#' @importFrom pbapply pbmapply
batch_mode_structure <-
  function(
           metadata,
           tmp_stru
           ){
    ## ---------------------------------------------------------------------- 
    ## collate metadata
    meta_stru <- dplyr::mutate(metadata,
                               stru_file = paste0(tmp_stru, "/", .id, ".svg"),
                               stru_check = file.exists(stru_file))
    meta_stru <- merge(meta_stru, .MCn.structure_set[, c(".id", "smiles")], by = ".id", all.x = T)
    meta_stru <- dplyr::filter(meta_stru, is.na(smiles) == F)
    cat("## STAT of structure set:",
        paste0(nrow(meta_stru), "(compounds with structure)", "/", nrow(metadata), "(all compounds)"),
        "\n")
    meta_stru <- dplyr::filter(meta_stru, stru_check == F)
    ## ---------------------------------------------------------------------- 
    if(nrow(meta_stru) > 0){
      pbapply::pbmapply(
                        base_vis_structure, # function
                        meta_stru$".id", # key_id
                        meta_stru$smiles, # smiles
                        MoreArgs = list( path = normalizePath(tmp_stru) ))
    }
  }
base_vis_structure <-
  function(
           key_id,
           smiles,
           path,
           to_file = paste0(path, "/", key_id, ".svg")
           ){
    ## openbabal. only support for linux
    ChemmineOB::convertToImage("SMI", "SVG", source = smiles, toFile = to_file)
    ## as transparent bg
    svg <- data.table::fread(file = to_file, sep = "", quote ="", header = F)
    svg$V1 = sub('fill="white"', 'fill="transparent"', svg$V1)
    write.table(x = svg, file = to_file, sep = "", col.names = F, row.names = F, quote = F)
    ## convert as cairo svg
    rsvg::rsvg_svg(to_file, to_file)
    return("Done")
  }

