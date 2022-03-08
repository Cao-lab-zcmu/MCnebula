compare_score <- 
  function(
           structure_df,
           formula_df,
           formula_zodiac_rank1,
           formula_info,
           fc
           ){
    if(is.na(fc) == F){
      ## get top score (of structure) formula
      top_struc_formula <- structure_df[1,]$molecularFormula
      ## -------------------------------------------------------------------
      ## ------------------- get score
      formula_structure_rank1 <-
        formula_df[ molecularFormula == top_struc_formula, formula_info, with = F]
      ## get score
      ## the data.frame is a data.table project
      score_rank1_zodiac <- as.numeric(formula_zodiac_rank1[1, ZodiacScore])
      score_rank1_structure <- as.numeric(formula_structure_rank1[1, ZodiacScore])
      ## -------------------------------------------------------------------
      ## -------------------------------------------------------------------
      ## ------------------- comparation
      if( score_rank1_zodiac >= (score_rank1_structure * fc) ){
        use_zodiac = T
        ## sometimes rank 1 zodiac formulae are plural, due to complex adduct type
        ## hence based on structure score to filter them
        structure_df <- structure_df[molecularFormula %in% formula_zodiac_rank1$molecularFormula, ]
      }else{
        use_zodiac = F
      }
    }else{
      use_zodiac = F
    }
    list <- list(structure_df, use_zodiac)
    return(list)
  }
