by_group_as_list <-
  function(
           df,
           colnames
           ){
    vector <- unique(df[[colnames]])
    list <- lapply(vector, by_group_as_list_select,
                   df = df,
                   colNames = colnames)
    names(list) <- vector
    return(list)
  }
by_group_as_list_select <- 
  function(
           KEY,
           df,
           colNames
           ){
    df <- df[which(df[[colNames]] == KEY), ]
    return(df)
  }
