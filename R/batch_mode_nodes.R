batch_mode_nodes <-
  function(
           metadata,
           tmp_ppcp,
           with_structure = 0,
           plot_ppcp = plot_ppcp,
           plot_ratio = F,
           ratio_df = NA,
           palette = .MCn.palette,
           palette_stat = .MCn.palette_stat,
           ...
           ){
    ## remove exist files
    lapply(list.files(tmp_ppcp, full.names = T), file.remove)
    ## ---------------------------------------------------------------------- 
    ## nodes color setting, which parallel to the nodes color in plot of visualize_child_nebula function
    meta_color <- dplyr::select(metadata, vis_class) %>%
      dplyr::distinct() %>%
      dplyr::arrange(vis_class)
    if(is.vector(attr(palette, "name"))){
      ## filter palette
      palette <- palette[which(names(palette) %in% meta_color$vis_class)]
      ## sort according to the order of 'vis_class'
      palette <- palette[order(names(palette), levels = meta_color$vis_class)]
      ## set color
      meta_color$nodes_color <- palette
    }else{
      meta_color$nodes_color <- palette[1:nrow(meta_color)]
    }
    ## gather color data
    meta_nodes <- merge(metadata, meta_color, by = "vis_class", all.x = T)
    ## ---------------------------------------------------------------------- 
    ## pick ppcp_dataset
    ppcp_dataset = .MCn.ppcp_dataset[which(names(.MCn.ppcp_dataset) %in% meta_nodes$".id")]
    ## sort data
    meta_nodes$".id" <- factor(meta_nodes$".id",
                              levels = names(ppcp_dataset))
    meta_nodes <- meta_nodes[order(meta_nodes$".id"), ]
    ## ---------------------------------------------------------------------- 
    ## ratio_df, extra peak area data
    if(plot_ratio){
      ratio_df <- dplyr::mutate(ratio_df, .id = as.character(.id))
      ratio_df <- merge(dplyr::select(meta_nodes, .id), ratio_df, all.x = T, by = ".id", sort = F)
      ## get list data
      ratio_df_list <- by_group_as_list(ratio_df, ".id")
    }else{
      ratio_df_list <- rep(0, nrow(meta_nodes))
    }
    ## ---------------------------------------------------------------------- 
    cat("## annotate_child_nebulae: batch_mode_nodes\n")
    pbapply::pbmapply(base_vis_nodes, # function
                      ppcp_dataset, # main 1
                      meta_nodes$nodes_color, # main 2
                      names(ppcp_dataset), # main 3, key_id
                      ratio_df_list, # main 4, draw pie diagram
                      MoreArgs = list(
                                      path = normalizePath(tmp_ppcp),
                                      with_structure = with_structure,
                                      plot_ppcp = plot_ppcp,
                                      plot_ratio = plot_ratio,
                                      palette_stat = palette_stat,
                                      ...))
  }
## function mannually draw nodes
base_vis_nodes <-
  function(
           ppcp, ## main 1
           nodes_color, ## main 2
           key_id = NA, ## main 3
           ratio_df = NA, ## main 4, draw pie diagram
           plot_ratio = F,
           plot_nodes_id = T,
           plot_ppcp = T,
           label_color = "black",
           with_structure = 0,
           path = ".",
           class_index = unique(.MCn.nebula_index$relativeIndex),
           palette_ppcp = colorRampPalette(.MCn.palette_ppcp)(length(class_index)),
           palette_stat = .MCn.palette_stat,
           size_adjust = 0.7
           ){
    ## ---------------------------------------------------------------------- 
    ## filter via class_index
    ppcp <- ppcp[ppcp$relativeIndex %in% class_index, ]
    ## as factor, for painting color
    ppcp$relativeIndex <- factor(ppcp$relativeIndex, levels = sort(ppcp$relativeIndex))
    ppcp$num <- seq(1, nrow(ppcp))
    if(plot_ppcp == F){
      ppcp$V1 = 0
    }
    ## ---------------------------------------------------------------------- 
    ## plot nodes 
    p <- ggplot(ppcp, aes(x = num, y = V1)) +
      ## nodes color
      geom_ribbon(fill = nodes_color,
                           aes(ymin = -5, ymax = 0,
                               x = ifelse(num == 1, 0,
                                          ifelse(num == nrow(ppcp), num + 1, num)))) +
      ## border color
      geom_ribbon(fill = "black",
                           aes(ymin = 0, ymax = 1.1,
                               x = ifelse(num == 1, 0,
                                          ifelse(num == nrow(ppcp), num + 1, num)))) +
        ## ppcp bar plot
        geom_col(alpha = 1, aes(fill = relativeIndex), color = "white", size = 0.25) +
        ## nodes border ratio
        ylim(-5, 1.3) +
        ## Polar coordinate transformation
        coord_polar() +
        theme_minimal() +
        theme(
              text = element_text(family="Times"),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              panel.grid = element_blank(),
              ## remove legend
              legend.position = "none",
              panel.border = element_blank(),
              plot.margin =unit(c(0,0,0,0),"cm"),
              panel.spacing =unit(c(0,0,0,0),"cm")
        )
    ## ---------------------------------------------------------------------- 
    ## draw pie diagram
    if(plot_ratio){
      ratio_df <- reshape2::melt(ratio_df, id.vars = ".id", variable.name = "group", value.name = "value")
      ## mutate NA as 0
      ratio_df <- dplyr::mutate(ratio_df, value = ifelse(is.na(value), 0, value))
      ## value stack
      ratio_df <- dplyr::mutate(ratio_df,
                                xend = stack_sum(ratio_df$value),
                                x = stack_sum(c(0, ratio_df$value[1:(nrow(ratio_df)-1)])))
      ## normalize x axis range and x value
      n_factor = (max(ppcp$num) + 1)/max(ratio_df$xend)
      ratio_df <- dplyr::mutate(ratio_df,
                         midd = (x + value/2) * n_factor,
                         width = value * n_factor)
      ## add pie plot into ggplot2 project
      p <- p + geom_tile(data = ratio_df, size = 0.2, color = "white",
                                  aes(y = -2.5, x = midd, width = width, height = 2.5, fill = group)) +
        ## add 'fill' palette
        scale_fill_manual(values = c(palette_ppcp, palette_stat[1:nrow(ratio_df)]))
    }else{
      ## add 'fill' palette
      p <- p + scale_fill_manual(values = palette_ppcp)
    }
    ## ---------------------------------------------------------------------- 
    ## generate Graphics Device
    savepath = paste0(path, "/", key_id, ".svg")
    svglite::svglite(savepath, bg = "transparent")
    ## print nodes 
    print(p)
    ## ---------------------------------------------------------------------- 
    ## print structure or not
    if(with_structure == 1){
      s_file = paste0(normalizePath(paste0(path, "/../structure")), "/", key_id, ".svg")
      if(file.exists(s_file)){
        ## via grImport2 import Cairo svg
        ps <- grImport2::readPicture(file = s_file)
        ## grid draw
        grImport2::grid.picture(ps, width = size_adjust, height = size_adjust)
      }
    }
    ## ---------------------------------------------------------------------- 
    ## grid nodes ID in nodes 
    if(plot_nodes_id == T){
      ## a grid object
      ps <- grid::textGrob(paste0("ID:", key_id),
                           y = 0.25,
                           gp = grid::gpar(fontfamily = "Times", fontsize = 20, col = label_color))
      grid::grid.draw(ps)
    }
    ## ---------------------------------------------------------------------- 
    dev.off()
    # as cairo svg
    rsvg::rsvg_svg(savepath, savepath)
    ## ---------------------------------------------------------------------- 
  }
## function read cairo svg
base_read_cairo <-
  function(
           key_id,
           path,
           suffix = ".svg"
           ){
    prefix <- c()
    svg <- grImport2::grobify(grImport2::readPicture(paste0(path, "/", key_id, suffix)))
    svg <- gridExtra::arrangeGrob(svg)
    return(svg)
  }
stack_sum <-
  function(
           vector
           ){
    stack <- c()
    for(i in 1:length(vector)){
      stack[i] <- sum(vector[1:i])
    }
    return(stack)
  }
