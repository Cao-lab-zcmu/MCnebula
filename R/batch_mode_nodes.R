#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param metadata PARAM_DESCRIPTION
#' @param tmp_ppcp PARAM_DESCRIPTION
#' @param with_structure PARAM_DESCRIPTION, Default: 0
#' @param plot_ppcp PARAM_DESCRIPTION, Default: plot_ppcp
#' @param plot_ratio PARAM_DESCRIPTION, Default: F
#' @param ratio_df PARAM_DESCRIPTION, Default: NULL
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{distinct}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{mutate}}
#'  \code{\link[pbapply]{pbapply}}
#' @rdname batch_mode_nodes
#' @importFrom dplyr select distinct arrange mutate
#' @importFrom pbapply pbmapply
batch_mode_nodes <-
  function(
           metadata,
           tmp_ppcp,
           with_structure = 0,
           plot_ppcp = plot_ppcp,
           plot_ratio = F,
           ratio_df = NULL,
           ...
           ){
    ## remove exist files
    lapply(list.files(tmp_ppcp, full.names = T), file.remove)
    ## ---------------------------------------------------------------------- 
    ## nodes color, which parallel to the nodes color in plot of visualize_child_nebula function
    meta_color <- dplyr::select(metadata, vis_class) %>%
      dplyr::distinct() %>%
      dplyr::arrange(vis_class)
    meta_color$nodes_color <- .MCn.palette[1:nrow(meta_color)]
    ## gather color data
    meta_ppcp <- merge(metadata, meta_color, by = "vis_class", all.x = T)
    ## ---------------------------------------------------------------------- 
    ## pick ppcp_dataset
    ppcp_dataset = .MCn.ppcp_dataset[which(names(.MCn.ppcp_dataset) %in% meta_ppcp$".id")]
    ## sort data
    meta_ppcp$".id" <- factor(meta_ppcp$".id",
                              levels = names(ppcp_dataset))
    meta_ppcp <- meta_ppcp[order(meta_ppcp$".id"), ]
    ## ---------------------------------------------------------------------- 
    ## ratio_df, extra peak area data
    if(plot_ratio == T){
      ratio_df <- dplyr::mutate(ratio_df, .id = as.character(.id))
      ratio_df <- merge(dplyr::select(meta_ppcp, .id), ratio_df, all.x = T, by = ".id", sort = F)
      assign("envir_nebula", environment(), parent.env(environment()))
      ## get list data
      ratio_df_list <- lapply(ratio_df$".id", by_group_for_list,
                              df = get("ratio_df", envir = get("envir_nebula")),
                              col = ".id")
    }else{
      ratio_df_list <- rep(0, nrow(meta_ppcp))
    }
    ## ---------------------------------------------------------------------- 
    cat("## annotate_child_nebulae: batch_mode_nodes\n")
    pbapply::pbmapply(
                      base_vis_nodes, # function
                      ppcp_dataset, # main 1
                      meta_ppcp$nodes_color, # main 2
                      names(ppcp_dataset), # main 3, key_id
                      ratio_df_list, # main 4, draw pie diagram
                      MoreArgs = list(
                                      path = normalizePath(tmp_ppcp),
                                      with_structure = with_structure,
                                      plot_ppcp = plot_ppcp,
                                      plot_ratio = plot_ratio,
                                      ...))
  }
## function mannually draw nodes
base_vis_nodes <-
  function(
           ppcp, ## main 1
           nodes_color, ## main 2
           key_id = NULL, ## main 3
           ratio_df = NULL, ## main 4, draw pie diagram
           plot_ratio = F,
           plot_nodes_id = T,
           plot_ppcp = T,
           label_color = "black",
           with_structure = 0,
           path = ".",
           class_index = unique(.MCn.nebula_index$relativeIndex),
           palette = colorRampPalette(c(.MCn.palette))(length(class_index)),
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
      ggplot2::geom_ribbon(fill = nodes_color,
                           aes(ymin = -5, ymax = 0,
                               x = ifelse(num == 1, 0,
                                          ifelse(num == nrow(ppcp), num + 1, num)))) +
      ## border color
      ggplot2::geom_ribbon(fill = "black",
                           aes(ymin = 0, ymax = 1.1,
                               x = ifelse(num == 1, 0,
                                          ifelse(num == nrow(ppcp), num + 1, num)))) +
        ## ppcp bar plot
        ggplot2::geom_col(alpha = 1, aes(fill = relativeIndex), color = "white", size = 0.25) +
        ## nodes border ratio
        ggplot2::ylim(-5, 1.3) +
        ## Polar coordinate transformation
        ggplot2::coord_polar() +
        ggplot2::theme_minimal() +
        ggplot2::theme(
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
    if(plot_ratio == T){
      ratio_df <- reshape2::melt(ratio_df, id.vars = ".id", variable.name = "group", value.name = "value")
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
      p <- p + ggplot2::geom_tile(data = ratio_df, size = 0.2, color = "white",
                                  aes(y = -2.5, x = midd, width = width, height = 2.5, fill = group)) +
        ## add 'fill' palette
        ggplot2::scale_fill_manual(values = c(palette, palette_stat[1:nrow(ratio_df)]))
    }else{
      ## add 'fill' palette
      p <- p + ggplot2::scale_fill_manual(values = palette)
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
