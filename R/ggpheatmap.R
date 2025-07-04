# pheatmap(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                         "RdYlBu")))(100),
#          kmeans_k = NA, breaks = NA, border_color = "grey60",
#          cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE,
#          cluster_cols = TRUE, clustering_distance_rows = "euclidean",
#          clustering_distance_cols = "euclidean", clustering_method = "complete",
#          clustering_callback = identity2, cutree_rows = NA, cutree_cols = NA,
#          treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows,
#                                  50, 0), treeheight_col = ifelse((class(cluster_cols) == "hclust") ||
#                                                                    cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA,
#          legend_labels = NA, annotation_row = NA, annotation_col = NA,
#          annotation = NA, annotation_colors = NA, annotation_legend = TRUE,
#          annotation_names_row = TRUE, annotation_names_col = TRUE,
#          drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA,
#          fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize,
#          angle_col = c("270", "0", "45", "90", "315"), display_numbers = F,
#          number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8
#          * fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL,
#          labels_col = NULL, filename = NA, width = NA, height = NA,
#          silent = FALSE, na_col = "#DDDDDD", ...
         
         

         
ggpheatmap = 
  function(mat,
           names.column = NULL,
           color.map = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(101),
           scale = "none",
           cluster.columns = TRUE,
           cluster.rows = TRUE,
           clustering.distance.rows = "euclidean",
           clustering.distance.columns = "euclidean",
           clustering.method = "complete",
           annotation.row = NULL,
           annotation.row.position = "right",
           annotation.column = NULL,
           annotation.column.position = "right",
           annotation.colors = NA,
           annotation.na.color = "gray",
           annotation.spacing = 0.25,
           annotation.thickness = annotation.spacing*3,
           annotation.border.color = "white",
           row.dendrogram.position = "right",
           row.dendrogram.color = "black",
           row.dendrogram.linewidth = 0.5,
           column.dendrogram.position = "top",
           column.dendrogram.color = "black",
           column.dendrogram.linewidth = 0.5,
           show.rownames = T,
           show.colnames = T,
           column.label.rotation.angle = 90,
           x.label = NULL,
           y.label = NULL,
           na.color = "gray",
           cell.border.size = 0.5,
           cell.border.color = NA,
           font.size = 10) {
    
    # libs
    require(ggplot2)
    require(ggh4x)
    # require(ggdendro)
    # require(ggnewscale)
    
    
    # check matrix format
    names.column = names.column[1]
    if (is.null(names.column)) {
      matrix = as.matrix(mat)
      if (!is.numeric(matrix)) {
        return(warning("Please provide a numeric table or indicate the name/number of the column corresponding to the rownames using the flag 'names.column = ...'."))
      }
    } else {
      if (is.numeric(names.column)) {
        if (names.column <= ncol(mat)) {
          matrix = as.matrix(mat[,-names.column])
          rownames(matrix) = as.character(mat[,names.column])
          if (!is.numeric(matrix)) {
            return(warning("Your data.frame contains columns that are not numeric (besides the names column).\nPlease provide a numeric table or indicate the name/number of the column corresponding to the rownames using the flag 'names.column = ...'."))
          }
        } else {
          return(warning(paste0("The `names.column` provided (", names.column,") is a number higher that the total number of columns of the mat.")))
        }
      } else if (is.character(names.column)) {
        if (names.column %in% colnames(mat)) {
          matrix = as.matrix(mat[,which(colnames(mat) != names.column)])
          rownames(matrix) = as.character(mat[,which(colnames(mat) == names.column)])
          if (!is.numeric(matrix)) {
            return(warning("Your data.frame contains columns that are not numeric (besides the names column).\nPlease provide a numeric table or indicate the name/number of the column corresponding to the rownames using the flag 'names.column = ...'."))
          }
        } else {
          return(warning(paste0("The `names.column` provided (", names.column,") is not among the mat colnames.")))
        }
      }
    }
    
    
    ## Scale matrix if required
    if (is.null(scale)) {scale = "none"}
    if (is.na(scale)) {scale = "none"}
    if (tolower(scale) != "none") {
      if (tolower(scale) %in% c("row", "rows")) {
        matrix = pheatmap:::scale_mat(mat = matrix, scale = "row")
      } else if (tolower(scale) %in% c("column", "columns")){
        matrix = pheatmap:::scale_mat(mat = matrix, scale = "column")
      } else {
        return(warning("THe parameter `scale` must be a value among: 'none', 'row'/'rows', 'column'/'columns'."))
      }
    }
    
    matrix[is.nan(matrix)] = NA
    
    
    
    # Row/column clustering (if required)
    if (cluster.columns == TRUE) {
      col_clust = hclust(dist(t(matrix), method = clustering.distance.columns), method = clustering.method)
    }
    
    if (cluster.rows == TRUE) {
      row_clust = hclust(dist(matrix, method = clustering.distance.rows), method = clustering.method)
    }
    
    
    # melt the matrix
    melted_mat = reshape2::melt(data = matrix,
                                value.name = "score",
                                varnames = c("y", "x"))
    
    
    # Generate the basic plot
    heatmap =
      ggplot(data = melted_mat,
            mapping = aes(x = x,
                          y = y,
                          fill = score)) +
      geom_tile(size = cell.border.size,
                color = cell.border.color) +
      xlab(x.label) +
      ylab(y.label) +
      scale_x_discrete(expand = c(0,0)) +
      scale_y_discrete(expand = c(0,0)) +
      theme_classic() +
      theme(axis.text = element_text(color = "black", size = font.size),
            axis.text.x = element_text(angle = column.label.rotation.angle, vjust = 0.5, hjust = 1),
            axis.ticks = element_blank(),
            axis.line = element_blank())
    
    
    # get gradient color limits
    if (scale != "none") {
      gradient_max = ceiling(max(abs(matrix), na.rm = T)*10)/10
      heatmap =
        heatmap +
        scale_fill_gradientn(colors = color.map,
                             breaks = c(-1,0,1)*gradient_max,
                             limits = c(-1,1)*gradient_max,
                             na.value = na.color)
    } else {
      heatmap = heatmap + scale_fill_gradientn(colors = color.map)
    }
    
    
    # Remove row/col names
    if (show.colnames) {
      heatmap = heatmap + theme(axis.text.x = element_text(color = "black", size = font.size))
    } else {
      heatmap = heatmap + theme(axis.text.x = element_blank())
    }
    
    if (show.rownames) {
      heatmap = heatmap + theme(axis.text.y = element_text(color = "black", size = font.size))
    } else {
      heatmap = heatmap + theme(axis.text.y = element_blank())
    }
    
    
    # Add dendrogram
    if (cluster.columns) {
      heatmap =
        heatmap +
        scale_x_dendrogram(hclust = col_clust,
                           expand = c(0,0),
                           position = column.dendrogram.position) +
        theme(axis.ticks.x = element_line(color = column.dendrogram.color,
                                          linewidth = column.dendrogram.linewidth,
                                          lineend = "round"))}
    
    if (cluster.rows) {
      heatmap =
        heatmap +
        scale_y_dendrogram(hclust = row_clust,
                           expand = c(0,0),
                           position = row.dendrogram.position) +
        theme(axis.ticks.y = element_line(color = row.dendrogram.color,
                                          linewidth = row.dendrogram.linewidth,
                                          lineend = "round"))}
    
    
    
    
    ### Add row annotations --------------------------------------------------------------
    if (!is.null(annotation.row) | !is.null(annotation.column)) {
      heatmap_build = ggplot_build(heatmap)
      row_positions = data.frame(label = row_clust$labels,
                                 ymin = row_clust$order - 0.5,
                                 ymax = row_clust$order + 0.5)
    }
    
    if (!is.null(annotation.row)) {
      # get location and extremity
      row.anno.start = ifelse(tolower(annotation.row.position) == "right",
                              yes = max(heatmap_build$layout$panel_params[[1]]$x.range),
                              no = min(heatmap_build$layout$panel_params[[1]]$x.range))
      
      row.anno.pos = row.anno.start
      
      for (i in 1:ncol(annotation.row)) {
        row.anno.pos = ifelse(tolower(annotation.row.position) == "right",
                              yes = row.anno.pos + annotation.spacing,
                              no = row.anno.pos - annotation.spacing)
        
        current_anno =
          dplyr::left_join(annotation.row %>%
                             dplyr::select(colnames(annotation.row)[i]) %>%
                             dplyr::mutate(label = rownames(annotation.row),
                                           xmin = row.anno.pos,
                                           xmax = ifelse(tolower(annotation.row.position) == "right",
                                                         yes = row.anno.pos + annotation.thickness,
                                                         no = row.anno.pos - annotation.thickness)),
                           row_positions,
                           by = "label")
        
        
        heatmap =
          heatmap +
          ggnewscale::new_scale_fill() +
          geom_rect(data = current_anno,
                    mapping = aes(xmin = xmin,
                                  xmax = xmax,
                                  ymin = ymin,
                                  ymax = ymax,
                                  fill = current_anno[,1]),
                    inherit.aes = F,
                    color = annotation.border.color,
                    show.legend = T) +
          guides(fill = guide_legend(title=colnames(annotation.row)[i]))
        
        # # modify scale color
        # if (colnames(annotation.row)[i] %in% names(annotation.colors)) {
        #   heatmap =
        #     heatmap +
        #     scale_fill_manual(values = annotation.colors[[which(names(annotation.colors) == colnames(annotation.row)[i])]],
        #                       guide = guide_legend(title = colnames(annotation.row)[i]),
        #                       na.value = annotation.na.color)
        # } else {
        #   heatmap =
        #     heatmap +
        #     scale_fill_manual(values = scales::hue_pal()(length(unique(current_anno[,1]))),
        #                       guide = guide_legend(title = colnames(annotation.row)[i]),
        #                       na.value = annotation.na.color)
        # }

        print(heatmap)
        
        
        # re-set start for a annotation new bar
        row.anno.pos = ifelse(tolower(annotation.row.position) == "right",
                              yes = unique(current_anno$xmax),
                              no = unique(current_anno$xmax))
      }
    }
    
    


  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ## Add annotaions
    #step5.draw your annotation
    annotation.rows.label.position <- ifelse(tree_position_cols=="bottom","right","left")
    annotation.col.label.position <- ifelse(tree_position_rows=="right","left","right")
    if(!is.null(annotation_rows)){
      annotation_rows <- rownames_to_column(annotation_rows,var = "none1")
      for(i in 1:(ncol(annotation_rows)-1)){
        annotation_rows$none <- factor(rep(names(annotation_rows)[i+1],nrow(annotation_rows)))
        rowlist <- list()
        rowanno <- ggplot()+
          geom_exec(geom_tile,data = annotation_rows,x="none1",y="none",fill=names(annotation_rows)[i+1])+
          scale_fill_manual(values =annotation_color[names(annotation_color)==names(annotation_rows)[i+1]][[1]])+
          theme(axis.title = element_blank(),axis.text.y = element_blank(),
                axis.ticks = element_blank(),panel.background = element_blank(),
                axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5))+
          theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
          coord_flip()+labs(fill=names(annotation_rows)[i+1])+
          scale_y_discrete(position = annotation.rows.label.position)
        if(annotation_position_rows=="left"){
          p <- p %> %insert_left(rowanno,width =annotation_width )
        }else{
          p <- p %>% insert_right(rowanno,width =annotation_width )
        }
      }
    }else{
      rowanno <- NULL
    }
    
    
    if(!is.null(annotation_cols)){
      annotation_cols <- rownames_to_column(annotation_cols,var = "none1")
      for(i in 1:(ncol(annotation_cols)-1)){
        annotation_cols$none <- factor(rep(names(annotation_cols)[i+1],nrow(annotation_cols)))
        collist <- list()
        colanno <- ggplot()+
          geom_exec(geom_tile,data = annotation_cols,x="none1",y="none",fill=names(annotation_cols)[i+1])+
          scale_fill_manual(values = annotation_color[names(annotation_color)==names(annotation_cols)[i+1]][[1]])+
          theme(axis.title = element_blank(),axis.text.x = element_blank(),
                axis.ticks = element_blank(),panel.background = element_blank())+
          theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
          labs(fill=names(annotation_cols)[i+1])+
          scale_y_discrete(position = annotation.col.label.position)
        if(annotation_position_cols=="top"){
          p <- p %>% insert_top(colanno,height =annotation_width)
        } else {
          p <- p %>% insert_bottom(colanno,height =annotation_width)
        }
      }
    } else {
      colanno <- NULL
    }
    
    
    
    
    
    
    
    
    
    
    
    return(heatmap)  
    
  } # END function

    
    
    
    
    
    
    
    
    