#' @description Plotting correlation matrices in ggplot
#' 
#'
#' @param corr.matrix A NxN data.frame or matrix with the correlation scores. Alternatively a Nx(N+1), where the extra column indicates the sample names. The function can guess automatically or the column number can be indicated in the 'row.names.column' parameter.
#' @param row.names.column Numeric value indicating the column containing the sample names (row names). Default \code{NULL}.
#' @param clustering.method Method passed to hclust for the clustering. Default \code{NULL}. Options (an unambiguous abbreviation of): one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param dendrogram.position String indicating where should be positioned the dendrogram. Default \code{"right"}. Options: "right" or "left".
#' @param palette Color vector indicating the colors to use for the correlation scores (cells filling). Default: \code{rev(viridis::viridis(100))}
#' @param dendrogram.color String indicating the color to use for the dendrogram. Default: \code{"black"}.
#' @param dendrogram.linewidth Numeric value indicating the thickness of the dendrogram lines. Default: \code{0.5}.
#' @param correlation.limits 2-length numeric vector indication minimum and maximum of the correlation limits. Default: \code{c(0,1)}; examples: c(NA, 0.5), c(0, NA), c(0.25, 0.75).
#' @param title String indicating the plot title.
#' @param add.numbers Logic values to indicate whether correlation values should be displayed on the plot cells.
#' @param numbers.digits Number of digits to use for the correlation values.
#' @param numbers.color String indicating the color to use for the correlation values. Default: \code{"white"}.
#' @param numbers.font.size Number indicating the font size of the correlation values. Default: \code{"NA"} (automatic).
#'
#' @retunr A list containing:
#'  \itemize{
#'   \item{corr.matrix}{The input correlation matrix;}
#'   \item{cluster.tree}{The hclust object used for the clustering and dendrogram;}
#'    \item{heatmap}{A ggplot object with the heatmap}
#'  }
#' 
#' 
#' @export plot.correlation



plot.correlation =
  function(corr.matrix,
           row.names.column = NULL,
           clustering.method = "complete",
           dendrogram.position = "right",
           palette = rev(viridis::viridis(100)),
           dendrogram.color = "black",
           dendrogram.linewidth = 0.5,
           correlation.limits = c(0,1),
           title = NULL,
           add.numbers = FALSE,
           numbers.digits = 2,
           numbers.color = "white",
           numbers.font.size = NA) {
    
    ## Required packages
    #for (i in c("data.table", "ggh4x", "ggplot2", "reshape2", "dplyr", "viridis")) {if (!require(i)) install.packages(i)}
    require(dplyr)
    
    
    
    ## Read correlation matrix
    if ("character" %in% class(corr.matrix)) {
      matrix = data.frame(data.table::fread(paste0(corr.matrix), stringsAsFactors = F),
                          check.names = F, stringsAsFactors = F)
    } else {
      matrix = as.data.frame(corr.matrix)
    }
    colnames(matrix) = gsub("'","",colnames(matrix))
    
    
    if (!is.null(row.names.column)) {
      matrix[,row.names.column] = gsub("'","",matrix[,row.names.column])
      colnames(matrix)[row.names.column] = "sample"
    } else if (nrow(matrix) == ncol(matrix)) {
      matrix$sample = colnames(matrix)
    } else {
      # guess sample column
      col.class = sapply(1:ncol(matrix), function(x){class(matrix[,x])})
      if ("character" %in% col.class) {
        if (length(col.class[col.class == "character"]) == 1) {
          sample.column = which(col.class == "character")
          matrix[,sample.column] = gsub("'","",matrix[,sample.column])
          colnames(matrix)[sample.column] = "sample"
        }
      } else {
        warning("You should provide a data.frame or matrix with equal n.col and n.row in the same order.\n  Alternatively it is possibile to indicate a number corrresponding to the row.name labels by the option: 'row.names.column'.")
        invisible(return(invisible()))
      }
    }
    sample.col.position = grep("sample", colnames(matrix))
    
    
    ### Compute clustering
    corr_clust = hclust(d = as.dist(1-matrix[,-sample.col.position]), method = clustering.method)
    corr_clust$dist.method = "1 - correlation"
    
    
    ### Generate basic plot
    corr_heatmap =
      ggplot(data =
               reshape2::melt(data = matrix,
                              value.name = "correlation",
                              id.vars = "sample") %>%
               dplyr::mutate(sample = factor(sample, levels = rev(corr_clust$labels[corr_clust$order]))) %>%
               dplyr::filter(sample != variable),
             aes(x = sample,
                 y = variable,
                 fill = correlation)) + 
      geom_tile(title) +
      xlab(NULL) +
      ylab(NULL) +
      ggh4x::scale_y_dendrogram(hclust = corr_clust,
                                expand = c(0,0),
                                position = dendrogram.position) +
      scale_x_discrete(expand = c(0,0)) +
      scale_fill_gradientn(colours = palette,
                           limits = correlation.limits,
                           na.value = palette[ncol(palette)]) +
      coord_fixed() +
      ggtitle(title) +
      theme(axis.text = element_text(color = "black"),
            axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.ticks = element_line(color = dendrogram.color, linewidth = dendrogram.linewidth),
            axis.ticks.x = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", linewidth = 1, fill = NA))
    
    
    
    if (add.numbers == T) {
      corr_heatmap =
        corr_heatmap +
        geom_text(aes(label = round(correlation, digits = numbers.digits)),
                  color = numbers.color,
                  size = numbers.font.size)
    }
    
    
    
    # Export results
    return(list(corr.matrix = corr.matrix,
                cluster.tree = corr_clust,
                heatmap = corr_heatmap))
  }