#' @title plot.correlation.heatmap
#'
#' @description Function to generate a correlation heatmap. Includes dendrogram and clustering of the data.
#'
#' @param DEprot.object An object of class \code{DEprot}.
#' @param palette Vector of colors corresponding to the palette to use for the heatmap color scale. Default: \code{viridis::mako(100, direction = -1)}.
#' @param correlation.method String indicating the clustering method to use to generate the correlation matrix. Possible options: 'pearson', 'spearman', 'kendall'. Default: \code{"pearson"}.
#' @param sample.subset Vector indicating the name of the columns (\code{column.id} in the metadata table) to use/subset for the correlation. Default: \code{NULL} (no subsetting).
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'imp'. Default: \code{"imputed"}.
#' @param correlation.scale.limits Two-elements vector to indicate lower and higher limits, respectively, to apply to the correlation coefficient color scale. Default: \code{c(0,1)}.
#' @param exclude.diagonal Logical value indicating whether the plot diagonal (y = x) should be omitted. Default: \code{FALSE}.
#' @param dendrogram.position String indicating the position of the dendrogram. One among: "top", "bottom", "left", "right". Default: \code{"left"}.
#' @param dendrogram.color String indicating the color of the dendrogram lines. Default: \code{"black"}.
#' @param dendrogram.linewidth Numeric value indicating the line.width of the dendrogram. Default: \code{"black"}.
#' @param display.values Logical value indicating whether the correlation coefficient should be displayed for each cell. Default: \code{TRUE}.
#' @param values.color String indicating the color to use for the correlation coefficient labels. Default: \code{"NULL"}: colors will be assigned automatically depending on the background contrast, white for dark colors and black for light ones.
#' @param values.decimals Numeric value indicating the number of decimals at which round the correlation coefficient labels. Default: \code{2}.
#' @param values.font.size Numeric value indicating the font size of the correlation coefficient labels. Default: \code{2}.
#' @param values.transparency Numeric value between 0-1 indicating the transparency (alpha) of the correlation coefficient labels. Default: \code{1}, full color.
#' @param plot.title String indicating the main title of the plot. Default: \code{paste(stringr::str_to_title(correlation.method), "correlation")}.
#' @param plot.subtitle String indicating the subtitle of the plot. Default: \code{NULL}.
#' @param clustering.method String indicating the clustering method to use. The value should be (an unambiguous abbreviation of) one among: 'ward.D', 'ward.D2', 'single', 'complete', 'average' (= UPGMA), 'mcquitty' (= WPGMA), 'median' (= WPGMC) or 'centroid' (= UPGMC).
#'
#' @return A \code{DEprot.correlation} with the correlation heatmap in ggplot format.
#'
#' @name plot.correlation.heatmap
#'
#' @import dplyr
#' @import ggplot2
#' @import legendry
#' @import ggdendro
#' @import viridis
#' @importFrom farver get_channel
#' @importFrom stats as.dist
#' @importFrom reshape2 melt
#' @import ggtext
#'
#' @author Sebastian Gregoricchio
#'
#' @export plot.correlation.heatmap



plot.correlation.heatmap =
  function(DEprot.object,
           correlation.method = "pearson", #c("pearson", "kendall", "spearman")
           sample.subset = NULL,
           which.data = "imputed",
           palette = viridis::mako(100, direction = -1),
           correlation.scale.limits = c(0,1),
           exclude.diagonal = FALSE,
           dendrogram.position = "left",
           dendrogram.color = "black",
           dendrogram.linewidth = 0.5,
           display.values = TRUE,
           values.color = NULL,
           values.decimals = 2,
           values.font.size = 2,
           values.transparency = 1,
           plot.title = paste0("**",stringr::str_to_title(correlation.method), " correlation**"),
           plot.subtitle = NULL,
           clustering.method = "complete") {

    # ### Libraries
    # require(dplyr)
    # require(ggplot2)
    # require(legendry)
    # # require(farver)

    ### Functions
    contrast = # to define the contrast to decide whether the numbers should be in white or black depending on the background
      function(colour) {
        out = rep("black", length(colour))
        light = farver::get_channel(colour, "l", space = "hcl")
        out[light < 50] = "white"
        out
      }

    autocontrast = aes(colour = after_scale(contrast(fill)))

    #############################################


    ### check object
    if (!("DEprot" %in% class(DEprot.object)) & !("DEprot.analyses" %in% class(DEprot.object))) {
      stop("The input must be an object of class 'DEprot'.")
      #return()
    }

    ### Check and extract table
    if (tolower(which.data) == "raw") {
      if (!is.null(DEprot.object@raw.counts)) {
        mat = DEprot.object@raw.counts
        data.used = "raw"
      } else {
        stop(paste0("Use of RAW counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        #return()
      }
    } else if (tolower(which.data) %in% c("norm", "normalized", "normal")) {
      if (!is.null(DEprot.object@norm.counts)) {
        mat = DEprot.object@norm.counts
        data.used = "normalized"
      } else {
        stop(paste0("Use of NORMALIZED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        #return()
      }
    } else if (tolower(which.data) %in% c("imputed", "imp", "impute")) {
      if (!is.null(DEprot.object@imputed.counts)) {
        mat = DEprot.object@imputed.counts
        data.used = "imputed"
      } else {
        stop(paste0("Use of IMPUTED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        #return()
      }
    } else {
      stop(paste0("The 'which.data' value is not recognized.\n",
                  "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
      #return()
    }


    ## Check parameters
    if (!(tolower(dendrogram.position) %in% c("top", "bottom", "left", "right", "rigth"))) {
      stop("The 'dendrogram.position' must be a value among: 'top', 'bottom', 'left', 'right'.")
    } else if (tolower(dendrogram.position) == "rigth") {
      dendrogram.position = "right"
    }


    if (!(clustering.method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))) {
      stop("The 'clustering.method' should be (an unambiguous abbreviation of) one of: 'ward.D', 'ward.D2', 'single', 'complete', 'average' (= UPGMA), 'mcquitty' (= WPGMA), 'median' (= WPGMC) or 'centroid' (= UPGMC).")
    }

    if (!(tolower(correlation.method) %in% c("pearson", "kendall", "spearman"))) {
      stop("The correlation method is not supported. Choose one among: 'pearson', 'kendall', 'spearman'.")
    }


    ### subset table
    if (!is.null(sample.subset)) {
      mat = mat[,which(colnames(mat) %in% sample.subset)]
      corr.meta = dplyr::filter(DEprot.object@metadata, column.id %in% sample.subset)
    } else {
      corr.meta = DEprot.object@metadata
    }


    ### Compute correlation
    corr.mat = cor(mat, method = tolower(correlation.method), use = "complete.obs")

    ## Define dendrogram
    distance = stats::as.dist(1-corr.mat)
    corr_clust = hclust(d = distance, method = tolower(clustering.method))


    ## Reshape table
    matrix =
      reshape2::melt(data = corr.mat,
                     value.name = "correlation") %>%
      dplyr::mutate(Var1 = factor(Var1, levels = rev(colnames(corr.mat)[corr_clust$order]))) %>%
      dplyr::mutate(Var2 = factor(Var2, levels = rev(colnames(corr.mat)[corr_clust$order])))

    if (exclude.diagonal == TRUE) {
      matrix = dplyr::filter(.data = matrix, Var1 != Var2)
    }


    ## Generate the basic plot
    corr_heatmap =
      ggplot(data = matrix,
             aes(x = Var1,
                 y = Var2,
                 fill = correlation)) +
      geom_tile() +
      xlab(NULL) +
      ylab(NULL) +
      coord_fixed() +
      ggtitle(label = plot.title, subtitle = plot.subtitle) +
      theme(axis.text = element_text(color = "black"),
            axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
            aspect.ratio = 1,
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5))


    if (tolower(dendrogram.position) %in% c("left", "right")) {
      corr_heatmap =
        corr_heatmap +
        scale_x_discrete(expand = c(0,0)) +
        legendry::scale_y_dendro(clust = corr_clust,
                                 expand = c(0,0),
                                 position = tolower(dendrogram.position)) +
        theme(axis.ticks.y = element_line(color = dendrogram.color,
                                          linewidth = dendrogram.linewidth),
              axis.ticks.x = element_blank())
    } else {
      corr_heatmap =
        corr_heatmap +
        scale_y_discrete(expand = c(0,0)) +
        legendry::scale_x_dendro(hclust = corr_clust,
                                 expand = c(0,0),
                                 position = tolower(dendrogram.position)) +
        theme(axis.ticks.x = element_line(color = dendrogram.color,
                                          linewidth = dendrogram.linewidth),
              axis.ticks.y = element_blank())
    }


    ## Add values if required
    if (display.values == TRUE) {
      if (!is.null(values.color)) {
        corr_heatmap =
          corr_heatmap +
          geom_text(aes(label = round(correlation, values.decimals)),
                    color = values.color,
                    size = values.font.size,
                    alpha = 1)
      } else {
        corr_heatmap =
          corr_heatmap +
          geom_text(aes(label = round(correlation, values.decimals), !!!autocontrast),
                    size = values.font.size,
                    alpha = 1)
      }
    }


    # Remake color map
    corr_heatmap =
      corr_heatmap +
      scale_fill_gradientn(name = "Correlation\ncoefficient",
                           colours = palette,
                           limits = correlation.scale.limits,
                           na.value = palette[ncol(palette)])


    # Build out object
    corr_clust$dist.method = "as.dist(1 - correlation.matrix)"
    corr_clust$call = paste0("hclust(d = as.dist(1 - correlation.matrix), method = ",clustering.method,")")


    DEprot.corr.object =
      new(Class = "DEprot.correlation",
          heatmap = corr_heatmap,
          sample.subset = sample.subset,
          data.used = which.data,
          corr.metadata = corr.meta,
          corr.matrix = corr.mat,
          distance = distance,
          cluster = corr_clust)

    return(DEprot.corr.object)
  } # END function


