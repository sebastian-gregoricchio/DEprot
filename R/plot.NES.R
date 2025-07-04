#' @title GSEA analyses NES plotter
#'
#' @description Plotting of a bar plot indicating the Normalized Enrichment Score (NES) of a gsea result. Bars will have different colors by positive or negative enrichment and transparency proportional to their p-value.
#'
#' @param gsea.object A \code{gseaResult} object as generate by \href{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}{clusterprofiler}.
#' @param pos.NES.label String for the label of sets enriched in the group with positive ranking feature scores. Default: \code{"+NES"}.
#' @param neg.NES.label String for the label of sets enriched in the group with negative ranking feature scores. Default: \code{"-NES"}.
#' @param pos.NES.color String for the color of sets enriched in the group with positive ranking feature scores. Default: \code{"steelblue"}.
#' @param neg.NES.color String for the color of sets enriched in the group with negative ranking feature scores. Default: \code{"orange"}.
#' @param string.pattern.to.remove String with a regular expression of a pattern to be removed from the set names. Default: \code{"HALLMARK|GOBP"}.
#' @param alpha.range Numeric vector of length 2 indicating minimum and maximum value for the transparency. Individual values must be a number between 0 and 1. Default: \code{c(0.3, 1)}.
#' @param add.counts Logic value to indicate whether add labels with the set size counts. Default: \code{"TRUE"}.
#' @param perc.bleeding.x Numeric value indicating the percentage of the full x.axis range to add on the left on the right. Useful when labels are falling outside for the x-max. Default: \code{8} (%).
#' @param axes.text.size Numeric value indicating the font size of the axis text. Default: \code{10}.
#' @param title String indicating the title of the plot. Default: \code{"NES enrichments"}.
#'
#' @return A ggplot object.
#' 
#' @export plot.NES 


plot.NES = function(gsea.object,
                    pos.NES.label = "+NES",
                    neg.NES.label = "-NES",
                    pos.NES.color = "steelblue",
                    neg.NES.color = "orange",
                    string.pattern.to.remove = "HALLMARK|GOBP",
                    alpha.range = c(0.3,1),
                    add.counts = TRUE,
                    perc.bleeding.x = 8,
                    axes.text.size = 10,
                    title = "NES enrichments") {
  
  # libraries
  require(dplyr)
  require(ggplot2)
  
  
  # extract results and clean
  result =
    gsea.object@result %>%
    dplyr::arrange(desc(NES)) %>%
    dplyr::mutate(alias = gsub("_", " ", gsub(string.pattern.to.remove, "", ID)),
                  dataset = ifelse(NES >= 0, yes = pos.NES.label, no = neg.NES.label)) %>%
    dplyr::mutate(alias = factor(alias, levels = rev(alias)))
  
  geneSets_size = data.frame(sapply(gsea.object@geneSets, length), stringsAsFactors = F)
  geneSets_size$ID = rownames(geneSets_size)
  colnames(geneSets_size)[1] = "geneSet_size"
  
  result = dplyr::left_join(result, geneSets_size, by = "ID")
  
  
  # define colors
  NES.colors = c(pos.NES.color, neg.NES.color)
  names(NES.colors) = c(pos.NES.label, neg.NES.label)
  
  max_abs.NES = ceiling(max(abs(result$NES), na.rm = T))
  
  # Generating the plot
  NES.plot =
    ggplot(data = result,
           aes(x = NES,
               y = alias,
               fill = dataset)) +
    geom_bar(aes(alpha = -log10(p.adjust)),
             stat = "identity",
             show.legend = T,
             width = 0.8) +
    scale_alpha_continuous(range = alpha.range) +
    scale_fill_manual(values = NES.colors, drop = F) +
    ylab(NULL) +
    ggtitle(title) +
    theme_classic() +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(color = "black",
                                   size = axes.text.size),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  
  # Add counts
  if (add.counts == TRUE) {
    if (nrow(result %>% dplyr::filter(NES >= 0)) > 0) {
      NES.plot =
        NES.plot +
        geom_text(data = result %>% dplyr::filter(NES >= 0),
                  aes(x = NES,
                      y = alias,
                      label = paste0(setSize,"/",geneSet_size)),
                  vjust = 0.5,
                  hjust = -0.2)
    }
    
    if (nrow(result %>% dplyr::filter(NES < 0)) > 0) {
      NES.plot =
        NES.plot +
        geom_text(data = result %>% dplyr::filter(NES < 0),
                  aes(x = NES,
                      y = alias,
                      label = paste0(setSize,"/",geneSet_size)),
                  vjust = 0.5,
                  hjust = 1.2)
    }
  } else {
    perc.bleeding.x = 0
  }
  
  NES.plot =
    NES.plot +
    #xlim(c(-1,1)*(max_abs.NES + (perc.bleeding.x/100)*max_abs.NES)) +
    scale_x_continuous(limits = c(-1,1)*(max_abs.NES + (perc.bleeding.x/100)*max_abs.NES),
                       expand = c(0,0))
  
  
  # export data
  return(NES.plot)
  
} #END plot.NES function