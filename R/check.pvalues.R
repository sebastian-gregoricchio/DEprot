#' @title check.pvalues
#'
#' @description Plots a volcano plot log2(FoldChange) x -log10(p.adjusted) of differential expression results
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast Number indicating the position of the contrast to use for the plotting. Default: \code{1}.
#' @param histogram.binwidth Numeric value indicating. Default: \code{0.05}
#' @param p.value.color String indicating the color to use for the p-values. Default: \code{"steelblue"}.
#' @param p.adjusted.color String indicating the color to use for the p-values adjusted. Default: \code{"indianred"}.
#'
#' @return A DEprot.pvalues object containing 3 slots: distribution histogram of p-values and p-values adjusted and, line plot of ranked p-values.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export check.pvalues


check.pvalues =
  function(DEprot.analyses.object,
           contrast = 1,
           histogram.binwidth = 0.05,
           p.value.color = "steelblue",
           p.adjusted.color = "indianred") {

    # ## Libraries
    # require(dplyr)
    # require(ggplot2)

    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      stop("The input must be an object of class 'DEprot.analyses'.")
      #return(invisible())
    }

    ### check and collect contrast
    if (is.numeric(contrast)) {
      if (contrast <= length(DEprot.analyses.object@analyses.result.list)) {
        data = DEprot.analyses.object@analyses.result.list[[contrast]]$results
        data =
          data %>%
          dplyr::arrange(p.value) %>%
          dplyr::mutate(rank = 1:nrow(data))
        # n.diff = DEprot.analyses.object@analyses.result.list[[contrast]]$n.diff
        contrasts.info = DEprot.analyses.object@contrasts[[contrast]]
      } else {
        stop("The 'contrast' indicated is not available.")
        #return(invisible())
      }
    } else {
      stop("The 'contrast' must be a numeric value.")
      #return(invisible())
    }


    ### get p-adjusted threshold
    padj_th = DEprot.analyses.object@differential.analyses.params$padj.th



    ### Plot p-value histogram
    histogram_pval =
      ggplot(data = data) +
      geom_histogram(mapping = aes(x = p.value),
                     binwidth = histogram.binwidth,
                     color = p.value.color,
                     fill = p.value.color,
                     alpha = 0.5) +
      geom_hline(yintercept = median(hist(data$p.value, breaks = 1/histogram.binwidth, plot = FALSE)$counts, na.rm = TRUE), linetype = 2, color = "gray30") +
      geom_vline(xintercept = padj_th, linetype = 3, color = "black") +
      scale_x_continuous(expand = c(0.01,0)) +
      scale_y_continuous(expand = c(0,0)) +
      xlab("*P*-value") +
      ylab("Count") +
      ggtitle(label = paste0("**",contrasts.info$var.1, "** *vs* **", contrasts.info$var.2, "**"),
              subtitle = "*P*-value distribution") +
      theme_classic() +
      theme(axis.text = element_text(color = "black"),
            axis.title = ggtext::element_markdown(),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5))




    ### Plot p-value ADJUSTED histogram
    histogram_padj =
      ggplot(data = data) +
      geom_histogram(mapping = aes(x = padj),
                     binwidth = histogram.binwidth,
                     color = p.adjusted.color,
                     fill = p.adjusted.color,
                     alpha = 0.5) +
      geom_hline(yintercept = median(hist(data$padj, breaks = 1/histogram.binwidth, plot = FALSE)$counts, na.rm = TRUE), linetype = 2, color = "gray30") +
      geom_vline(xintercept = padj_th, linetype = 3, color = "black") +
      scale_x_continuous(expand = c(0.01,0)) +
      scale_y_continuous(expand = c(0,0)) +
      xlab("*P*~adjusted~") +
      ylab("Count") +
      ggtitle(label = paste0("**",contrasts.info$var.1, "** *vs* **", contrasts.info$var.2, "**"),
              subtitle = paste0("*P*~adjusted~ distribution<br>(adjust method: ",
                                DEprot.analyses.object@differential.analyses.params$padj.method,")")) +
      theme_classic() +
      theme(axis.text = element_text(color = "black"),
            axis.title = ggtext::element_markdown(),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5))




    ### Plot ranked p-values
    rank_plot =
      ggplot(data = data,
             aes(x = rank)) +
      geom_line(mapping = aes(y = p.value),
                color = p.value.color) +
      geom_line(mapping = aes(y = padj),
                color = p.adjusted.color) +
      geom_hline(yintercept = padj_th, linetype = 2) +
      scale_x_continuous(expand = c(0.01,0)) +
      scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
      annotate(label = "\u2014 P-value", geom = "text", x = -Inf, y = +Inf, color = p.value.color, hjust = -0.2, vjust = 2.5) +
      annotate(label = "\u2014 P-adjusted", geom = "text", x = -Inf, y = +Inf, color = p.adjusted.color, hjust = -0.16, vjust = 5) +
      annotate(label = paste0("Threshold: ", padj_th), geom = "text", x = -Inf, y = padj_th, color = "black", hjust = -0.1, vjust = -0.5) +
      ggtitle("Sorted *P*-values") +
      ggtitle(label = paste0("**",contrasts.info$var.1, "** *vs* **", contrasts.info$var.2, "**"),
              subtitle = "Sorted *P*-values") +
      ylab("*P*-value") +
      xlab("Rank") +
      theme_classic() +
      theme(axis.text = element_text(color = "black"),
            axis.title = ggtext::element_markdown(),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5))


    ### Return pvalues object
    DEprot.pvalues.object =
      new(Class = "DEprot.pvalues",
          pvalue.distribution = histogram_pval,
          padjusted.distribution = histogram_padj,
          pvalue.rank = rank_plot)

    return(DEprot.pvalues.object)
  } # END function

