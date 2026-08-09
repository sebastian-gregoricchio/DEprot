#' @title divergent.enrichment
#'
#' @description Combines two enrichment analyses in a single divergent (back-to-back) bar plot: the genesets enriched in the first element of the list are drawn as positive values, the ones enriched in the second element as negative values. It is meant to display together the two sides of a differential expression analyses, for instance the ORA of the up- and down-regulated proteins of the same contrast. The two enrichments can come from different sources (DEprot or clusterProfiler) and can be mixed.
#'
#' @param enrichment.list Named list of two enrichments. Each element can be an object of class \code{DEprot.enrichResult} (DEprot), \code{enrichResult} or \code{gseaResult} (clusterProfiler), or a data.frame with the same structure of the tables returned by clusterProfiler. The first element is plotted on the positive side, the second one on the negative side. No default.
#' @param value String indicating the metric used for the length of the bars. One among: \code{"FoldEnrichment"}, \code{"GeneRatio"}, \code{"Count"}, \code{"NES"}, \code{"padj"} (-log10 of the adjusted p-value). Default: \code{"FoldEnrichment"}.
#' @param top.n Numeric value indicating the maximum number of genesets to display for each side. Default: \code{10}.
#' @param terms Character vector indicating specific geneset IDs to display. When provided, \code{top.n} is ignored. Default: \code{NULL}.
#' @param padj.cutoff Numeric value indicating the adjusted p-value threshold used to consider a geneset as significantly enriched. Default: \code{0.05}.
#' @param top.by String indicating the metric used to select the top genesets of each side. One among: \code{"significance"}, \code{"value"}. Default: \code{"significance"}.
#' @param pos.color String indicating the color of the bars of the first enrichment (positive side). Default: \code{"steelblue"}.
#' @param neg.color String indicating the color of the bars of the second enrichment (negative side). Default: \code{"orange"}.
#' @param alpha.range Numeric vector of length 2 indicating minimum and maximum value for the transparency, which is proportional to the significance. Individual values must be a number between 0 and 1. Default: \code{c(0.3, 1)}.
#' @param add.counts Logic value indicating whether labels with the counts of proteins ('count/geneset size') should be added at the end of each bar. Default: \code{TRUE}.
#' @param string.pattern.to.remove String with a regular expression of a pattern to be removed from the geneset names. Default: \code{NULL} (no changes).
#' @param max.term.length Numeric value indicating the maximal number of characters of the geneset names displayed on the y-axis (longer names are truncated). Default: \code{60}.
#' @param perc.bleeding.x Numeric value indicating the percentage of the full x-axis range to add on the left and on the right. Useful when labels are falling outside the x-max. Default: \code{8} (\%).
#' @param axes.text.size Numeric value indicating the font size of the axis text. Default: \code{10}.
#' @param bar.width Numeric value indicating the width of the bars. Default: \code{0.8}.
#' @param title String indicating the title of the plot, markdown formatting is supported. Default: \code{NULL} (built from the names of the list).
#'
#' @return A list with two elements: \code{results} (data.frame combining the results of the two enrichments) and \code{divergent.plot} (ggplot object).
#'
#' @seealso \code{\link{geneset.enrichment}}, \code{\link{combine.enrichments}}, \code{\link{NES.plot}}
#'
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#' @importFrom utils head
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' # Plot together the ORA of the two sides of the same contrast
#' \dontrun{
#' divergent <- divergent.enrichment(enrichment.list = list(`up-regulated` = ora.up,
#'                                                          `down-regulated` = ora.down),
#'                                   value = "FoldEnrichment",
#'                                   top.n = 10)
#' divergent$divergent.plot
#' }
#'
#' @export divergent.enrichment

divergent.enrichment =
  function(enrichment.list,
           value = "FoldEnrichment",
           top.n = 10,
           terms = NULL,
           padj.cutoff = 0.05,
           top.by = "significance",
           pos.color = "steelblue",
           neg.color = "orange",
           alpha.range = c(0.3, 1),
           add.counts = TRUE,
           string.pattern.to.remove = NULL,
           max.term.length = 60,
           perc.bleeding.x = 8,
           axes.text.size = 10,
           bar.width = 0.8,
           title = NULL) {

    ######################################################################################
    ### Checks
    ######################################################################################

    if (length(enrichment.list) != 2) {
      stop("The 'enrichment.list' must contain exactly two enrichments: use `combine.enrichments()` to display more than two discoveries.")
    }

    if (!(tolower(top.by) %in% c("significance", "value"))) {
      stop("The 'top.by' must be either 'significance' or 'value'.")
    }

    ## the name of the column corresponding to the metric asked
    value.column = switch(tolower(value),
                          "foldenrichment" = "FoldEnrichment",
                          "fold.enrichment" = "FoldEnrichment",
                          "generatio" = "GeneRatio.numeric",
                          "gene.ratio" = "GeneRatio.numeric",
                          "count" = "Count",
                          "nes" = "NES",
                          "padj" = "log10.padj",
                          "p.adjust" = "log10.padj",
                          NULL)

    if (is.null(value.column)) {
      stop("The 'value' must be one among: 'FoldEnrichment', 'GeneRatio', 'Count', 'NES', 'padj'.")
    }

    value.label = switch(value.column,
                         "FoldEnrichment" = "fold enrichment",
                         "GeneRatio.numeric" = "gene ratio",
                         "Count" = "protein count",
                         "NES" = "NES",
                         "log10.padj" = "-log~10~(P~adj~)")


    ######################################################################################
    ### Collection of the results
    ######################################################################################

    results = .collect.enrichments(enrichment.list = enrichment.list,
                                   string.pattern.to.remove = string.pattern.to.remove)

    if (length(levels(results$discovery)) != 2) {
      stop("The results of both the enrichments are required: one of the two could not be used.")
    }

    results$log10.padj = -log10(results$p.adjust)
    results$value = results[,value.column]

    ## the sign is what makes the plot divergent: first element up, second element down
    results$side.sign = ifelse(as.numeric(results$discovery) == 1, yes = 1, no = -1)

    ## the NES is already signed: only its magnitude is used, the side being given by the discovery
    results$signed.value = abs(results$value) * results$side.sign

    if (all(is.na(results$value))) {
      stop(paste0("The metric '", value, "' is not available for these enrichments."))
    }


    ######################################################################################
    ### Selection of the genesets to display
    ######################################################################################

    signif = dplyr::filter(.data = results, p.adjust <= padj.cutoff, !is.na(value))

    if (nrow(signif) == 0) {
      warning("No geneset passed the significance threshold: the plot is not generated.")
      return(list(results = results,
                  divergent.plot = NULL))
    }

    if (!is.null(terms)) {
      bar.tb = dplyr::filter(.data = signif, ID %in% terms | alias %in% terms)

      if (nrow(bar.tb) == 0) {
        stop("None of the 'terms' provided could be found in the enrichments.")
      }
    } else {
      ## the top genesets are taken independently for each side, so that both are represented
      ## even when one of the two is much more enriched than the other
      bar.tb =
        dplyr::bind_rows(lapply(split(signif, signif$discovery),
                                function(x) {
                                  if (tolower(top.by) == "value") {
                                    utils::head(x[order(-abs(x$value)),], top.n)
                                  } else {
                                    utils::head(x[order(x$p.adjust),], top.n)
                                  }
                                }))
    }


    ######################################################################################
    ### Order of the genesets on the y-axis
    ######################################################################################

    ## genesets found on both sides keep a single row, in which the two bars face each other:
    ## the sum of the signed values places them in the middle of the plot
    term.order = dplyr::summarise(dplyr::group_by(bar.tb, alias),
                                  order.value = sum(signed.value, na.rm = TRUE),
                                  .groups = "drop")
    term.order = term.order[order(term.order$order.value),]$alias

    bar.tb$alias = factor(bar.tb$alias, levels = term.order)


    ######################################################################################
    ### Divergent bar plot
    ######################################################################################

    ### colors
    side.colors = c(pos.color, neg.color)
    names(side.colors) = levels(results$discovery)

    if (is.null(title)) {
      title = paste0("**", levels(results$discovery)[1], "** *vs* **", levels(results$discovery)[2], "**")
    }


    divergent.plot =
      ggplot2::ggplot(data = bar.tb,
                      ggplot2::aes(x = signed.value,
                                   y = alias,
                                   fill = discovery)) +
      ggplot2::geom_bar(ggplot2::aes(alpha = log10.padj),
                        stat = "identity",
                        show.legend = TRUE,
                        width = bar.width) +
      ggplot2::geom_vline(xintercept = 0, color = "black") +
      ggplot2::scale_alpha_continuous(range = alpha.range, name = "-log~10~(P~adj~)") +
      ggplot2::scale_fill_manual(values = side.colors, name = NULL, drop = FALSE)


    ### counts written at the end of each bar, outside of it
    if (add.counts == TRUE) {
      bar.tb$count.label = ifelse(is.na(bar.tb$set.size),
                                  yes = as.character(bar.tb$Count),
                                  no = paste0(bar.tb$Count, "/", bar.tb$set.size))

      if (nrow(dplyr::filter(.data = bar.tb, signed.value >= 0)) > 0) {
        divergent.plot =
          divergent.plot +
          ggplot2::geom_text(data = dplyr::filter(.data = bar.tb, signed.value >= 0),
                             ggplot2::aes(x = signed.value, y = alias, label = count.label),
                             inherit.aes = FALSE,
                             size = axes.text.size / 3.5,
                             vjust = 0.5,
                             hjust = -0.2)
      }

      if (nrow(dplyr::filter(.data = bar.tb, signed.value < 0)) > 0) {
        divergent.plot =
          divergent.plot +
          ggplot2::geom_text(data = dplyr::filter(.data = bar.tb, signed.value < 0),
                             ggplot2::aes(x = signed.value, y = alias, label = count.label),
                             inherit.aes = FALSE,
                             size = axes.text.size / 3.5,
                             vjust = 0.5,
                             hjust = 1.2)
      }
    } else {
      perc.bleeding.x = 0
    }


    ### the axis is symmetric and the negative side is labeled with absolute values
    max.abs.value = max(abs(bar.tb$signed.value), na.rm = TRUE)

    divergent.plot =
      divergent.plot +
      ggplot2::scale_x_continuous(limits = c(-1,1) * (max.abs.value + (perc.bleeding.x/100) * max.abs.value),
                                  labels = function(x){abs(x)},
                                  expand = c(0,0)) +
      ggplot2::scale_y_discrete(labels = function(x) {
        ifelse(nchar(x) > max.term.length, paste0(substr(x, 1, max.term.length), "..."), x)
      }) +
      ggplot2::labs(x = value.label,
                    y = NULL,
                    title = title) +
      ggpubr::theme_pubr(legend = "right") +
      ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.text = ggplot2::element_text(color = "black", size = axes.text.size),
                     axis.title.x = ggtext::element_markdown(),
                     legend.title = ggtext::element_markdown(),
                     panel.background = ggplot2::element_blank(),
                     plot.title = ggtext::element_markdown(hjust = 0.5))


    ######################################################################################
    ### Output
    ######################################################################################

    return(list(results = results,
                divergent.plot = divergent.plot))

  } # END function
