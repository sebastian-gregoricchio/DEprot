#' @title combine.enrichments
#'
#' @description Combines the results of multiple enrichment analyses, performed independently, in a single dotplot: the discoveries are displayed on the x-axis and the enriched genesets on the y-axis. The size of each dot reflects the enrichment (fold enrichment or gene ratio), its color the significance, and the number written inside the dot the count of proteins found in that geneset. The input enrichments can come from different sources (DEprot or clusterProfiler) and can be mixed.
#'
#' @param enrichment.list Named list of enrichments. Each element can be an object of class \code{DEprot.enrichResult} (DEprot), \code{DEprot.timecourse.enrichment} (DEprot), \code{enrichResult} or \code{gseaResult} (clusterProfiler), or a data.frame with the same structure of the tables returned by clusterProfiler. The names of the list are used as labels of the x-axis and define their order. No default.
#' @param dotplot.n Numeric value indicating the maximum number of genesets to keep for each discovery. The genesets displayed correspond to the union of these top hits. Default: \code{5}.
#' @param terms Character vector indicating specific geneset IDs to display. When provided, \code{dotplot.n} is ignored. Default: \code{NULL}.
#' @param padj.cutoff Numeric value indicating the adjusted p-value threshold used to consider a geneset as significantly enriched. Default: \code{0.05}.
#' @param size.by String indicating the metric used for the size of the dots. One among: \code{"FoldEnrichment"}, \code{"GeneRatio"}. Default: \code{"FoldEnrichment"}.
#' @param order.by String indicating how the genesets are sorted on the y-axis. One among: \code{"discovery"} (by the discovery in which each geneset is the most significant, generating a diagonal pattern), \code{"significance"} (by the best adjusted p-value), \code{"clustering"} (hierarchical clustering of the enrichment values), \code{"alphabetical"}. Default: \code{"discovery"}.
#' @param show.non.significant Logic value indicating whether the genesets that were tested but did not pass the \code{padj.cutoff} in a given discovery should be displayed as grey dots. Default: \code{TRUE}.
#' @param show.numbers Logic value indicating whether the number of proteins should be written inside the dots. Default: \code{TRUE}.
#' @param string.pattern.to.remove String with a regular expression of a pattern to be removed from the geneset names. Default: \code{NULL} (no changes).
#' @param size.range Numeric vector of length 2 indicating the minimum and maximum size of the dots. Default: \code{c(5, 14)}.
#' @param number.size Numeric value indicating the font size of the numbers written inside the dots. Default: \code{2.7}.
#' @param max.term.length Numeric value indicating the maximal number of characters of the geneset names displayed on the y-axis (longer names are truncated). Default: \code{60}.
#' @param viridis.option String indicating the viridis palette used for the significance color scale. Default: \code{"rocket"}.
#' @param title String indicating the title of the plot, markdown formatting is supported. Default: \code{"**Combined enrichments**"}.
#'
#' @return A list with two elements: \code{results} (data.frame combining the results of all the discoveries) and \code{dotplot} (ggplot object).
#'
#' @seealso \code{\link{geneset.enrichment}}, \code{\link{divergent.enrichment}}, \code{\link{timecourse.enrichment}}
#'
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#' @import viridis
#' @importFrom reshape2 acast
#' @importFrom stats dist hclust
#' @importFrom utils head
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' # Combine the ORA of two different contrasts
#' \dontrun{
#' combined <- combine.enrichments(enrichment.list = list(`6h E2` = ora.6h,
#'                                                        `24h E2` = ora.24h),
#'                                 dotplot.n = 5,
#'                                 size.by = "FoldEnrichment")
#' combined$dotplot
#' }
#'
#' @export combine.enrichments

combine.enrichments =
  function(enrichment.list,
           dotplot.n = 5,
           terms = NULL,
           padj.cutoff = 0.05,
           size.by = "FoldEnrichment",
           order.by = "discovery",
           show.non.significant = TRUE,
           show.numbers = TRUE,
           string.pattern.to.remove = NULL,
           size.range = c(5, 14),
           number.size = 2.7,
           max.term.length = 60,
           viridis.option = "rocket",
           title = "**Combined enrichments**") {

    ######################################################################################
    ### Checks
    ######################################################################################

    if (!(tolower(size.by) %in% c("foldenrichment", "fold.enrichment", "generatio", "gene.ratio"))) {
      stop("The 'size.by' must be either 'FoldEnrichment' or 'GeneRatio'.")
    } else {
      size.by = ifelse(grepl("fold", tolower(size.by)), yes = "FoldEnrichment", no = "GeneRatio")
    }

    if (!(tolower(order.by) %in% c("discovery", "significance", "clustering", "alphabetical"))) {
      stop("The 'order.by' must be one among: 'discovery', 'significance', 'clustering', 'alphabetical'.")
    }


    ######################################################################################
    ### Collection of the results
    ######################################################################################

    results = .collect.enrichments(enrichment.list = enrichment.list,
                                   string.pattern.to.remove = string.pattern.to.remove)

    results$significant = results$p.adjust <= padj.cutoff
    results$size.value = if (size.by == "FoldEnrichment") {results$FoldEnrichment} else {results$GeneRatio.numeric}

    ## the fold enrichment is not defined for a GSEA: those discoveries would silently disappear
    if (size.by == "FoldEnrichment" & any(is.na(results$size.value))) {
      missing.discoveries = unique(as.character(results$discovery[is.na(results$size.value)]))
      warning(paste0("The FoldEnrichment is not available for the following discoveries: ",
                     paste(missing.discoveries, collapse = ", "),
                     ". Use size.by = 'GeneRatio' to include them."))
    }


    ######################################################################################
    ### Selection of the genesets to display
    ######################################################################################

    signif = dplyr::filter(.data = results, significant == TRUE)

    if (nrow(signif) == 0) {
      warning("No geneset passed the significance threshold: the dotplot is not generated.")
      return(list(results = results,
                  dotplot = NULL))
    }

    if (!is.null(terms)) {
      selected.terms = unique(results$ID[results$ID %in% terms | results$alias %in% terms])

      if (length(selected.terms) == 0) {
        stop("None of the 'terms' provided could be found in the enrichments.")
      }
    } else {
      ## top genesets of EACH discovery, then their union: this keeps the plot readable while
      ## still showing, for every geneset, in which other discoveries it also came up
      selected.terms = unique(unlist(lapply(split(signif, signif$discovery),
                                            function(x) {utils::head(x[order(x$p.adjust),]$ID, dotplot.n)})))
    }

    dot.tb = dplyr::filter(.data = results, ID %in% selected.terms)

    if (show.non.significant == FALSE) {
      dot.tb = dplyr::filter(.data = dot.tb, significant == TRUE)
    }

    dot.tb = dplyr::filter(.data = dot.tb, !is.na(size.value))


    ######################################################################################
    ### Order of the genesets on the y-axis
    ######################################################################################

    ## the ordering is always computed on the significant hits: the grey dots must not move the rows
    order.tb = dplyr::filter(.data = dot.tb, significant == TRUE)

    if (tolower(order.by) == "discovery") {
      term.order = dplyr::summarise(dplyr::group_by(order.tb, alias),
                                    best.discovery = as.numeric(discovery)[which.min(p.adjust)],
                                    best.padj = min(p.adjust, na.rm = TRUE),
                                    .groups = "drop")
      term.order = term.order[order(-term.order$best.discovery, term.order$best.padj),]$alias

    } else if (tolower(order.by) == "significance") {
      term.order = dplyr::summarise(dplyr::group_by(order.tb, alias),
                                    best.padj = min(p.adjust, na.rm = TRUE),
                                    .groups = "drop")
      term.order = term.order[order(-term.order$best.padj),]$alias

    } else if (tolower(order.by) == "clustering" & length(unique(dot.tb$alias)) > 2) {
      ## genesets showing a similar behavior among the discoveries are placed next to each other
      enrichment.matrix = reshape2::acast(data = dot.tb,
                                          formula = alias ~ discovery,
                                          value.var = "size.value",
                                          fun.aggregate = function(x){mean(x, na.rm = TRUE)},
                                          fill = 0)
      enrichment.matrix[is.na(enrichment.matrix)] = 0

      term.order = rownames(enrichment.matrix)[stats::hclust(stats::dist(enrichment.matrix))$order]

    } else {
      term.order = rev(sort(unique(dot.tb$alias)))
    }

    dot.tb$alias = factor(dot.tb$alias, levels = term.order)
    dot.tb = dplyr::filter(.data = dot.tb, !is.na(alias))

    dot.tb$log10.padj = -log10(dot.tb$p.adjust)


    ######################################################################################
    ### Dotplot
    ######################################################################################

    ## the numbers must stay readable on both the dark and the light end of the palette
    mid.point = mean(range(dot.tb$log10.padj[dot.tb$significant == TRUE], na.rm = TRUE))
    dot.tb$label.color = ifelse(dot.tb$significant == FALSE, yes = "grey40",
                                no = ifelse(dot.tb$log10.padj >= mid.point, yes = "white", no = "black"))

    dotplot.combined =
      ggplot2::ggplot(data = dot.tb,
                      ggplot2::aes(x = discovery, y = alias))

    ## the genesets tested but not enriched are shown as empty dots, so that a missing dot
    ## keeps the meaning of "the geneset was not tested/returned at all"
    if (nrow(dplyr::filter(.data = dot.tb, significant == FALSE)) > 0) {
      dotplot.combined =
        dotplot.combined +
        ggplot2::geom_point(data = dplyr::filter(.data = dot.tb, significant == FALSE),
                            ggplot2::aes(size = size.value),
                            shape = 21,
                            fill = "grey92",
                            color = "grey70",
                            stroke = 0.3)
    }

    dotplot.combined =
      dotplot.combined +
      ggplot2::geom_point(data = dplyr::filter(.data = dot.tb, significant == TRUE),
                          ggplot2::aes(size = size.value, fill = log10.padj),
                          shape = 21,
                          color = "grey30",
                          stroke = 0.3)

    ## the count of proteins is written in the middle of each dot
    if (show.numbers == TRUE) {
      dotplot.combined =
        dotplot.combined +
        ggplot2::geom_text(ggplot2::aes(label = Count, color = label.color),
                           size = number.size,
                           fontface = "bold") +
        ggplot2::scale_color_identity()
    }

    dotplot.combined =
      dotplot.combined +
      ggplot2::scale_size_continuous(range = size.range, name = size.by) +
      viridis::scale_fill_viridis(option = viridis.option, direction = -1, begin = 0.15, end = 0.9,
                                  name = "-log~10~(P~adj~)") +
      ggplot2::scale_y_discrete(labels = function(x) {
        ifelse(nchar(x) > max.term.length, paste0(substr(x, 1, max.term.length), "..."), x)
      }) +
      ggplot2::labs(x = NULL,
                    y = NULL,
                    title = title,
                    caption = paste0("The number inside each dot is the count of proteins in the geneset (P~adj~ \u2264 ",
                                     padj.cutoff, ")")) +
      ggpubr::theme_pubr(legend = "right") +
      ggplot2::theme(plot.title = ggtext::element_markdown(hjust = 0.5),
                     plot.caption = ggtext::element_markdown(size = 7, color = "grey30"),
                     legend.title = ggtext::element_markdown(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "gray80"),
                     panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1),
                     axis.line = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))


    ######################################################################################
    ### Output
    ######################################################################################

    return(list(results = results,
                dotplot = dotplot.combined))

  } # END function
