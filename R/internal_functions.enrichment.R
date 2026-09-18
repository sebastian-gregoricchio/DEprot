##########################################
### INTERNAL FUNCTIONS :: ENRICHMENTS  ###
##########################################
#
#   .parse.ratio()           converts the 'x/y' strings of clusterProfiler into a numeric ratio
#   .filter.enrichment()     applies the significance thresholds to the results of an enrichment
#   .get.enrichment.table()  extracts an harmonized results table from any enrichment object
#   .collect.enrichments()   applies .get.enrichment.table() over a named list of enrichments
#   .scale.alpha.padj()      builds the transparency scale used for the adjusted p-values
#
##########################################


# ----------------------------------------------------------------------------------------

#' @title .parse.ratio
#'
#' @description Internal. Converts the 'x/y' strings returned by clusterProfiler into a numeric ratio.
#'
#' @param x Vector of strings in the form 'x/y'.
#' @param position String indicating which value to return: \code{"ratio"} (x/y), \code{"numerator"} (x) or \code{"denominator"} (y). Default: \code{"ratio"}.
#'
#' @return A numeric vector.
#'
#' @keywords internal

.parse.ratio =
  function(x,
           position = "ratio") {

    split.values = strsplit(as.character(x), "/")

    if (tolower(position) %in% c("num", "numerator")) {
      values = vapply(split.values, function(v){as.numeric(v[1])}, FUN.VALUE = numeric(1))
    } else if (tolower(position) %in% c("den", "denominator")) {
      values = vapply(split.values, function(v){as.numeric(v[2])}, FUN.VALUE = numeric(1))
    } else {
      values = vapply(split.values, function(v){as.numeric(v[1]) / as.numeric(v[2])}, FUN.VALUE = numeric(1))
    }

    return(values)

  } # END .parse.ratio



# ----------------------------------------------------------------------------------------

#' @title .filter.enrichment
#'
#' @description Internal. Applies the significance thresholds to the results of an enrichment discovery. The GSEA function of clusterProfiler does not provide a 'qvalueCutoff' option and, depending on the version installed, it does not always apply the 'pvalueCutoff' to the adjusted p-values either. The thresholds are therefore re-applied here, so that the results returned by DEprot do not depend on the version of clusterProfiler used.
#'
#' @param enrichment An object of class \code{gseaResult} or \code{enrichResult} (clusterProfiler), or \code{NULL} when the analyses could not be performed.
#' @param pvalueCutoff Numeric value indicating the threshold applied to both the uncorrected and the adjusted p-values. Default: \code{0.05}.
#' @param qvalueCutoff Numeric value indicating the threshold applied to the q-values. Default: \code{0.05}.
#'
#' @return The input object with a filtered results table. Genesets for which the q-value could not be estimated (NA) are kept, while \code{NULL} inputs are returned as such.
#'
#' @keywords internal

.filter.enrichment =
  function(enrichment,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05) {

    ### the discovery is NULL when the enrichment failed: there is nothing to filter
    if (is.null(enrichment)) {return(NULL)}

    if (!(methods::is(enrichment, "gseaResult") | methods::is(enrichment, "enrichResult"))) {
      return(enrichment)
    }

    results = enrichment@result

    if (is.null(results)) {return(enrichment)}
    if (nrow(results) == 0) {return(enrichment)}


    ### each threshold is applied only if the corresponding column is available
    keep = rep(TRUE, nrow(results))

    if ("pvalue" %in% colnames(results)) {
      keep = keep & !is.na(results$pvalue) & results$pvalue <= pvalueCutoff
    }

    if ("p.adjust" %in% colnames(results)) {
      keep = keep & !is.na(results$p.adjust) & results$p.adjust <= pvalueCutoff
    }

    ## the q-value is NA when its estimation fails (too few genesets tested): these are kept
    if ("qvalue" %in% colnames(results)) {
      keep = keep & (is.na(results$qvalue) | results$qvalue <= qvalueCutoff)
    }


    ### base subsetting is used to keep the rownames, which clusterProfiler uses as geneset IDs
    enrichment@result = results[keep, , drop = FALSE]

    return(enrichment)

  } # END .filter.enrichment



# ----------------------------------------------------------------------------------------

#' @title .get.enrichment.table
#'
#' @description Internal. Extracts the results of a single enrichment discovery and returns them with a common set of column names, independently of the tool and of the type of analyses used to generate them. GSEA and ORA results do not share the same columns: the leading edge of a GSEA is used as the equivalent of the ORA 'Count', and the size of the geneset as its background.
#'
#' @param enrichment An object of class \code{DEprot.enrichResult} (DEprot), \code{enrichResult} or \code{gseaResult} (clusterProfiler), or a data.frame with the same structure of the ones returned by clusterProfiler.
#' @param name String indicating the name of the discovery, used to fill the 'discovery' column. Default: \code{NA}.
#' @param string.pattern.to.remove String with a regular expression of a pattern to be removed from the geneset names when generating the 'alias' column. Default: \code{NULL} (no changes).
#'
#' @return A data.frame with the columns: discovery, enrichment.type, ID, alias, Count, set.size, GeneRatio.numeric, FoldEnrichment, NES, pvalue, p.adjust. \code{NULL} when the enrichment is empty or of an unsupported class.
#'
#' @keywords internal

.get.enrichment.table =
  function(enrichment,
           name = NA,
           string.pattern.to.remove = NULL) {

    ### the DEprot wrappers are unwrapped first, so that only the clusterProfiler
    ### objects and the data.frames must be handled hereafter
    if ("DEprot.enrichResult" %in% class(enrichment)) {
      enrichment = enrichment@enrichment.discovery
    }

    if (is.null(enrichment)) {
      warning(paste0("The enrichment '", name, "' is empty: it will be skipped."))
      return(NULL)
    }


    ### collect the table and the type of analyses
    if ("gseaResult" %in% class(enrichment)) {
      tb = as.data.frame(enrichment@result)
      enrichment.type = "GSEA"
    } else if ("enrichResult" %in% class(enrichment)) {
      tb = as.data.frame(enrichment@result)
      enrichment.type = "ORA"
    } else if ("data.frame" %in% class(enrichment)) {
      tb = as.data.frame(enrichment)
      enrichment.type = ifelse("NES" %in% colnames(tb), yes = "GSEA", no = "ORA")
    } else {
      warning(paste0("The enrichment '", name, "' is of class '", paste(class(enrichment), collapse = "/"),
                     "', which is not supported: it will be skipped."))
      return(NULL)
    }


    if (nrow(tb) == 0) {
      warning(paste0("The enrichment '", name, "' does not contain any geneset: it will be skipped."))
      return(NULL)
    }

    if (!all(c("ID", "p.adjust") %in% colnames(tb))) {
      warning(paste0("The results of the enrichment '", name, "' do not contain the columns 'ID' and 'p.adjust': it will be skipped."))
      return(NULL)
    }

    if (!("pvalue" %in% colnames(tb))) {tb$pvalue = NA}


    ### harmonization of the metrics
    if (enrichment.type == "GSEA") {
      ## the leading edge is the closest equivalent of the ORA 'Count'
      tb$Count = lengths(strsplit(as.character(tb$core_enrichment), "/"))
      tb$set.size = tb$setSize
      tb$GeneRatio.numeric = tb$Count / tb$setSize
      tb$FoldEnrichment = NA

    } else {
      if (!("GeneRatio.numeric" %in% colnames(tb))) {
        tb$GeneRatio.numeric = .parse.ratio(tb$GeneRatio)
      }

      ## the older versions of clusterProfiler do not return the FoldEnrichment column
      if (!("FoldEnrichment" %in% colnames(tb))) {
        tb$FoldEnrichment = tb$GeneRatio.numeric / .parse.ratio(tb$BgRatio)
      }

      if ("BgRatio" %in% colnames(tb)) {
        tb$set.size = .parse.ratio(tb$BgRatio, position = "numerator")
      } else {
        tb$set.size = NA
      }

      tb$NES = NA
    }


    ### the geneset names are cleaned only for the display: the 'ID' is kept as key
    tb$alias = tb$ID

    if (!is.null(string.pattern.to.remove)) {
      tb$alias = gsub("_", " ", gsub(string.pattern.to.remove, "", tb$alias))
    }

    tb$discovery = name
    tb$enrichment.type = enrichment.type


    results = tb[,c("discovery", "enrichment.type", "ID", "alias", "Count", "set.size",
                    "GeneRatio.numeric", "FoldEnrichment", "NES", "pvalue", "p.adjust")]
    rownames(results) = NULL

    return(results)

  } # END .get.enrichment.table



# ----------------------------------------------------------------------------------------

#' @title .collect.enrichments
#'
#' @description Internal. Loops \code{.get.enrichment.table} over a named list of enrichments and combines the results in a single table. Objects of class \code{DEprot.timecourse.enrichment} carry several discoveries at once: each of their clusters is expanded into an independent discovery.
#'
#' @param enrichment.list Named list of enrichment objects.
#' @param string.pattern.to.remove String with a regular expression of a pattern to be removed from the geneset names. Default: \code{NULL} (no changes).
#'
#' @return A data.frame in which the 'discovery' column is a factor following the order of the input list.
#'
#' @import dplyr
#'
#' @keywords internal

.collect.enrichments =
  function(enrichment.list,
           string.pattern.to.remove = NULL) {

    if (!("list" %in% class(enrichment.list))) {
      stop("The 'enrichment.list' must be a (possibly named) list of enrichment objects.")
    }

    if (length(enrichment.list) == 0) {
      stop("The 'enrichment.list' is empty.")
    }

    ### unnamed elements would collapse into a single column of the plot: they get a generic name
    discovery.names = names(enrichment.list)

    if (is.null(discovery.names)) {
      discovery.names = paste0("enrichment.", 1:length(enrichment.list))
      warning("The 'enrichment.list' is not named: generic names have been assigned.")
    } else if (any(discovery.names %in% c("", NA))) {
      missing.names = which(discovery.names %in% c("", NA))
      discovery.names[missing.names] = paste0("enrichment.", missing.names)
      warning("Some elements of the 'enrichment.list' are not named: generic names have been assigned.")
    }


    ### collection of the single tables
    tables.list = list()

    for (i in 1:length(enrichment.list)) {

      enrichment = enrichment.list[[i]]

      if ("DEprot.timecourse.enrichment" %in% class(enrichment)) {

        if (is.null(enrichment@results)) {
          warning(paste0("The enrichment '", discovery.names[i], "' is empty: it will be skipped."))
          next
        }

        for (k in sort(unique(enrichment@results$cluster))) {
          cluster.name = paste0(discovery.names[i], ".cluster.", k)

          tables.list[[cluster.name]] =
            .get.enrichment.table(enrichment = enrichment@results[enrichment@results$cluster == k,],
                                  name = cluster.name,
                                  string.pattern.to.remove = string.pattern.to.remove)
        }

      } else {
        tables.list[[discovery.names[i]]] =
          .get.enrichment.table(enrichment = enrichment,
                                name = discovery.names[i],
                                string.pattern.to.remove = string.pattern.to.remove)
      }
    }


    if (length(tables.list) == 0) {
      stop("None of the elements of the 'enrichment.list' could be used.")
    }

    results = dplyr::bind_rows(tables.list)

    ## the order of the input list defines the order of the discoveries in the plots
    results$discovery = factor(results$discovery, levels = names(tables.list))
    rownames(results) = NULL

    return(results)

  } # END .collect.enrichments



# ----------------------------------------------------------------------------------------

#' @title .scale.alpha.padj
#'
#' @description Internal. Builds the ggplot2 transparency scale used to display the adjusted p-values. The transparency is mapped on the \code{-log10} of the adjusted p-value, hence the breaks are spaced logarithmically and the most opaque bars correspond to the most significant genesets, but the labels of the legend show the p-values themselves and not their logarithm.
#'
#' @param padj Numeric vector of the adjusted p-values displayed, used to define the breaks.
#' @param alpha.range Numeric vector of length 2 indicating minimum and maximum value for the transparency. Default: \code{c(0.3, 1)}.
#' @param name String indicating the title of the legend. Markdown is supported as long as the theme of the plot defines \code{legend.title = ggtext::element_markdown()}. Default: \code{"P~adj~"}.
#' @param max.breaks Numeric value indicating the maximum number of breaks displayed in the legend. Default: \code{5}.
#'
#' @return A ggplot2 continuous alpha scale, to be combined with an \code{aes(alpha = -log10(p.adjust))} mapping.
#'
#' @import ggplot2
#'
#' @keywords internal

.scale.alpha.padj =
  function(padj,
           alpha.range = c(0.3,1),
           name = "P~adj~",
           max.breaks = 5) {

    ### only positive and finite values can be placed on a log scale
    values = padj[is.finite(padj) & padj > 0]

    if (length(values) == 0) {
      return(ggplot2::scale_alpha_continuous(range = alpha.range, name = name))
    }

    log.range = range(-log10(values), na.rm = TRUE)


    ### the breaks are searched among 'round' p-values: the decades alone when the values
    ### span more than ~1.5 orders of magnitude, their halves and fifths otherwise
    decades = 10^-(0:20)

    if (diff(log.range) >= 1.5) {
      candidates = decades
    } else {
      candidates = sort(unique(as.vector(outer(c(1,5,2), decades))), decreasing = TRUE)
      candidates = candidates[candidates <= 1]
    }

    breaks = candidates[-log10(candidates) >= log.range[1] & -log10(candidates) <= log.range[2]]

    if (length(breaks) > max.breaks) {
      breaks = breaks[seq(1, length(breaks), by = ceiling(length(breaks)/max.breaks))]
    }

    ## no round p-value falls in the range when all the values are very close to each other:
    ## in this case the range itself is split in equally spaced points on the log scale
    if (length(breaks) < 2) {
      breaks = signif(10^-seq(log.range[1], log.range[2], length.out = 3), 2)
    }

    breaks = sort(unique(breaks), decreasing = TRUE)


    ### the scale is built on the -log10 of the p-values, the labels show the p-values
    labels = vapply(X = breaks,
                    FUN = function(x){
                      if (x >= 1e-4) {
                        format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
                      } else {
                        format(x, scientific = TRUE, digits = 1, trim = TRUE)
                      }},
                    FUN.VALUE = character(1))

    ## the legend is reversed to show the most significant (most opaque) values on the top
    alpha.scale =
      ggplot2::scale_alpha_continuous(range = alpha.range,
                                      breaks = -log10(breaks),
                                      labels = labels,
                                      name = name,
                                      guide = ggplot2::guide_legend(reverse = TRUE))

    return(alpha.scale)

  } # END .scale.alpha.padj
