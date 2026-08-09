##########################################
### INTERNAL FUNCTIONS :: ENRICHMENTS  ###
##########################################
#
#   .parse.ratio()           converts the 'x/y' strings of clusterProfiler into a numeric ratio
#   .get.enrichment.table()  extracts an harmonized results table from any enrichment object
#   .collect.enrichments()   applies .get.enrichment.table() over a named list of enrichments
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
