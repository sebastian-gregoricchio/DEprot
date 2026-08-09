#' @title missingness.diagnostic
#'
#' @description Diagnostics of the missing-value structure of a proteomics dataset, aimed at
#' discriminating values missing not-at-random (MNAR, left-censored, below detection limit) from values
#' missing (completely) at random (MCAR), before/while choosing an imputation strategy.
#' The classification of the missing values follows the same logic of the \code{DEprot} "double-imputation"
#' strategy: values missing in at least \code{percentage.missing}\% of the replicates of a group are considered
#' MNAR (and would be replaced by \link{randomize.missing.values} using the bottom \code{tail.percentage}\% of the
#' distribution), while the remaining missing values are considered MCAR (and would be replaced by \link{impute.counts}).
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param group.column String indicating the metadata column defining the groups of replicates within which
#' the missing values are counted. Default: \code{NULL} (retrieved from the randomization parameters stored in the
#' object; if not available, from the first contrast of a \code{DEprot.analyses} object).
#' @param which.data String indicating which counts should be used. One among: 'auto', 'raw', 'normalized',
#' 'randomized', 'imputed'. When \code{'auto'} (default) the \emph{lowest} available level is used, following the
#' priority raw > normalized > randomized > imputed. Default: \code{"auto"}.
#' @param percentage.missing Numeric value between 0 and 100 indicating the minimal percentage of missing values
#' per group required to consider a protein MNAR-like in that group. Default: \code{NULL} (retrieved from the
#' randomization parameters if available, otherwise \code{100}, the default of \link{randomize.missing.values}).
#' @param tail.percentage Numeric value between 0 and 100 (excluded) indicating the bottom percentage of the
#' distribution used by \link{randomize.missing.values} to sample the random values. Used here only to display the
#' corresponding intensity threshold on the diagnostic plots. Default: \code{NULL} (retrieved from the randomization
#' parameters if available, otherwise \code{3}).
#' @param contrasts Vector of contrast indexes (or names) for which a contrast-level diagnostic should be computed.
#' Only used when a \code{DEprot.analyses} object is provided. Use \code{NULL} for all the contrasts available and
#' \code{NA}/\code{"none"} to skip the contrast-level analyses. Default: \code{NULL} (all contrasts).
#' @param sample.subset String vector indicating the samples (\code{column.id}) to keep. Default: \code{NULL} (all samples).
#' @param min.detected Integer indicating the minimal number of samples in which a protein must be detected to be
#' included in the intensity-dependent diagnostics (density and dropout curve). Default: \code{1}.
#' @param max.proteins.heatmap Integer indicating the maximum number of proteins displayed in the missingness heatmap
#' (randomly sampled among the ones showing at least one missing value, in order to limit the clustering time).
#' Default: \code{2500}.
#' @param cluster.method String indicating the agglomeration method used by \code{stats::hclust} for the rows of the
#' missingness heatmap and for the clustering of the samples in the Jaccard similarity heatmap. Default: \code{"ward.D2"}.
#' @param colors Named string vector with the colors used for the missing-value classes. Names required:
#' 'detected', 'MNAR', 'MCAR', 'all.missing'. Default: \code{c(detected = "grey85", MNAR = "indianred", MCAR = "steelblue", all.missing = "grey25")}.
#' @param jaccard.palette Vector of colors corresponding to the palette to use for the color scale of the
#' sample-vs-sample Jaccard similarity heatmap. Default: \code{viridis::mako(100, direction = -1)}.
#' @param dendrogram.color String indicating the color of the dendrogram lines of the Jaccard similarity heatmap.
#' Default: \code{"black"}.
#' @param convert.log2 Logical value to define whether counts should be converted to log2. Default: \code{TRUE}.
#' @param seed Numeric value indicating the random seed used to subsample the proteins in the heatmap. Default: \code{20240101}.
#' @param verbose Logical value to define whether messages should be printed. Default: \code{TRUE}.
#'
#' @return An object of class \code{\link{DEprot.missingness}}.
#'
#' @details
#' Six families of diagnostics are computed:
#' \describe{
#'   \item{\code{detection.density}: }{densities of the average protein intensity, split by proteins fully quantified
#'   and proteins showing at least one missing value. Superimposed densities support MCAR, a left-shift of the
#'   incomplete proteins supports MNAR (left-censoring).}
#'   \item{\code{dropout.curve}: }{per-protein fraction of missing values as a function of the average protein
#'   intensity, with the fitted logistic dropout model. The intensity at which 50\% of the values are missing is
#'   reported as \code{LOD50}, i.e. an empirical estimate of the limit of detection.}
#'   \item{\code{missingness.heatmap}: }{binary protein-by-sample map in which each missing cell is colored according
#'   to the strategy that would be applied to it (randomization of the bottom distribution for MNAR-like cells,
#'   imputation for MCAR-like cells).}
#'   \item{\code{missing.per.sample} / \code{detection.frequency}: }{completeness summaries per sample and per protein.}
#'   \item{\code{sample.similarity}: }{Jaccard similarity of the detection patterns between samples (clustered, with
#'   dendrogram); a structure that follows the experimental groups/batches indicates that the missingness is not
#'   completely at random.}
#'   \item{\code{pattern.summary} / \code{upset}: }{number of proteins per missing-value class, globally and per group,
#'   and intersections of the MNAR-like patterns between groups.}
#' }
#'
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#' @importFrom reshape2 melt
#' @importFrom viridis viridis mako
#' @importFrom stats quantile dist as.dist hclust glm binomial median wilcox.test sd coef predict setNames
#' @importFrom colorspace lighten
#' @importFrom legendry scale_y_dendro
#' @importFrom grDevices pdf dev.off
#' @importFrom utils modifyList
#' @importFrom methods new slot
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' dpo <- load.counts2(counts = DEprot::unimputed.counts,
#'                     metadata = DEprot::sample.config,
#'                     log.base = 2,
#'                     data.type = "norm")
#'
#' miss <- missingness.diagnostic(DEprot.object = dpo,
#'                                group.column = "combined.id")
#' miss
#' plot(miss)
#'
#' @export missingness.diagnostic

missingness.diagnostic =
  function(DEprot.object,
           group.column = NULL,
           which.data = "auto",
           percentage.missing = NULL,
           tail.percentage = NULL,
           contrasts = NULL,
           sample.subset = NULL,
           min.detected = 1,
           max.proteins.heatmap = 2500,
           cluster.method = "ward.D2",
           colors = c(detected = "grey85",
                      MNAR = "indianred",
                      MCAR = "steelblue",
                      all.missing = "grey25"),
           jaccard.palette = viridis::mako(100, direction = -1),
           dendrogram.color = "black",
           convert.log2 = TRUE,
           seed = 20240101,
           verbose = TRUE) {


    ### Check object
    if (!any(c("DEprot", "DEprot.analyses") %in% class(DEprot.object))) {
      stop("The input must be an object of class 'DEprot' or 'DEprot.analyses'.")
    }

    ### Check colors
    required.colors = c("detected", "MNAR", "MCAR", "all.missing")
    if (!all(required.colors %in% names(colors))) {
      stop(paste0("The 'colors' vector must be a named vector including: ",
                  paste0("'", required.colors, "'", collapse = ", "), "."))
    }


    ##--------------------------------------------------------------------##
    ##  Selection of the counts table (lowest level available)            ##
    ##--------------------------------------------------------------------##
    counts.priority = c(raw = "raw.counts",
                        normalized = "norm.counts",
                        randomized = "random.counts",
                        imputed = "imputed.counts")

    available = names(counts.priority)[!vapply(X = counts.priority,
                                               FUN = function(s) {.deprot_slot_is_empty(methods::slot(DEprot.object, s))},
                                               FUN.VALUE = logical(1))]

    if (length(available) == 0) {
      stop("The 'DEprot' object does not contain any counts table.")
    }

    if (identical(available, "imputed")) {
      stop(paste0("Only IMPUTED counts are available in this object: no missing value is left to diagnose.\n",
                  "       Re-run the analyses starting from a 'DEprot' object containing raw and/or normalized counts."))
    }

    if (tolower(which.data) %in% c("auto", "automatic", "lowest")) {
      data.used = available[1]  # counts.priority is already ordered raw > norm > random > imputed
    } else if (tolower(which.data) %in% c("raw", "r")) {
      data.used = "raw"
    } else if (tolower(which.data) %in% c("norm", "normalized", "normal", "n")) {
      data.used = "normalized"
    } else if (tolower(which.data) %in% c("random", "randomized")) {
      data.used = "randomized"
    } else if (tolower(which.data) %in% c("imputed", "imp", "impute")) {
      data.used = "imputed"
    } else {
      stop(paste0("The 'which.data' value is not recognized.\n",
                  "       Please indicate a count type among 'auto', 'raw', 'normalized', 'randomized', 'imputed'."))
    }

    if (!(data.used %in% available)) {
      stop(paste0("Use of ", toupper(data.used), " counts was required, but not available.\n",
                  "       Counts available: ", paste(available, collapse = ", "), "."))
    }

    if (data.used %in% c("randomized", "imputed")) {
      warning(paste0("The ", toupper(data.used), " counts are used: part (or all) of the missing values have already been ",
                     "replaced and the diagnostic will underestimate the missingness."), call. = FALSE)
    }

    mat = as.matrix(methods::slot(DEprot.object, counts.priority[data.used]))

    if (verbose == TRUE) {
      message(paste0("Counts available: ", paste(available, collapse = ", "),
                     " | counts used: ", data.used, "."))
    }


    ##--------------------------------------------------------------------##
    ##  Metadata, samples and groups                                      ##
    ##--------------------------------------------------------------------##
    meta = as.data.frame(DEprot.object@metadata)

    if (!is.null(sample.subset)) {
      if (!all(sample.subset %in% meta$column.id)) {
        stop("Some of the samples indicated in 'sample.subset' are not available in the metadata table ('column.id').")
      }
      meta = meta[meta$column.id %in% sample.subset, , drop = FALSE]
      mat = mat[, meta$column.id, drop = FALSE]
    } else {
      meta = meta[match(colnames(mat), meta$column.id), , drop = FALSE]
    }


    ##--------------------------------------------------------------------##
    ##  Retrieval of the randomization parameters                         ##
    ##--------------------------------------------------------------------##
    rand.params = DEprot.object@randomization.method
    params.source = "user-defined"

    if (!.deprot_slot_is_empty(rand.params) & is.list(rand.params)) {
      params.source = "randomization"
      if (is.null(group.column))        {group.column = rand.params$group.column}
      if (is.null(percentage.missing))  {percentage.missing = rand.params$percentage.missing}
      if (is.null(tail.percentage))     {tail.percentage = rand.params$tail.percentage}
    }

    ## defaults of randomize.missing.values()
    if (is.null(percentage.missing)) {percentage.missing = 100 ; params.source = paste0(params.source, " (defaults)")}
    if (is.null(tail.percentage))    {tail.percentage = 3}

    ## group.column fallback: first contrast of a DEprot.analyses object
    if (is.null(group.column)) {
      if (!.deprot_slot_is_empty(DEprot.object@contrasts)) {
        group.column = DEprot.object@contrasts[[1]]$metadata.column
        if (verbose == TRUE) {
          message(paste0("No 'group.column' provided: the metadata column of the first contrast ('",
                         group.column, "') will be used."))
        }
      } else {
        stop(paste0("Impossible to define the groups of replicates.\n",
                    "       Provide a 'group.column', or run 'randomize.missing.values()' beforehand."))
      }
    }

    if (!(group.column %in% colnames(meta))) {
      stop("The 'group.column' provided is not available in the metadata column names.")
    }

    if (percentage.missing < 0 | percentage.missing > 100) {
      stop("The 'percentage.missing' must be a number between 0 and 100 (included).")
    }
    if (tail.percentage <= 0 | tail.percentage >= 100) {
      stop("The 'tail.percentage' must be a number between 0 and 100 (both excluded).")
    }

    if (verbose == TRUE) {
      message(paste0("MNAR definition: >= ", percentage.missing, "% of missing values within the groups of '",
                     group.column, "' (parameters: ", params.source, ")."))
    }


    ##--------------------------------------------------------------------##
    ##  Log2 conversion (only for the intensity-based diagnostics)        ##
    ##--------------------------------------------------------------------##
    log.base = DEprot.object@log.base

    if (convert.log2 == TRUE) {
      if (!is.numeric(log.base) | all(is.na(log.base))) {
        warning("The log.base is not numeric, linear counts are assumed: the matrix will be converted to log2(x+1).", call. = FALSE)
        mat.log2 = log2(mat + 1)
      } else if (as.numeric(log.base) != 2) {
        warning("The log.base is not 2: the counts will be converted to log2 values.", call. = FALSE)
        mat.log2 = log2(log.base^mat)
      } else {
        mat.log2 = mat
      }
    } else {
      mat.log2 = mat
    }


    ##--------------------------------------------------------------------##
    ##  Cell-level classification (detected / MNAR / MCAR)                ##
    ##--------------------------------------------------------------------##
    groups = unique(as.character(meta[, group.column]))

    cell.class = .missingness.cell.classes(mat = mat.log2,
                                           meta = meta,
                                           group.column = group.column,
                                           percentage.missing = percentage.missing)

    missing.matrix = is.na(mat.log2)

    ## intensity threshold used by the randomization to sample the bottom values
    tail.threshold = as.numeric(stats::quantile(x = as.vector(mat.log2), probs = tail.percentage/100, na.rm = TRUE))


    ##--------------------------------------------------------------------##
    ##  Per-protein statistics                                            ##
    ##--------------------------------------------------------------------##
    protein.stats = .missingness.protein.stats(mat.log2 = mat.log2,
                                               cell.class = cell.class,
                                               meta = meta,
                                               group.column = group.column,
                                               percentage.missing = percentage.missing)


    ##--------------------------------------------------------------------##
    ##  Per-sample statistics                                             ##
    ##--------------------------------------------------------------------##
    sample.stats =
      data.frame(column.id = colnames(mat.log2),
                 group = as.character(meta[, group.column]),
                 n.proteins = nrow(mat.log2),
                 n.detected = colSums(!missing.matrix),
                 n.missing = colSums(missing.matrix),
                 n.MNAR = colSums(cell.class == "MNAR"),
                 n.MCAR = colSums(cell.class == "MCAR"),
                 median.intensity = apply(mat.log2, 2, stats::median, na.rm = TRUE),
                 row.names = NULL) %>%
      dplyr::mutate(perc.missing = 100 * n.missing / n.proteins,
                    column.id = factor(column.id, levels = colnames(mat.log2)),
                    group = factor(group, levels = groups))


    ##--------------------------------------------------------------------##
    ##  Group-level summary                                               ##
    ##--------------------------------------------------------------------##
    group.summary = .missingness.group.summary(cell.class = cell.class,
                                               missing.matrix = missing.matrix,
                                               meta = meta,
                                               group.column = group.column,
                                               groups = groups)


    ##--------------------------------------------------------------------##
    ##  Dropout model (LOD estimation) and global statistics              ##
    ##--------------------------------------------------------------------##
    dropout = .missingness.dropout.model(protein.stats = protein.stats, min.detected = min.detected)

    complete.means = protein.stats$mean.intensity[protein.stats$missing.class == "complete"]
    incomplete.means = protein.stats$mean.intensity[protein.stats$missing.class %in% c("MNAR", "MCAR")]

    shift.test = tryCatch(expr = stats::wilcox.test(x = incomplete.means, y = complete.means, alternative = "less"),
                          error = function(e) {NULL})

    pattern.summary =
      data.frame(missing.class = factor(c("complete", "MCAR", "MNAR", "all.missing"),
                                        levels = c("complete", "MCAR", "MNAR", "all.missing"))) %>%
      dplyr::mutate(n.proteins = as.numeric(table(factor(protein.stats$missing.class,
                                                         levels = c("complete", "MCAR", "MNAR", "all.missing")))),
                    percentage = 100 * n.proteins / nrow(protein.stats))

    global.stats =
      list(n.proteins = nrow(mat.log2),
           n.samples = ncol(mat.log2),
           n.values = length(mat.log2),
           n.missing = sum(missing.matrix),
           perc.missing = 100 * sum(missing.matrix) / length(mat.log2),
           n.cells.MNAR = sum(cell.class == "MNAR"),
           n.cells.MCAR = sum(cell.class == "MCAR"),
           perc.missing.MNAR = ifelse(sum(missing.matrix) > 0,
                                      yes = 100 * sum(cell.class == "MNAR") / sum(missing.matrix),
                                      no = NA),
           proteins.per.class = pattern.summary,
           LOD50 = dropout$LOD50,
           dropout.slope = dropout$slope,
           dropout.pvalue = dropout$p.value,
           intensity.shift.pvalue = ifelse(is.null(shift.test), yes = NA, no = shift.test$p.value),
           median.intensity.complete = stats::median(complete.means, na.rm = TRUE),
           median.intensity.incomplete = stats::median(incomplete.means, na.rm = TRUE),
           tail.threshold = tail.threshold)


    ##--------------------------------------------------------------------##
    ## Plots                                                              ##
    ##--------------------------------------------------------------------##
    plots.out = .missingness.plots(mat.log2 = mat.log2,
                                   cell.class = cell.class,
                                   missing.matrix = missing.matrix,
                                   protein.stats = protein.stats,
                                   sample.stats = sample.stats,
                                   group.summary = group.summary,
                                   pattern.summary = pattern.summary,
                                   dropout = dropout,
                                   meta = meta,
                                   group.column = group.column,
                                   groups = groups,
                                   tail.threshold = tail.threshold,
                                   tail.percentage = tail.percentage,
                                   percentage.missing = percentage.missing,
                                   colors = colors,
                                   jaccard.palette = jaccard.palette,
                                   dendrogram.color = dendrogram.color,
                                   max.proteins.heatmap = max.proteins.heatmap,
                                   cluster.method = cluster.method,
                                   log.base = ifelse(convert.log2 == TRUE, yes = 2, no = log.base),
                                   seed = seed)


    ##--------------------------------------------------------------------##
    ##  Contrast-level diagnostics                                        ##
    ##--------------------------------------------------------------------##
    contrast.stats = NULL

    if ("DEprot.analyses" %in% class(DEprot.object) &
        !.deprot_slot_is_empty(DEprot.object@contrasts) &
        !(length(contrasts) == 1 && (identical(contrasts, NA) | identical(tolower(as.character(contrasts)), "none")))) {

      contrast.list = DEprot.object@contrasts

      if (is.null(contrasts)) {
        contrast.ids = seq_along(contrast.list)
      } else if (is.numeric(contrasts)) {
        if (!all(contrasts %in% seq_along(contrast.list))) {
          stop(paste0("The contrasts requested are not available: only ", length(contrast.list), " contrasts are stored in this object."))
        }
        contrast.ids = contrasts
      } else {
        if (!all(contrasts %in% names(contrast.list))) {
          stop("Some of the contrast names requested are not available in the object.")
        }
        contrast.ids = match(contrasts, names(contrast.list))
      }

      contrast.stats =
        lapply(X = contrast.ids,
               FUN = function(i) {
                 .missingness.by.contrast(contrast = contrast.list[[i]],
                                          contrast.name = names(contrast.list)[i],
                                          mat.log2 = mat.log2,
                                          meta = meta,
                                          percentage.missing = percentage.missing,
                                          tail.threshold = tail.threshold,
                                          tail.percentage = tail.percentage,
                                          colors = colors)
               })

      names(contrast.stats) = names(contrast.list)[contrast.ids]
    }


    ##--------------------------------------------------------------------##
    ##  Export                                                            ##
    ##--------------------------------------------------------------------##
    parameters = list(which.data = which.data,
                      data.used = data.used,
                      counts.available = available,
                      group.column = group.column,
                      percentage.missing = percentage.missing,
                      tail.percentage = tail.percentage,
                      tail.threshold = tail.threshold,
                      parameters.source = params.source,
                      min.detected = min.detected,
                      max.proteins.heatmap = max.proteins.heatmap,
                      cluster.method = cluster.method,
                      jaccard.palette = jaccard.palette,
                      convert.log2 = convert.log2,
                      seed = seed)

    result = methods::new(Class = "DEprot.missingness",
                          data.used = data.used,
                          counts.available = available,
                          metadata = meta,
                          group.column = group.column,
                          missing.matrix = missing.matrix,
                          imputation.map = cell.class,
                          protein.stats = protein.stats,
                          sample.stats = sample.stats,
                          group.summary = group.summary,
                          pattern.summary = pattern.summary,
                          global.stats = global.stats,
                          dropout.model = dropout$model,
                          plots = plots.out$plots,
                          jaccard.matrix = plots.out$jaccard.matrix,
                          jaccard.cluster = plots.out$jaccard.cluster,
                          contrast.stats = contrast.stats,
                          parameters = parameters)

    return(result)
  } ## END function






##=======================================================================##
##                 INTERNAL FUNCTIONS for MISSINGNESS                   ##
##=======================================================================##


#' @title .missingness.cell.classes
#' @description Classifies each cell of the counts matrix as 'detected', 'MNAR' or 'MCAR' following the
#' \code{DEprot} double-imputation logic (group-wise missingness threshold).
#' @param mat Numeric matrix of counts.
#' @param meta Metadata data.frame.
#' @param group.column String indicating the grouping column.
#' @param percentage.missing Numeric threshold (percentage).
#' @return A character matrix with the same dimensions of \code{mat}.
#' @keywords internal

.missingness.cell.classes =
  function(mat,
           meta,
           group.column,
           percentage.missing) {

    cell.class = matrix(data = "detected", nrow = nrow(mat), ncol = ncol(mat),
                        dimnames = dimnames(mat))

    groups = unique(as.character(meta[, group.column]))

    for (g in groups) {
      samples.in.group = meta$column.id[as.character(meta[, group.column]) == g]
      idx = match(samples.in.group, colnames(mat))
      sub.na = is.na(mat[, idx, drop = FALSE])
      n.NA = rowSums(sub.na)

      ## same rule as randomize.missing.values(), but requiring at least one NA
      threshold = max(1, floor((percentage.missing/100) * length(idx)))
      is.MNAR.group = (n.NA >= threshold)

      class.sub = ifelse(sub.na, yes = "MCAR", no = "detected")
      class.sub[is.MNAR.group & sub.na] = "MNAR"
      cell.class[, idx] = class.sub
    }

    return(cell.class)
  }




#' @title .missingness.protein.stats
#' @description Builds the per-protein missingness table.
#' @param mat.log2 Numeric matrix of log2 counts.
#' @param cell.class Character matrix from \link{.missingness.cell.classes}.
#' @param meta Metadata data.frame.
#' @param group.column String indicating the grouping column.
#' @param percentage.missing Numeric threshold (percentage).
#' @return A data.frame.
#' @keywords internal

.missingness.protein.stats =
  function(mat.log2,
           cell.class,
           meta,
           group.column,
           percentage.missing) {

    missing.matrix = is.na(mat.log2)
    groups = unique(as.character(meta[, group.column]))

    ## per-group number of missing values and MNAR status
    group.NA = matrix(data = NA_real_, nrow = nrow(mat.log2), ncol = length(groups),
                      dimnames = list(rownames(mat.log2), groups))
    group.MNAR = group.NA

    for (g in groups) {
      idx = match(meta$column.id[as.character(meta[, group.column]) == g], colnames(mat.log2))
      group.NA[, g] = rowSums(missing.matrix[, idx, drop = FALSE])
      threshold = max(1, floor((percentage.missing/100) * length(idx)))
      group.MNAR[, g] = as.numeric(group.NA[, g] >= threshold)
    }

    n.missing = rowSums(missing.matrix)
    n.samples = ncol(mat.log2)
    n.MNAR.groups = rowSums(group.MNAR)

    missing.class =
      dplyr::case_when(n.missing == 0 ~ "complete",
                       n.missing == n.samples ~ "all.missing",
                       n.MNAR.groups > 0 ~ "MNAR",
                       TRUE ~ "MCAR")

    protein.stats =
      data.frame(prot.id = rownames(mat.log2),
                 n.samples = n.samples,
                 n.detected = n.samples - n.missing,
                 n.missing = n.missing,
                 freq.missing = n.missing / n.samples,
                 n.cells.MNAR = rowSums(cell.class == "MNAR"),
                 n.cells.MCAR = rowSums(cell.class == "MCAR"),
                 n.groups.MNAR = n.MNAR.groups,
                 groups.MNAR = apply(group.MNAR, 1, function(x) {paste(groups[x == 1], collapse = ",")}),
                 mean.intensity = rowMeans(mat.log2, na.rm = TRUE),
                 sd.intensity = apply(mat.log2, 1, stats::sd, na.rm = TRUE),
                 min.intensity = suppressWarnings(apply(mat.log2, 1, min, na.rm = TRUE)),
                 max.intensity = suppressWarnings(apply(mat.log2, 1, max, na.rm = TRUE)),
                 missing.class = factor(missing.class, levels = c("complete", "MCAR", "MNAR", "all.missing")),
                 row.names = NULL,
                 stringsAsFactors = FALSE)

    ## non-finite values generated by fully missing proteins
    protein.stats$mean.intensity[!is.finite(protein.stats$mean.intensity)] = NA
    protein.stats$min.intensity[!is.finite(protein.stats$min.intensity)] = NA
    protein.stats$max.intensity[!is.finite(protein.stats$max.intensity)] = NA

    ## append the per-group counts of missing values
    colnames(group.NA) = paste0("n.missing_", colnames(group.NA))
    protein.stats = cbind(protein.stats, as.data.frame(group.NA, row.names = NULL))

    return(protein.stats)
  }




#' @title .missingness.group.summary
#' @description Summarizes the missing values per group of replicates.
#' @param cell.class Character matrix from \link{.missingness.cell.classes}.
#' @param missing.matrix Logical matrix of missing values.
#' @param meta Metadata data.frame.
#' @param group.column String indicating the grouping column.
#' @param groups Vector of group names.
#' @return A data.frame.
#' @keywords internal

.missingness.group.summary =
  function(cell.class,
           missing.matrix,
           meta,
           group.column,
           groups) {

    summary.list =
      lapply(X = groups,
             FUN = function(g) {
               idx = match(meta$column.id[as.character(meta[, group.column]) == g], colnames(missing.matrix))
               sub.na = missing.matrix[, idx, drop = FALSE]
               sub.class = cell.class[, idx, drop = FALSE]
               n.NA = rowSums(sub.na)

               data.frame(group = g,
                          n.samples = length(idx),
                          n.proteins = nrow(sub.na),
                          n.values = length(sub.na),
                          n.missing = sum(sub.na),
                          perc.missing = 100 * sum(sub.na) / length(sub.na),
                          n.cells.MNAR = sum(sub.class == "MNAR"),
                          n.cells.MCAR = sum(sub.class == "MCAR"),
                          proteins.complete = sum(n.NA == 0),
                          proteins.MNAR = sum(rowSums(sub.class == "MNAR") > 0),
                          proteins.MCAR = sum(rowSums(sub.class == "MNAR") == 0 & n.NA > 0),
                          stringsAsFactors = FALSE)
             })

    summary.tb = do.call(rbind, summary.list)
    summary.tb$group = factor(summary.tb$group, levels = groups)

    return(summary.tb)
  }




#' @title .missingness.dropout.model
#' @description Fits a logistic dropout model (fraction of missing values ~ average intensity) and estimates the LOD50.
#' @param protein.stats Data.frame generated by \link{.missingness.protein.stats}.
#' @param min.detected Minimal number of samples in which a protein must be detected.
#' @return A list with the model, the LOD50, the slope, the p-value and the prediction table.
#' @keywords internal

.missingness.dropout.model =
  function(protein.stats,
           min.detected = 1) {

    tb = protein.stats[protein.stats$n.detected >= min.detected & !is.na(protein.stats$mean.intensity), , drop = FALSE]

    model = tryCatch(expr = stats::glm(formula = cbind(n.missing, n.detected) ~ mean.intensity,
                                       family = stats::binomial(link = "logit"),
                                       data = tb),
                     error = function(e) {NULL},
                     warning = function(w) {
                       suppressWarnings(stats::glm(formula = cbind(n.missing, n.detected) ~ mean.intensity,
                                                   family = stats::binomial(link = "logit"),
                                                   data = tb))})

    if (is.null(model)) {
      return(list(model = NULL, LOD50 = NA, slope = NA, p.value = NA, prediction = NULL, data = tb))
    }

    coefs = stats::coef(model)
    slope = unname(coefs[2])
    LOD50 = ifelse(is.na(slope) | slope == 0, yes = NA, no = unname(-coefs[1]/coefs[2]))
    p.value = tryCatch(expr = summary(model)$coefficients[2, 4], error = function(e) {NA})

    ## prediction over the observed intensity range
    x.range = range(tb$mean.intensity, na.rm = TRUE)
    new.data = data.frame(mean.intensity = seq(from = x.range[1], to = x.range[2], length.out = 200))
    new.data$freq.missing = stats::predict(object = model, newdata = new.data, type = "response")

    ## LOD50 is meaningful only if it falls within the observed range
    if (!is.na(LOD50)) {
      if (LOD50 < x.range[1] | LOD50 > x.range[2]) {LOD50 = NA}
    }

    return(list(model = model,
                LOD50 = LOD50,
                slope = slope,
                p.value = p.value,
                prediction = new.data,
                data = tb))
  }




#' @title .missingness.render.check
#' @description Verifies that a plot can effectively be drawn. Some plots (in particular the ones built by
#' \code{ComplexUpset}) are evaluated only at printing time, therefore a package-version incompatibility would
#' raise an error when the user prints the object rather than when the object is generated. The plot is drawn on
#' a null device and, if the rendering fails, \code{NULL} is returned instead of a broken plot.
#' @param plot A ggplot/patchwork object (or \code{NULL}).
#' @param plot.name String used in the warning message. Default: \code{"plot"}.
#' @return The input plot if it can be rendered, \code{NULL} otherwise.
#' @importFrom grDevices pdf dev.off
#' @keywords internal

.missingness.render.check =
  function(plot,
           plot.name = "plot") {

    if (is.null(plot)) {return(NULL)}

    tmp.file = tempfile(fileext = ".pdf")

    rendered =
      tryCatch(expr = {
        grDevices::pdf(file = tmp.file)
        on.exit(expr = {grDevices::dev.off(); unlink(tmp.file)}, add = TRUE)
        print(plot)
        TRUE},
        error = function(e) {
          warning(paste0("The ", plot.name, " could not be rendered and was not included in the output ",
                         "(", conditionMessage(e), ")."), call. = FALSE)
          return(FALSE)})

    if (isTRUE(rendered)) {return(plot)} else {return(NULL)}
  }




#' @title .missingness.plots
#' @description Generates all the diagnostic plots.
#' @param mat.log2 Numeric matrix of log2 counts.
#' @param cell.class Character matrix of the cell classes.
#' @param missing.matrix Logical matrix of missing values.
#' @param protein.stats,sample.stats,group.summary,pattern.summary Tables generated by the internal functions.
#' @param dropout Output of \link{.missingness.dropout.model}.
#' @param meta Metadata data.frame.
#' @param group.column,groups Grouping column and group names.
#' @param tail.threshold,tail.percentage,percentage.missing Randomization parameters.
#' @param colors Named vector of colors.
#' @param jaccard.palette Vector of colors for the scale of the Jaccard similarity heatmap.
#' @param dendrogram.color String indicating the color of the dendrogram lines.
#' @param max.proteins.heatmap,cluster.method,log.base,seed Plotting parameters.
#'
#' @importFrom ggpubr theme_pubr
#'
#' @return A list containing the ggplot objects (\code{plots}), the Jaccard similarity matrix (\code{jaccard.matrix})
#' and the corresponding \code{hclust} object (\code{jaccard.cluster}).
#'
#' @keywords internal

.missingness.plots =
  function(mat.log2,
           cell.class,
           missing.matrix,
           protein.stats,
           sample.stats,
           group.summary,
           pattern.summary,
           dropout,
           meta,
           group.column,
           groups,
           tail.threshold,
           tail.percentage,
           percentage.missing,
           colors,
           jaccard.palette,
           dendrogram.color,
           max.proteins.heatmap,
           cluster.method,
           log.base,
           seed) {

    intensity.lab = ifelse(is.na(log.base),
                           yes = "Mean intensity",
                           no = paste0("Mean ", ifelse(log.base == exp(1), yes = "ln", no = paste0("log<sub>", log.base, "</sub>")), "(Intensity)"))

    base.theme =
      ggpubr::theme_pubr(legend = "right") +
      theme(axis.text = element_text(color = "black"),
            axis.title.x = ggtext::element_markdown(color = "black"),
            axis.title.y = ggtext::element_markdown(color = "black"),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5),
            legend.title = ggtext::element_markdown(),
            strip.background = element_blank())


    ##--------------------------------------------------------##
    ## A. Detection density (MCAR vs MNAR key diagnostic)     ##
    ##--------------------------------------------------------##
    density.tb =
      protein.stats %>%
      dplyr::filter(!is.na(mean.intensity)) %>%
      dplyr::mutate(completeness = factor(ifelse(missing.class == "complete",
                                                 yes = "fully quantified",
                                                 no = "with missing value(s)"),
                                          levels = c("fully quantified", "with missing value(s)")))

    detection.density =
      ggplot(data = density.tb,
             mapping = aes(x = mean.intensity, fill = completeness)) +
      geom_density(alpha = 0.4, color = NA) +
      geom_vline(xintercept = tail.threshold, linetype = 2, color = "black") +
      annotate(geom = "text", x = tail.threshold, y = Inf, hjust = -0.05, vjust = 1.5, size = 3,
               label = paste0("randomization pool: \u2264 ", round(tail.threshold, 2),
                              " (bottom ", tail.percentage, "%)")) +
      scale_fill_manual(values = c("fully quantified" = "grey50", "with missing value(s)" = unname(colors["MNAR"])), name = NULL) +
      xlab(intensity.lab) +
      ylab("Density") +
      ggtitle("**Detection *vs* abundance**",
              subtitle = "*superimposed: MCAR-like &nbsp;&nbsp;|&nbsp;&nbsp; left-shifted: MNAR-like*") +
      base.theme +
      theme(aspect.ratio = 0.7)


    ##--------------------------------------------------------##
    ## B. Dropout curve (missingness vs intensity)            ##
    ##--------------------------------------------------------##
    dropout.curve =
      ggplot(data = dropout$data,
             mapping = aes(x = mean.intensity, y = 100 * freq.missing)) +
      geom_point(mapping = aes(color = missing.class), alpha = 0.35, stroke = NA, size = 1.5) +
      scale_color_manual(values = stats::setNames(unname(colors[c("detected", "MCAR", "MNAR", "all.missing")]),
                                                  c("complete", "MCAR", "MNAR", "all.missing")),
                         name = "Class",
                         drop = FALSE) +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      xlab(intensity.lab) +
      ylab("Missing values (%)") +
      ggtitle("**Dropout curve**") +
      base.theme +
      theme(aspect.ratio = 0.7)

    if (!is.null(dropout$prediction)) {
      dropout.curve =
        dropout.curve +
        geom_line(data = dropout$prediction,
                  mapping = aes(x = mean.intensity, y = 100 * freq.missing),
                  color = "black", linewidth = 0.8, inherit.aes = FALSE)
    }

    if (!is.na(dropout$LOD50)) {
      dropout.curve =
        dropout.curve +
        geom_vline(xintercept = dropout$LOD50, linetype = 2, color = "black") +
        labs(subtitle = paste0("*LOD<sub>50</sub> = ", round(dropout$LOD50, 2), "*"))
    }


    ##--------------------------------------------------------##
    ## C. Missingness heatmap (imputation map)                ##
    ##--------------------------------------------------------##
    prot.with.NA = rownames(mat.log2)[rowSums(missing.matrix) > 0]

    if (length(prot.with.NA) > 1) {
      set.seed(seed)
      if (length(prot.with.NA) > max.proteins.heatmap) {
        prot.with.NA = sample(x = prot.with.NA, size = max.proteins.heatmap, replace = FALSE)
      }

      bin.mat = matrix(as.numeric(!missing.matrix[prot.with.NA, , drop = FALSE]),
                       nrow = length(prot.with.NA),
                       dimnames = list(prot.with.NA, colnames(missing.matrix)))

      row.order =
        tryCatch(expr = {hc = stats::hclust(stats::dist(bin.mat, method = "binary"), method = cluster.method)
        prot.with.NA[hc$order]},
        error = function(e) {prot.with.NA[order(rowSums(bin.mat))]})

      heat.tb =
        suppressMessages(reshape2::melt(cell.class[row.order, , drop = FALSE])) %>%
        dplyr::rename(prot.id = "Var1", column.id = "Var2", status = "value") %>%
        dplyr::mutate(prot.id = factor(prot.id, levels = row.order),
                      column.id = factor(as.character(column.id), levels = colnames(mat.log2)),
                      status = factor(status, levels = c("detected", "MCAR", "MNAR")),
                      group = factor(meta[match(as.character(column.id), meta$column.id), group.column], levels = groups))

      missingness.heatmap =
        ggplot(data = heat.tb,
               mapping = aes(x = column.id, y = prot.id, fill = status)) +
        geom_tile() +
        facet_grid(cols = vars(group), scales = "free_x", space = "free_x") +
        scale_fill_manual(values = colors[c("detected", "MCAR", "MNAR")], name = NULL, drop = FALSE,
                          labels = c(detected = "detected",
                                     MCAR = "missing: MCAR-like (imputed)",
                                     MNAR = "missing: MNAR-like (randomized)")) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        xlab(NULL) +
        ylab(paste0("Proteins with \u22651 missing value (n = ", length(row.order), ")")) +
        ggtitle("**Missingness pattern**") +
        base.theme +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_blank(),
              panel.spacing = unit(3, "pt"),
              panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
        )
    } else {
      missingness.heatmap = NULL
    }


    ##--------------------------------------------------------##
    ## D. Missing values per sample                           ##
    ##--------------------------------------------------------##
    missing.per.sample =
      ggplot(data = sample.stats,
             mapping = aes(x = column.id, y = perc.missing, fill = group)) +
      geom_col(color = NA, width = 0.8) +
      geom_hline(yintercept = mean(sample.stats$perc.missing), linetype = 2, color = "grey30") +
      scale_fill_manual(values = stats::setNames(viridis::viridis(n = length(groups), option = "D", end = 0.9), groups), name = group.column) +
      xlab(NULL) +
      ylab("Missing values (%)") +
      ggtitle("**Missingness per sample**") +
      base.theme +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            axis.ticks.x = element_blank(),
            aspect.ratio = 10/nrow(sample.stats))


    ##--------------------------------------------------------##
    ## E. Detection frequency per protein                     ##
    ##--------------------------------------------------------##
    detection.frequency =
      ggplot(data = protein.stats,
             mapping = aes(x = n.detected)) +
      geom_bar(fill = "grey40", color = NA, width = 0.8) +
      scale_x_continuous(breaks = 0:ncol(mat.log2)) +
      xlab("Number of samples in which a protein is detected") +
      ylab("Number of proteins") +
      ggtitle("**Detection frequency**") +
      base.theme +
      theme(aspect.ratio = 0.6)


    ##--------------------------------------------------------##
    ## F. Pattern summary (global + per group)                ##
    ##--------------------------------------------------------##
    global.bar =
      pattern.summary %>%
      dplyr::mutate(group = "all samples") %>%
      dplyr::select(group, missing.class, n.proteins)

    group.bar =
      group.summary %>%
      dplyr::select(group, proteins.complete, proteins.MCAR, proteins.MNAR) %>%
      dplyr::rename(complete = "proteins.complete", MCAR = "proteins.MCAR", MNAR = "proteins.MNAR") %>%
      reshape2::melt(id.vars = "group", variable.name = "missing.class", value.name = "n.proteins") %>%
      dplyr::mutate(group = as.character(group))

    bar.tb =
      rbind(global.bar, group.bar) %>%
      dplyr::mutate(group = factor(group, levels = c("all samples", groups)),
                    missing.class = factor(as.character(missing.class),
                                           levels = c("complete", "MCAR", "MNAR", "all.missing")))

    pattern.barplot =
      ggplot(data = bar.tb,
             mapping = aes(x = group, y = n.proteins, fill = missing.class)) +
      geom_col(color = NA, width = 0.8) +
      scale_fill_manual(values = stats::setNames(unname(colors[c("detected", "MCAR", "MNAR", "all.missing")]),
                                                 c("complete", "MCAR", "MNAR", "all.missing")),
                        name = NULL, drop = FALSE) +
      xlab(NULL) +
      ylab("Number of proteins") +
      ggtitle("**Missing-value classes**",
              subtitle = paste0("*MNAR: \u2265", percentage.missing, "% missing within a group of '", group.column, "'*")) +
      base.theme +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            axis.ticks.x = element_blank(),
            aspect.ratio = 0.8)


    ##--------------------------------------------------------##
    ## G. Sample-sample Jaccard similarity of detection       ##
    ##--------------------------------------------------------##
    detected = !missing.matrix
    intersection = t(detected) %*% detected
    total = matrix(colSums(detected), nrow = ncol(detected), ncol = ncol(detected))
    union.mat = total + t(total) - intersection
    jaccard = intersection / union.mat
    dimnames(jaccard) = list(colnames(mat.log2), colnames(mat.log2))

    ## Clustering of the samples on the Jaccard dissimilarity
    jaccard.distance = stats::as.dist(1 - jaccard)
    jaccard.cluster = stats::hclust(d = jaccard.distance, method = cluster.method)
    sample.order = colnames(jaccard)[jaccard.cluster$order]

    jaccard.tb =
      suppressMessages(reshape2::melt(jaccard)) %>%
      dplyr::rename(sample.x = "Var1", sample.y = "Var2", value = "value") %>%
      dplyr::mutate(sample.x = factor(as.character(sample.x), levels = rev(sample.order)),
                    sample.y = factor(as.character(sample.y), levels = rev(sample.order)))

    sample.similarity =
      ggplot(data = jaccard.tb,
             mapping = aes(x = sample.x, y = sample.y, fill = value)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradientn(colors = jaccard.palette, name = "Jaccard<br>(detection)") +
      scale_x_discrete(expand = c(0, 0)) +
      legendry::scale_y_dendro(clust = jaccard.cluster,
                               expand = c(0, 0),
                               position = "left") +
      xlab(NULL) + ylab(NULL) +
      ggtitle("**Similarity of the detection patterns**",
              subtitle = "*a group/batch structure indicates non-random missingness*") +
      base.theme +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = dendrogram.color, linewidth = 0.5),
            panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
            aspect.ratio = 1)

    ## keep track of how the dendrogram was built (as for DEprot.correlation objects)
    jaccard.cluster$dist.method = "as.dist(1 - jaccard.matrix)"
    jaccard.cluster$call = paste0("hclust(d = as.dist(1 - jaccard.matrix), method = ", cluster.method, ")")


    ##--------------------------------------------------------##
    ## H. UpSet of the MNAR-like patterns between groups      ##
    ##--------------------------------------------------------##
    upset.plot =
      tryCatch(expr = {
        obs = as.data.frame(sapply(groups, function(g) {grepl(paste0("(^|,)", g, "(,|$)"), protein.stats$groups.MNAR)}))
        colnames(obs) = groups
        obs = obs[rowSums(obs) > 0, , drop = FALSE]

        if (nrow(obs) > 0 & length(groups) > 1) {

          ## ComplexUpset expects real 'theme' objects: a list of theme components
          ## (e.g. list(theme_minimal(), theme(...))) is not recognized as a valid
          ## theme since ggplot2 4.0.0. All the panel keys are provided explicitly,
          ## so that the (stale) internal defaults of ComplexUpset are never used.
          upset.base.theme =
            theme_minimal() +
            theme(axis.text.x.bottom = ggtext::element_markdown(color = "black"),
                  axis.text.y.left = ggtext::element_markdown(color = "black"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  panel.grid.major.x = element_blank())

          upset.themes =
            list("intersections_matrix" = upset.base.theme,
                 "overall_sizes" = upset.base.theme + theme(panel.grid.major.y = element_blank()),
                 "Intersection size" = upset.base.theme + theme(axis.text.x = element_blank()),
                 "default" = upset.base.theme + theme(axis.text.x = element_blank()))

          ComplexUpset::upset(data = obs,
                              intersect = groups,
                              name = NULL,
                              base_annotations =
                                list("Intersection size" =
                                       ComplexUpset::intersection_size(fill = unname(colors["MNAR"]),
                                                                       counts = TRUE,
                                                                       width = 0.75) +
                                       scale_y_continuous(expand = c(0, 0)) +
                                       ggtitle(label = "**MNAR-like proteins shared between groups**",
                                               subtitle = paste0("*missing in \u2265", percentage.missing, "% of the replicates of a group*")) +
                                       theme(axis.text.x.bottom = ggtext::element_markdown(color = "black"),
                                             axis.text.y.left = ggtext::element_markdown(color = "black"),
                                             panel.grid.major.x = element_blank(),
                                             axis.line.y = element_line(colour = "black"),
                                             axis.ticks.y = element_line(colour = "black"),
                                             plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
                                             plot.subtitle = ggtext::element_markdown(color = "black", hjust = 0.5))),
                              set_sizes =
                                ComplexUpset::upset_set_size(geom = geom_bar(fill = "black", width = 0.4, show.legend = FALSE)) +
                                scale_y_reverse(expand = c(0, 0)) +
                                theme(axis.line.x = element_line(colour = "black"),
                                      axis.ticks.x = element_line(colour = "black")),
                              themes = upset.themes,
                              sort_intersections_by = "cardinality",
                              sort_sets = FALSE,
                              min_size = 1) }
        else {NULL}
      },
      error = function(e) {NULL})

    ## ComplexUpset objects are evaluated only when printed: check the rendering
    upset.plot = .missingness.render.check(plot = upset.plot,
                                           plot.name = "UpSet plot of the MNAR-like patterns")


    ##--------------------------------------------------------##
    return(list(plots = list(detection.density = detection.density,
                             dropout.curve = dropout.curve,
                             missingness.heatmap = missingness.heatmap,
                             missing.per.sample = missing.per.sample,
                             detection.frequency = detection.frequency,
                             pattern.barplot = pattern.barplot,
                             sample.similarity = sample.similarity,
                             upset = upset.plot),
                jaccard.matrix = jaccard,
                jaccard.cluster = jaccard.cluster))
  }




#' @title .missingness.by.contrast
#' @description Computes the missingness diagnostics restricted to the samples of a specific contrast.
#' @param contrast List describing the contrast (as stored in \code{DEprot.analyses@@contrasts}).
#' @param contrast.name String indicating the name of the contrast.
#' @param mat.log2 Numeric matrix of log2 counts.
#' @param meta Metadata data.frame.
#' @param percentage.missing Numeric threshold (percentage).
#' @param tail.threshold Numeric value of the bottom-distribution threshold (intensity).
#' @param tail.percentage Numeric value of the bottom percentage of the distribution used by the randomization.
#' @param colors Named vector of colors.
#'
#' @importFrom ggpubr theme_pubr
#'
#' @return A list with the tables and the plots of the contrast.
#'
#' @keywords internal

.missingness.by.contrast =
  function(contrast,
           contrast.name,
           mat.log2,
           meta,
           percentage.missing,
           tail.threshold,
           tail.percentage,
           colors) {

    metadata.column = contrast$metadata.column
    var.1 = contrast$var.1
    var.2 = contrast$var.2

    group.1 = if (!is.null(contrast$group.1)) {contrast$group.1} else {meta$column.id[as.character(meta[, metadata.column]) == var.1]}
    group.2 = if (!is.null(contrast$group.2)) {contrast$group.2} else {meta$column.id[as.character(meta[, metadata.column]) == var.2]}

    group.1 = intersect(group.1, colnames(mat.log2))
    group.2 = intersect(group.2, colnames(mat.log2))

    if (length(group.1) == 0 | length(group.2) == 0) {return(NULL)}

    sub.mat = mat.log2[, c(group.1, group.2), drop = FALSE]
    sub.meta = meta[match(c(group.1, group.2), meta$column.id), , drop = FALSE]

    NA.1 = rowSums(is.na(mat.log2[, group.1, drop = FALSE]))
    NA.2 = rowSums(is.na(mat.log2[, group.2, drop = FALSE]))

    th.1 = max(1, floor((percentage.missing/100) * length(group.1)))
    th.2 = max(1, floor((percentage.missing/100) * length(group.2)))

    MNAR.1 = NA.1 >= th.1
    MNAR.2 = NA.2 >= th.2
    n.NA = NA.1 + NA.2

    contrast.class =
      dplyr::case_when(n.NA == 0 ~ "complete",
                       MNAR.1 & MNAR.2 ~ "all.missing",
                       MNAR.1 & !MNAR.2 ~ paste0("MNAR in ", var.1),
                       MNAR.2 & !MNAR.1 ~ paste0("MNAR in ", var.2),
                       TRUE ~ "MCAR")

    class.levels = c("complete", "MCAR", paste0("MNAR in ", var.1), paste0("MNAR in ", var.2), "all.missing")

    protein.stats =
      data.frame(prot.id = rownames(mat.log2),
                 n.missing.group1 = NA.1,
                 n.missing.group2 = NA.2,
                 n.samples.group1 = length(group.1),
                 n.samples.group2 = length(group.2),
                 mean.intensity = rowMeans(sub.mat, na.rm = TRUE),
                 mean.intensity.group1 = rowMeans(mat.log2[, group.1, drop = FALSE], na.rm = TRUE),
                 mean.intensity.group2 = rowMeans(mat.log2[, group.2, drop = FALSE], na.rm = TRUE),
                 missing.class = factor(contrast.class, levels = class.levels),
                 testable = !(contrast.class == "all.missing"),
                 row.names = NULL,
                 stringsAsFactors = FALSE)

    protein.stats$mean.intensity[!is.finite(protein.stats$mean.intensity)] = NA
    protein.stats$mean.intensity.group1[!is.finite(protein.stats$mean.intensity.group1)] = NA
    protein.stats$mean.intensity.group2[!is.finite(protein.stats$mean.intensity.group2)] = NA

    summary.tb =
      data.frame(missing.class = factor(class.levels, levels = class.levels)) %>%
      dplyr::mutate(n.proteins = as.numeric(table(factor(protein.stats$missing.class, levels = class.levels))),
                    percentage = 100 * n.proteins / nrow(protein.stats))

    class.colors = stats::setNames(c(unname(colors["detected"]), unname(colors["MCAR"]),
                                     unname(colors["MNAR"]), colorspace::lighten(unname(colors["MNAR"]), 0.4),
                                     unname(colors["all.missing"])),
                                   class.levels)

    contrast.title = paste0("**", var.1, "** *vs* **", var.2, "**")

    barplot =
      ggplot(data = summary.tb,
             mapping = aes(x = missing.class, y = n.proteins, fill = missing.class)) +
      geom_col(color = NA, width = 0.8, show.legend = FALSE) +
      geom_text(mapping = aes(label = n.proteins), vjust = -0.4, size = 3) +
      scale_fill_manual(values = class.colors, drop = FALSE) +
      xlab(NULL) +
      ylab("Number of proteins") +
      ggtitle(contrast.title,
              subtitle = paste0("*missing-value classes (MNAR: \u2265", percentage.missing, "% per group)*")) +
      ggpubr::theme_pubr(legend = "right") +
      theme(axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 30, hjust = 1),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5),
            axis.ticks.x = element_blank(),
            aspect.ratio = 0.8)

    density.plot =
      ggplot(data = dplyr::filter(protein.stats, !is.na(mean.intensity), missing.class != "all.missing"),
             mapping = aes(x = mean.intensity,
                           fill = factor(ifelse(missing.class == "complete", "fully quantified", "with missing value(s)"),
                                         levels = c("fully quantified", "with missing value(s)")))) +
      geom_density(alpha = 0.4, linewidth = 0.4, color = NA) +
      geom_vline(xintercept = tail.threshold, linetype = 2, color = "black") +
      annotate(geom = "text", x = tail.threshold, y = Inf, hjust = -0.05, vjust = 1.5, size = 3,
               label = paste0("randomization pool: \u2264 ", round(tail.threshold, 2),
                              " (bottom ", tail.percentage, "%)")) +
      coord_cartesian(clip = "off") +
      scale_fill_manual(values = c("fully quantified" = "grey50", "with missing value(s)" = unname(colors["MNAR"])), name = NULL) +
      xlab("Mean log<sub>2</sub>(Intensity)") +
      ylab("Density") +
      ggtitle(contrast.title, subtitle = "*detection vs abundance*") +
      ggpubr::theme_pubr(legend = "right") +
      theme(axis.text = element_text(color = "black"),
            axis.title.x = ggtext::element_markdown(color = "black"),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5),
            aspect.ratio = 0.7)

    return(list(contrast.id = contrast.name,
                metadata.column = metadata.column,
                var.1 = var.1,
                var.2 = var.2,
                samples = list(group.1 = group.1, group.2 = group.2),
                protein.stats = protein.stats,
                summary = summary.tb,
                plots = list(pattern.barplot = barplot,
                             detection.density = density.plot)))
  }
