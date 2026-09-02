##########################################################################################
###                           DEprot :: TIME-COURSE MODULE                             ###
##########################################################################################
#
#   analyze.timecourse()          spline-based trend testing + soft clustering
#   rank.timecourse()             re-ranking without refitting
#   get.timecourse.results()      results getter, with cluster subsetting
#   plot.timecourse.protein()     trend of individual proteins
#   plot.timecourse.profiles()    cluster profile panels
#   heatmap.timecourse()          heatmap of the trending proteins, split by cluster
#   timecourse.enrichment()       over-representation analyses (ORA) cluster by cluster
#   .simulate.timecourse()        timecourse data simulator
#
#   Time is always handled as a NUMERIC covariate.
#
##########################################################################################



# ----------------------------------------------------------------------------------------
###                                INTERNAL FUNCTIONS                                   ###
# ----------------------------------------------------------------------------------------

#' @title .tc.get.counts
#'
#' @description
#' Extracts the requested counts table from a DEprot object, accepting the same aliases used elsewhere in the package.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'imp', 'randomized', 'random'.
#'
#' @return A list with two elements: \code{mat} (the counts matrix) and \code{data.used} (the normalized name of the table used).
#'
#' @keywords internal

.tc.get.counts =
  function(DEprot.object,
           which.data) {

    ## select the slot corresponding to the requested counts (aliases allowed)
    if (tolower(which.data) == "raw") {
      mat = DEprot.object@raw.counts
      data.used = "raw"
    } else if (tolower(which.data) %in% c("norm", "normalized", "normal")) {
      mat = DEprot.object@norm.counts
      data.used = "normalized"
    } else if (tolower(which.data) %in% c("imputed", "imp", "impute")) {
      mat = DEprot.object@imputed.counts
      data.used = "imputed"
    } else if (tolower(which.data) %in% c("randomized", "random")) {
      mat = DEprot.object@random.counts
      data.used = "randomized"
    } else {
      stop(paste0("The 'which.data' value is not recognized.\n",
                  "       Please indicate a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
    }

    ## the slot exists but has never been filled
    if (is.null(mat)) {
      stop(paste0("Use of ", toupper(data.used), " counts was required, but not available.\n",
                  "       Please indicate a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
    }

    return(list(mat = as.matrix(mat),
                data.used = data.used))
  } # END function






#' @title .tc.transform.time
#'
#' @description
#' Applies a transformation to the numeric time vector. The log transformations are computed as log(t + 1) so that time zero is allowed.
#'
#' @param time Numeric vector of the timepoints.
#' @param time.transform String indicating the transformation. One among: 'none', 'log2', 'log10', 'log1p', 'log', 'ln', 'sqrt'.
#'
#' @return A numeric vector of the transformed timepoints.
#'
#' @keywords internal

.tc.transform.time =
  function(time,
           time.transform) {

    tt = switch(tolower(time.transform),
                "none"  = time,
                "log2"  = log2(time + 1),
                "log10" = log10(time + 1),
                "log1p" = log1p(time),
                "log"   = log1p(time),
                "ln"    = log1p(time),
                "sqrt"  = sqrt(time),
                stop("'time.transform' must be one of: 'none', 'log2', 'log10', 'log1p', 'sqrt'."))

    return(tt)
  } # END function






#' @title .tc.backtransform.time
#'
#' @description
#' Inverse of \code{.tc.transform.time}: brings the (transformed) prediction grid back to the original time scale, so that the fitted curves can be plotted against the real timepoints.
#'
#' @param tt Numeric vector of transformed time values.
#' @param time.transform String indicating the transformation that was applied.
#'
#' @return A numeric vector on the original time scale.
#'
#' @keywords internal

.tc.backtransform.time =
  function(tt,
           time.transform) {

    ## the tolower() must match the one of .tc.transform.time, otherwise a capitalized
    ## value would be transformed but not brought back, and the switch would return NULL
    time = switch(tolower(time.transform),
                  "none"  = tt,
                  "log2"  = 2^tt - 1,
                  "log10" = 10^tt - 1,
                  "log1p" = exp(tt) - 1,
                  "log"   = exp(tt) - 1,
                  "ln"    = exp(tt) - 1,
                  "sqrt"  = tt^2,
                  stop("'time.transform' must be one of: 'none', 'log2', 'log10', 'log1p', 'sqrt'."))

    return(time)
  } # END function






#' @title .tc.mestimate
#'
#' @description
#' Estimates the fuzzifier \code{m} of the c-means algorithm from the size and the dimensionality of the data. The classical default (m = 2) makes the memberships of high-dimensional data converge towards 1/k, collapsing all the clusters onto the global centroid; this empirical relation keeps the partition informative.
#'
#' @param n Numeric value indicating the number of objects (proteins) clustered.
#' @param d Numeric value indicating the number of dimensions (timepoints).
#'
#' @return A numeric value, the estimated fuzzifier.
#'
#' @references \enc{Schwämmle}{Schwammle} V. & Jensen O.N. (2010). A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. \emph{Bioinformatics}, \strong{26}(22), 2841-2848. \doi{10.1093/bioinformatics/btq534}
#'
#' @keywords internal

.tc.mestimate =
  function(n,
           d) {

    m = 1 +
      (1418 / n + 22.05) * d^(-2) +
      (12.33 / n + 0.243) * d^(-0.0406 * log(n) - 0.1134)

    return(round(m, 3))
  } # END function






#' @title .tc.zscore
#'
#' @description
#' Row-wise Z-score of a matrix. The rows with a null variance are returned as zeros instead of NaN.
#'
#' @param x A numeric matrix.
#'
#' @return A matrix of the same dimensions.
#'
#' @importFrom stats sd
#'
#' @keywords internal

.tc.zscore =
  function(x) {

    z = t(apply(x, 1,
                function(r) {
                  s = stats::sd(r, na.rm = TRUE)
                  if (is.na(s) | s == 0) {return(rep(0, length(r)))}
                  return((r - mean(r, na.rm = TRUE)) / s)
                }))

    rownames(z) = rownames(x)
    return(z)
  } # END function






#' @title .tc.descriptors
#'
#' @description
#' Computes the kinetic descriptors of each fitted trajectory: amplitude, initial slope, time to half-maximum, time of maximal change and shape class. These summarize a continuous curve into the numbers that a user actually wants to read in a results table.
#'
#' @param curves Numeric matrix (proteins x grid points) of the fitted trajectories.
#' @param grid.time Numeric vector of the grid positions, on the original time scale.
#' @param slope.window Numeric value indicating the fraction of the grid used to compute the initial slope. Default: \code{0.05}.
#' @param flat.tolerance Numeric value below which an amplitude is considered null. Default: \code{1e-3}.
#'
#' @return A data.frame with the columns: amplitude, initial.slope, t.half, peak.time, trend.shape.
#'
#' @keywords internal

.tc.descriptors =
  function(curves,
           grid.time,
           slope.window = 0.05,
           flat.tolerance = 1e-3) {

    n.grid = ncol(curves)

    ## amplitude = full range spanned by the fitted curve (log2 units)
    amplitude = apply(curves, 1, function(x) {max(x) - min(x)})

    ## initial slope = rate of change over the first fraction of the time range, that is,
    ## how fast the protein reacts (per unit of ORIGINAL time)
    i.slope = max(2, ceiling(n.grid * slope.window))
    initial.slope = (curves[,i.slope] - curves[,1]) / (grid.time[i.slope] - grid.time[1])

    ## deviation from the baseline value, used by all the descriptors below
    deviation = curves - curves[,1]

    ## peak.time = when the curve is the furthest away from its baseline value
    peak.time = grid.time[apply(abs(deviation), 1, which.max)]

    ## t.half = when the curve reaches 50% of its OWN maximal deviation, restricted to the
    ## rising portion. Contrary to the initial slope, which scales with the size of the
    ## response, this is a pure timing parameter: it can be compared between proteins of
    ## very different abundances, and it is what orders sequential events.
    t.half = vapply(seq_len(nrow(curves)),
                    function(i) {
                      if (amplitude[i] < flat.tolerance) {return(NA_real_)}

                      i.max = which.max(abs(deviation[i,]))
                      if (i.max == 1) {return(NA_real_)}

                      crossing = which(abs(deviation[i,1:i.max]) >= 0.5 * abs(deviation[i,i.max]))
                      if (length(crossing) == 0) {return(NA_real_)}

                      return(grid.time[crossing[1]])
                    },
                    FUN.VALUE = numeric(1))

    ## shape = read from the sign changes of the numerical derivative:
    ## no change -> monotone, one change -> transient, more -> complex
    d = t(apply(curves, 1, diff))
    if (is.null(dim(d))) {d = matrix(d, nrow = nrow(curves))}

    shape = vapply(seq_len(nrow(curves)),
                   function(i) {
                     if (amplitude[i] < flat.tolerance) {return("flat")}

                     s = sign(d[i,])
                     s = s[abs(d[i,]) > (flat.tolerance / n.grid)]
                     if (length(s) == 0) {return("flat")}

                     n.changes = sum(diff(s) != 0)
                     i.max = which.max(abs(deviation[i,]))
                     direction = ifelse(deviation[i,n.grid] >= 0, "up", "down")

                     if (n.changes == 0) {
                       return(paste0("monotone.", direction))
                     } else if (n.changes == 1) {
                       ## transient: the extreme is reached before the last timepoint
                       return(paste0("transient.", ifelse(deviation[i,i.max] > 0, "up", "down")))
                     } else {
                       return("complex")
                     }
                   },
                   FUN.VALUE = character(1))

    return(data.frame(amplitude = amplitude,
                      initial.slope = initial.slope,
                      t.half = t.half,
                      peak.time = peak.time,
                      trend.shape = shape,
                      stringsAsFactors = FALSE))
  } # END function







#' @title .tc.check.values
#'
#' @description
#' Normalizes the \code{values} argument of the time-course plots, accepting aliases.
#'
#' @param values String indicating the quantity to plot. One among: 'zscore', 'log2FC', 'counts' (aliases accepted).
#'
#' @return One among "zscore", "log2FC", "counts".
#'
#' @keywords internal

.tc.check.values =
  function(values) {

    v = switch(tolower(values[1]),
               "zscore"      = "zscore",
               "z-score"     = "zscore",
               "z"           = "zscore",
               "scaled"      = "zscore",
               "log2fc"      = "log2FC",
               "fc"          = "log2FC",
               "foldchange"  = "log2FC",
               "fold.change" = "log2FC",
               "counts"      = "counts",
               "count"       = "counts",
               "log2"        = "counts",
               "lfq"         = "counts",
               stop("'values' must be one of: 'zscore', 'log2FC', 'counts'."))

    return(v)
  } # END function






#' @title .tc.ref.index
#'
#' @description
#' Identifies the timepoint used as baseline for the log2(FC) computation, and its position on the prediction grid.
#'
#' @param time.grid Numeric vector of the grid positions (original time scale).
#' @param timepoints Numeric vector of the observed timepoints.
#' @param reference.time Numeric value indicating the baseline timepoint. Default: \code{NULL} (the earliest timepoint).
#'
#' @return A list with the elements \code{time} and \code{index}.
#'
#' @keywords internal

.tc.ref.index =
  function(time.grid,
           timepoints,
           reference.time = NULL) {

    ## by default the baseline is the first timepoint of the experiment
    if (is.null(reference.time)) {reference.time = min(timepoints)}

    ## if the value provided is not an actual timepoint, it is snapped to the closest one
    if (!(reference.time %in% timepoints)) {
      closest = timepoints[which.min(abs(timepoints - reference.time))]
      warning(paste0("The 'reference.time' provided (", reference.time, ") is not among the timepoints (",
                     paste(timepoints, collapse = ", "), "): ", closest, " is used instead."))
      reference.time = closest
    }

    return(list(time = reference.time,
                index = which.min(abs(time.grid - reference.time))))
  } # END function






#' @title .tc.ylab
#'
#' @description
#' Builds the y-axis label matching the quantity plotted by the time-course functions.
#'
#' @param values String, output of \code{.tc.check.values}.
#' @param ref.time Numeric value of the baseline timepoint (used only for log2FC). Default: \code{NULL}.
#' @param time.column String indicating the name of the time column. Default: \code{"time"}.
#' @param fitted Logic value indicating whether the values plotted are fitted rather than measured. Default: \code{FALSE}.
#'
#' @return A string.
#'
#' @keywords internal

.tc.ylab =
  function(values,
           ref.time = NULL,
           time.column = "time",
           fitted = FALSE) {

    ## the subscripts are rendered by ggtext::element_markdown()
    lab = switch(values,
                 "zscore" = "Z-score",
                 "log2FC" = paste0("log~2~(FC) vs ", ref.time, " ", time.column),
                 "counts" = "log~2~(intensity)")

    if (isTRUE(fitted)) {lab = paste0(lab, " - fitted")}

    return(lab)
  } # END function





#' @title .tc.descriptor.column
#'
#' @description
#' Retrieves a kinetic descriptor from the results table, whatever the design. With a single group the columns carry no suffix; with several groups every descriptor is suffixed by the group name, and the first group level is then used as reference.
#'
#' @param results The results data.frame of a \code{DEprot.timecourse} object.
#' @param descriptor String indicating the descriptor to retrieve (e.g., 't.half', 'peak.time', 'initial.slope').
#' @param group.levels Character vector of the group levels of the analysis. Default: \code{NULL}.
#'
#' @return A numeric vector.
#'
#' @keywords internal

.tc.descriptor.column =
  function(results,
           descriptor,
           group.levels = NULL) {

    ## single group: the column carries no suffix
    if (descriptor %in% colnames(results)) {return(results[[descriptor]])}

    ## several groups: the descriptors are suffixed by the group name
    candidates = colnames(results)[startsWith(colnames(results), paste0(descriptor, "."))]

    if (length(candidates) == 0) {
      stop(paste0("The descriptor '", descriptor, "' is not available in this object."))
    }

    ## the first group level is used as reference, when it can be identified
    if (!is.null(group.levels)) {
      preferred = paste0(descriptor, ".", group.levels[1])
      if (preferred %in% candidates) {candidates = preferred}
    }

    return(results[[candidates[1]]])
  } # END function





# ----------------------------------------------------------------------------------------
###                                 MAIN FUNCTION                                       ###
# ----------------------------------------------------------------------------------------

#' @title analyze.timecourse
#'
#' @description Tests each protein for a temporal trend treating time as a NUMERIC covariate. A natural-spline basis of the (optionally transformed) time is fitted by \code{limma} and all the spline coefficients are tested jointly by a moderated F-test: this results in a single test per protein, hence no contrast has to be defined. The proteins showing a trend are then soft-clustered on the shape of their fitted trajectory.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param time.column String indicating a column of the metadata containing NUMERIC time values (e.g., hours). The values must be coercible to numeric.
#' @param group.column String indicating a column of the metadata defining sample groups (e.g., treatment). When provided, the group x time interaction is tested, meaning that the question becomes 'do the trajectories differ between the groups?'. Default: \code{NULL} (a single trend is tested).
#' @param replicate.column String indicating a column of the metadata containing the replicate IDs. Used only when \code{include.rep.model} is \code{TRUE}. Default: \code{NULL}.
#' @param include.rep.model Logic value indicating whether the intra-replicate correlation should be estimated by \code{limma::duplicateCorrelation} and included in the fit. It requires a \code{replicate.column}. Default: \code{FALSE}.
#' @param time.transform String indicating a transformation to apply to the time values. One among: 'none', 'log2', 'log10', 'log1p', 'log', 'ln', 'sqrt'. The log transformations are computed as log(t + 1), hence time zero is allowed. Strongly recommended ('log2') when the timepoints are log-spaced (e.g., 0, 1, 2, 6, 24 h), otherwise the last timepoint dominates the fit. Default: \code{"none"}.
#' @param spline.df Numeric value indicating the degrees of freedom of the natural-spline basis. When \code{NULL}, it is set to min(3, n.timepoints - 2); with less than 4 timepoints a linear trend (df = 1) is fitted instead. Default: \code{NULL}.
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'imp', 'randomized', 'random'. Default: \code{"imputed"}.
#' @param padj.method String indicating the multiple-testing correction, any method accepted by \code{stats::p.adjust}. Default: \code{"BH"}.
#' @param padj.th Numeric value indicating the adjusted p-value threshold used to call a trend. Default: \code{0.05}.
#' @param log2.amplitude.th Numeric value indicating the minimal amplitude (max - min of the fitted curve, in log2 units) required to call a trend. The proteins passing \code{padj.th} but not this threshold are labelled as 'unresponsive'. Default: \code{0}.
#' @param n.clusters Numeric value indicating the number of clusters. When \code{NULL}, it is estimated by silhouette over \code{k.range}; when \code{0}, the clustering is skipped. Default: \code{NULL}.
#' @param k.range Numeric vector of the candidate numbers of clusters tested when \code{n.clusters} is \code{NULL}. Default: \code{2:8}.
#' @param clustering.method String indicating the clustering algorithm. One among: 'cmeans' (soft, the default), 'pam' (hard, exact) or 'clara' (hard, approximated). CLARA runs PAM on subsamples instead of building the full dissimilarity matrix: it is only useful when the trending proteins are too many for an exact PAM, roughly above 10000, an exact PAM being preferable whenever it runs. Default: \code{"cmeans"}.
#' @param fuzzifier Numeric value indicating the fuzzifier 'm' used by c-means. When \code{NULL}, it is estimated with the \enc{Schwämmle}{Schwammle} & Jensen relation. Default: \code{NULL}.
#' @param clara.samples Numeric value indicating the number of subsamples drawn by CLARA. Used only when \code{clustering.method} is 'clara'. Default: \code{50}.
#' @param clara.sampsize Numeric value indicating the size of each CLARA subsample. When \code{NULL}, it is set to 140 + 2*k, the default of \code{cluster::clara} being too small to be representative. Default: \code{NULL}.
#' @param rank.by String indicating the metric used to rank the proteins within each cluster. One among: 'score' (amplitude * -log10(padj), consistent with the differential score used elsewhere in DEprot), 'membership' (how prototypical a protein is for its cluster, c-means only), 'amplitude', 'padj', or 't.half' (time to half-maximum, which orders the proteins by when they respond). Default: \code{"score"}.
#' @param grid.n Numeric value indicating the number of points at which the fitted curves are evaluated. Default: \code{100}.
#' @param seed Numeric value indicating the seed used for the clustering. Default: \code{1234}.
#' @param verbose Logic value indicating whether the messages should be printed. Default: \code{TRUE}.
#'
#' @details
#' The model spends (\code{spline.df} + 1) coefficients for each level of the \code{group.column}, hence
#' every group must provide at least two samples more than that: a group fitted with as many samples as
#' coefficients is saturated, its residual variance is not estimable, and the moderated F-test would rest
#' entirely on the empirical Bayes prior rather than on the data. When this is not the case the function
#' stops: lowering the \code{spline.df} is usually preferable to dropping timepoints.
#'
#' With less than four timepoints there is no room for any curvature and a linear trend is fitted instead of
#' a spline (a warning is raised). In that case, as well as with a single group, the joint test reduces to a
#' single coefficient: the \code{F.statistic} column then reports the square of the moderated t-statistic,
#' which is the F statistic of that test (1 numerator degree of freedom) and returns exactly the same p-value.
#'
#' @return An object of class \code{DEprot.timecourse}.
#'
#' @references
#' Ritchie M.E., Phipson B., Wu D., \emph{et al.} (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research}, \strong{43}(7), e47. \doi{10.1093/nar/gkv007}
#'
#' \enc{Schwämmle}{Schwammle} V. & Jensen O.N. (2010). A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. \emph{Bioinformatics}, \strong{26}(22), 2841-2848. \doi{10.1093/bioinformatics/btq534}
#'
#' Rousseeuw P.J. (1987). Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. \emph{Journal of Computational and Applied Mathematics}, \strong{20}, 53-65. \doi{10.1016/0377-0427(87)90125-7}
#'
#' Schubert E. & Rousseeuw P.J. (2021). Fast and eager k-medoids clustering: O(k) runtime improvement of the PAM, CLARA, and CLARANS algorithms. \emph{Information Systems}, \strong{101}, 101804. \doi{10.1016/j.is.2021.101804}
#'
#' @seealso \code{\link{rank.timecourse}}, \code{\link{get.timecourse.results}}, \code{\link{plot.timecourse.protein}}, \code{\link{plot.timecourse.profiles}}
#'
#' @importFrom splines ns
#' @importFrom stats model.matrix predict setNames dist
#' @importFrom methods new
#' @importFrom cluster pam clara
#' @importFrom e1071 cmeans
#' @import limma
#'
#' @author Sebastian Gregoricchio
#'
#' @export analyze.timecourse

analyze.timecourse =
  function(DEprot.object,
           time.column,
           group.column = NULL,
           replicate.column = NULL,
           include.rep.model = FALSE,
           time.transform = "none",
           spline.df = NULL,
           which.data = "imputed",
           padj.method = "BH",
           padj.th = 0.05,
           log2.amplitude.th = 0,
           n.clusters = NULL,
           k.range = 2:8,
           clustering.method = "cmeans",
           fuzzifier = NULL,
           clara.samples = 50,
           clara.sampsize = NULL,
           rank.by = "score",
           grid.n = 100,
           seed = 1234,
           verbose = TRUE) {

    ######################################################################################

    ### check object
    if (!(methods::is(DEprot.object, "DEprot.analyses"))) {
      if (!(methods::is(DEprot.object, "DEprot"))) {
        stop("The input must be an object of class 'DEprot' or 'DEprot.analyses'.")
      }
    }

    ### check metadata columns
    meta = DEprot.object@metadata

    if (!(time.column %in% colnames(meta))) {
      stop(paste0("The 'time.column' is not present in the metadata of the object provided.\n",
                  "       Available column IDs: ", paste0(colnames(meta), collapse = ", ")))
    }

    if (!is.null(group.column)) {
      if (!(group.column %in% colnames(meta))) {
        stop(paste0("The 'group.column' is not present in the metadata of the object provided.\n",
                    "       Available column IDs: ", paste0(colnames(meta), collapse = ", ")))
      }
    }


    ### convert the time to numeric: this function does NOT accept categorical timepoints
    time.raw = suppressWarnings(as.numeric(as.character(meta[[time.column]])))

    if (all(is.na(time.raw))) {
      stop(paste0("The column '", time.column, "' cannot be converted to numeric.\n",
                  "       analyze.timecourse() requires numeric timepoints (e.g., 0, 1, 2, 6, 24)."))
    }

    ## the samples without a valid time cannot enter the model and are dropped
    keep = !is.na(time.raw)

    if (any(!keep)) {
      warning(paste0(sum(!keep), " sample(s) have a non-numeric/NA time and are excluded: ",
                     paste(meta$column.id[!keep], collapse = ", ")))
    }

    meta = meta[keep,,drop = FALSE]
    time.raw = time.raw[keep]


    ### extract and subset the counts table
    counts.list = .tc.get.counts(DEprot.object = DEprot.object, which.data = which.data)
    mat = counts.list$mat[,as.character(meta$column.id), drop = FALSE]
    data.used = counts.list$data.used

    timepoints = sort(unique(time.raw))
    n.tp = length(timepoints)

    if (n.tp < 3) {
      stop(paste0("Only ", n.tp, " distinct timepoint(s) found: a time course requires at least 3."))
    }



    ######################################################################################
    ### Build the time basis
    ######################################################################################

    ## the transformation re-spaces the timepoints: with log-spaced designs it prevents
    ## the last timepoint from carrying almost all the leverage of the fit
    tt = .tc.transform.time(time = time.raw, time.transform = time.transform)

    ## a spline needs at least (df + 2) distinct timepoints, otherwise the model becomes
    ## saturated and is strictly equivalent to treating the time as a factor
    if (n.tp < 4) {
      ## below four timepoints there is no room for any curvature, whatever the df asked
      ## for: the fallback to a linear trend is signalled every time the function is the
      ## one taking the decision, that is, unless a linear trend was explicitly required
      if (is.null(spline.df)) {
        warning(paste0("Only ", n.tp, " timepoints available: a linear trend (spline.df = 1) is fitted instead of a spline."))
      } else if (spline.df > 1) {
        warning(paste0("The 'spline.df' provided (", spline.df, ") is too high for ", n.tp, " timepoints: a linear trend (spline.df = 1) is fitted instead of a spline."))
      }

      spline.df = 1

    } else {
      if (is.null(spline.df)) {spline.df = min(3, n.tp - 2)}

      if (spline.df > (n.tp - 2)) {
        warning(paste0("The 'spline.df' provided (", spline.df, ") is too high for ", n.tp, " timepoints: it has been reduced to ", n.tp - 2, ".\n",
                       "       A saturated spline is equivalent to treating the time as a factor."))
        spline.df = n.tp - 2
      }
    }

    ## '.basis.at' re-evaluates the basis on new time values keeping the SAME knots,
    ## which is mandatory to predict the fitted curves on the plotting grid
    if (spline.df == 1) {
      basis = matrix(tt, ncol = 1, dimnames = list(NULL, "1"))
      .basis.at = function(x) {matrix(x, ncol = 1, dimnames = list(NULL, "1"))}
      model.label = "linear"
    } else {
      basis = splines::ns(tt, df = spline.df)
      .basis.at = function(x) {stats::predict(basis, newx = x)}
      model.label = "natural spline"
    }



    ######################################################################################
    ### Design matrix and limma fit
    ######################################################################################

    if (is.null(group.column)) {
      ## single trend: all the spline coefficients are tested together
      design = stats::model.matrix(~ basis)
      test.coefs = 2:ncol(design)
      group.levels = "all"
      group.vector = rep("all", nrow(meta))
    } else {
      ## several groups: only the interaction terms are tested, meaning that the question
      ## becomes whether the SHAPE of the trajectory differs between the groups
      group.vector = factor(as.character(meta[[group.column]]))
      group.levels = levels(group.vector)
      design = stats::model.matrix(~ group.vector * basis)
      test.coefs = grep(":", colnames(design))

      if (length(test.coefs) == 0) {
        stop("No interaction coefficient could be built: please check the 'time.column' and 'group.column' provided.")
      }
    }

    colnames(design) = make.names(colnames(design))
    rownames(design) = as.character(meta$column.id)

    ## The model spends (spline.df + 1) coefficients per group. A group fitted with as many
    ## samples as coefficients is saturated: its residual variance is not estimable and the
    ## moderated F-test would rely on the empirical Bayes prior only, hence at least two
    ## residual degrees of freedom are required within each arm of the design.
    coefs.per.group = spline.df + 1
    n.per.group = table(factor(as.character(group.vector), levels = group.levels))

    if (nrow(design) <= ncol(design) | any(n.per.group < (coefs.per.group + 2))) {
      stop(paste0("Not enough samples (", nrow(design), ") for the ", ncol(design), " coefficients of this model: ",
                  "each group requires at least ", coefs.per.group + 2, " samples (",
                  paste0(names(n.per.group), ": ", as.numeric(n.per.group), collapse = ", "), ").\n",
                  "       Please reduce the 'spline.df', add replicates, or remove the 'group.column'."))
    }

    if (verbose == TRUE) {
      message(paste0("Fitting a ", model.label, " model (df = ", spline.df, ") on ",
                     nrow(mat), " proteins and ", ncol(mat), " samples."))
    }

    ## the replicate can optionally be modelled as a random block (paired designs)
    if (include.rep.model == TRUE & !is.null(replicate.column)) {
      block = as.character(meta[[replicate.column]])
      dupcor = limma::duplicateCorrelation(mat, design = design, block = block)
      fit = limma::lmFit(mat, design = design, block = block, correlation = dupcor$consensus.correlation)
    } else {
      fit = limma::lmFit(mat, design = design)
    }

    fit = limma::eBayes(fit)

    ## a single moderated F-test per protein over all the time coefficients;
    ## 'sort.by = "none"' keeps the row order of the counts table
    top = limma::topTable(fit,
                          coef = test.coefs,
                          number = Inf,
                          sort.by = "none",
                          adjust.method = padj.method)

    ## with a single coefficient (linear trend without groups, or a single interaction term)
    ## topTable returns a moderated t-test and no 'F' column: the F statistic of a test with
    ## one numerator degree of freedom is t^2, and the p-value is identical
    F.statistic = if ("F" %in% colnames(top)) {top$F} else {top$t^2}



    ######################################################################################
    ### Fitted trajectories
    ######################################################################################

    ## dense grid on the transformed scale (where the model lives), then brought back to
    ## the real time scale for the plots and for the kinetic descriptors
    grid.tt = seq(min(tt), max(tt), length.out = grid.n)
    grid.time = .tc.backtransform.time(tt = grid.tt, time.transform = time.transform)

    coefs = fit$coefficients

    ## rebuilds, for a given group, the same columns as the design matrix but evaluated
    ## on the prediction grid instead of on the samples
    .grid.design =
      function(g) {
        b = .basis.at(grid.tt)

        if (is.null(group.column)) {
          X = cbind(1, b)
        } else {
          gm = stats::model.matrix(~ factor(rep(g, grid.n), levels = group.levels))[,-1, drop = FALSE]
          ## the interaction columns must follow the same order as model.matrix(~ g * basis)
          inter = do.call(cbind, lapply(seq_len(ncol(b)), function(j) {gm * b[,j]}))
          X = cbind(1, gm, b, inter)
        }

        colnames(X) = colnames(design)
        return(X)
      }

    fitted.curves =
      stats::setNames(lapply(group.levels,
                             function(g) {
                               m = coefs %*% t(.grid.design(g))
                               rownames(m) = rownames(mat)
                               return(m)
                             }),
                      group.levels)

    ## measured mean at each timepoint, kept for the plots and for the users' own summaries
    observed.means =
      stats::setNames(lapply(group.levels,
                             function(g) {
                               idx.g = if (is.null(group.column)) {rep(TRUE, ncol(mat))} else {as.character(group.vector) == g}

                               m = vapply(timepoints,
                                          function(tp) {
                                            cols = which(idx.g & (time.raw == tp))
                                            if (length(cols) == 0) {return(rep(NA_real_, nrow(mat)))}
                                            return(rowMeans(mat[,cols, drop = FALSE], na.rm = TRUE))
                                          },
                                          FUN.VALUE = numeric(nrow(mat)))

                               colnames(m) = as.character(timepoints)
                               rownames(m) = rownames(mat)
                               return(m)
                             }),
                      group.levels)



    ######################################################################################
    ### Results table
    ######################################################################################

    desc.list =
      stats::setNames(lapply(group.levels,
                             function(g) {.tc.descriptors(curves = fitted.curves[[g]], grid.time = grid.time)}),
                      group.levels)

    results = data.frame(prot.id = rownames(mat),
                         basemean.log2 = rowMeans(mat, na.rm = TRUE),
                         F.statistic = F.statistic,
                         p.value = top$P.Value,
                         padj = top$adj.P.Val,
                         stringsAsFactors = FALSE)

    if (length(group.levels) == 1) {
      results = cbind(results, desc.list[[1]])
    } else {
      ## with several groups the descriptors are group-specific and get a suffix, while the
      ## overall amplitude is defined as the largest one across the groups
      for (g in group.levels) {
        d = desc.list[[g]]
        colnames(d) = paste0(colnames(d), ".", g)
        results = cbind(results, d)
      }

      results$amplitude = apply(do.call(cbind, lapply(desc.list, function(d) {d$amplitude})), 1, max)
    }

    ## 'unresponsive' = significant, but with an amplitude too small to be worth reporting
    results$trend.status = ifelse(results$padj < padj.th & results$amplitude > log2.amplitude.th,
                                  "trending",
                                  ifelse(results$padj < padj.th, "unresponsive", "null"))

    ## same logic as the differential score of DEprot: effect size weighted by significance
    results$score = results$amplitude * -log10(results$padj)



    ######################################################################################
    ### Clustering of the trending proteins
    ######################################################################################

    clusters = NULL
    results$cluster = NA_integer_
    results$membership = NA_real_

    trending = which(results$trend.status == "trending")

    if (!is.null(n.clusters) && n.clusters == 0) {
      if (verbose == TRUE) {message("Clustering skipped ('n.clusters' set to 0).")}

    } else if (length(trending) < 10) {
      warning("Less than 10 trending proteins found: the clustering has been skipped.")

    } else {
      ## the clustering runs on the Z-scored FITTED curves, and not on the raw timepoint
      ## means: this makes it robust to unequal time spacing and to missing values, and it
      ## groups the proteins by the SHAPE of their response rather than by their abundance.
      ## With several groups the trajectories are pasted side by side.
      z = do.call(cbind,
                  lapply(group.levels,
                         function(g) {.tc.zscore(fitted.curves[[g]][trending,, drop = FALSE])}))


      ## the number of clusters, when not provided, is estimated on a subsample (speed)
      if (is.null(n.clusters)) {
        if (verbose == TRUE) {message("Estimating the number of clusters by silhouette...")}

        set.seed(seed)
        sub = z[sample(seq_len(nrow(z)), min(2000, nrow(z))),, drop = FALSE]

        ## the dissimilarity matrix is computed once and reused for every candidate k:
        ## rebuilding it inside the loop is by far the slowest step of the function
        sub.diss = stats::dist(sub)

        sil = vapply(k.range,
                     function(k) {cluster::pam(sub.diss, k = k, diss = TRUE, pamonce = 5)$silinfo$avg.width},
                     FUN.VALUE = numeric(1))

        n.clusters = k.range[which.max(sil)]

        if (verbose == TRUE) {
          message(paste0("   -> k = ", n.clusters, " (average silhouette = ", round(max(sil), 3), ")"))
        }
      }


      set.seed(seed)

      if (tolower(clustering.method) %in% c("cmeans", "c-means", "soft")) {
        ## m = 2 (the classical default) collapses high-dimensional data onto a single
        ## centroid, hence the fuzzifier is estimated from the data when not provided
        if (is.null(fuzzifier)) {fuzzifier = .tc.mestimate(n = nrow(z), d = n.tp)}

        cm = e1071::cmeans(z, centers = n.clusters, m = fuzzifier, iter.max = 500)

        results$cluster[trending] = as.integer(cm$cluster)
        results$membership[trending] = apply(cm$membership, 1, max)

        clusters = list(membership = cm$membership,
                        centroids = cm$centers,
                        method = "cmeans",
                        k = n.clusters,
                        fuzzifier = fuzzifier)

      } else if (tolower(clustering.method) %in% c("pam", "hard")) {
        ## 'pamonce = 5' is the FastPAM reformulation: same medoids, much faster
        pm = cluster::pam(z, k = n.clusters, pamonce = 5)

        results$cluster[trending] = as.integer(pm$clustering)
        ## hard clustering: no membership degree is available
        results$membership[trending] = NA_real_

        clusters = list(membership = NULL,
                        centroids = pm$medoids,
                        method = "pam",
                        k = n.clusters,
                        fuzzifier = NA_real_)

      } else if (tolower(clustering.method) == "clara") {
        ## CLARA is an approximation of PAM: it runs PAM on several subsamples and keeps
        ## the best partition, which avoids the n x n dissimilarity matrix. It is only
        ## worth it on very large sets of trending proteins, an exact PAM being preferable
        ## whenever it runs. Two defaults of the function are deliberately overridden: the
        ## historical 'sampsize' is far too small to be representative, and 'pamLike =
        ## FALSE' makes the swap phase use a criterion different from the one of pam().
        ## 'rngR = TRUE' draws the subsamples from the RNG of R, so that the 'seed' above
        ## makes the result reproducible.
        if (is.null(clara.sampsize)) {clara.sampsize = min(nrow(z), 40 + 2 * n.clusters + 100)}

        cl = cluster::clara(z,
                            k = n.clusters,
                            samples = clara.samples,
                            sampsize = clara.sampsize,
                            pamLike = TRUE,
                            rngR = TRUE)

        results$cluster[trending] = as.integer(cl$clustering)
        ## hard clustering: no membership degree is available
        results$membership[trending] = NA_real_

        clusters = list(membership = NULL,
                        centroids = cl$medoids,
                        method = "clara",
                        k = n.clusters,
                        fuzzifier = NA_real_)

      } else {
        stop("The 'clustering.method' must be one among: 'cmeans', 'pam', 'clara'.")
      }
    }



    ######################################################################################
    ### Build the output object
    ######################################################################################

    tc = methods::new(Class = "DEprot.timecourse",
                      results = results,
                      fitted.curves = fitted.curves,
                      time.grid = grid.time,
                      observed.means = observed.means,
                      timepoints = timepoints,
                      counts.used = mat,
                      tc.metadata = meta,
                      sample.subset = as.character(meta$column.id),
                      data.used = data.used,
                      design = design,
                      clusters = clusters,
                      params = list(time.column = time.column,
                                    group.column = group.column,
                                    replicate.column = replicate.column,
                                    include.rep.model = include.rep.model,
                                    time.transform = time.transform,
                                    spline.df = spline.df,
                                    model = model.label,
                                    padj.method = padj.method,
                                    padj.th = padj.th,
                                    log2.amplitude.th = log2.amplitude.th,
                                    clustering.method = clustering.method,
                                    rank.by = rank.by,
                                    group.levels = group.levels,
                                    grid.n = grid.n,
                                    seed = seed),
                      profile.plot = NULL)

    ## ranking of the proteins, globally and within each cluster
    tc = rank.timecourse(DEprot.timecourse.object = tc, rank.by = rank.by)

    ## the profile plot is stored in the object, as done for the other DEprot classes
    if (!is.null(clusters)) {
      tc@profile.plot = tryCatch(plot.timecourse.profiles(DEprot.timecourse.object = tc),
                                 error = function(e) {return(NULL)})
    }

    return(tc)
  } # END function






# ----------------------------------------------------------------------------------------
###                              RANKING AND GETTERS                                    ###
# ----------------------------------------------------------------------------------------

#' @title rank.timecourse
#'
#' @description Recomputes the ranking of the proteins, globally and within each cluster, without refitting the model (analogous to \code{reapply.thresholds} for the differential analyses).
#'
#' @param DEprot.timecourse.object An object of class \code{DEprot.timecourse}.
#' @param rank.by String indicating the metric used for the ranking. One among: 'score' (amplitude * -log10(padj)), 'membership' (c-means only), 'amplitude', 'padj', or 't.half' (time to half-maximum: the fastest proteins come first, which orders the proteins by WHEN they respond rather than by how much). Default: \code{"score"}.
#'
#' @return The same object, with updated \code{rank.overall} and \code{rank.in.cluster} columns.
#'
#' @importFrom stats na.omit
#'
#' @author Sebastian Gregoricchio
#'
#' @export rank.timecourse

rank.timecourse =
  function(DEprot.timecourse.object,
           rank.by = "score") {

    ### check object
    if (!("DEprot.timecourse" %in% class(DEprot.timecourse.object))) {
      stop("The input must be an object of class 'DEprot.timecourse'.")
    }

    if (!(tolower(rank.by) %in% c("score", "membership", "amplitude", "padj", "t.half"))) {
      stop("The 'rank.by' must be one among: 'score', 'membership', 'amplitude', 'padj', 't.half'.")
    }

    rank.by = tolower(rank.by)
    res = DEprot.timecourse.object@results
    group.levels = DEprot.timecourse.object@params$group.levels

    ## the membership exists only when a soft clustering was performed
    if (rank.by == "membership" & all(is.na(res$membership))) {
      warning("No membership available (hard clustering was used): the ranking is performed by 'score'.")
      rank.by = "score"
    }

    ## for 'score', 'membership' and 'amplitude' the best protein has the HIGHEST value;
    ## for 'padj' and 't.half' it has the LOWEST one, hence the sign is inverted
    metric = switch(rank.by,
                    "score"      = res$score,
                    "membership" = res$membership,
                    "amplitude"  = res$amplitude,
                    "padj"       = -res$padj,
                    "t.half"     = -.tc.descriptor.column(results = res,
                                                          descriptor = "t.half",
                                                          group.levels = group.levels))

    ## a flat curve has no half-maximum: those proteins would end up last whatever happens
    if (rank.by == "t.half" & all(is.na(metric))) {
      warning("No 't.half' could be computed (all the curves are flat): the ranking is performed by 'score'.")
      rank.by = "score"
      metric = res$score
    }

    res$rank.overall = NA_integer_
    res$rank.in.cluster = NA_integer_

    trending = which(res$trend.status == "trending")

    if (length(trending) > 0) {
      res$rank.overall[trending] = rank(-metric[trending], ties.method = "min", na.last = TRUE)

      ## the ranking is then recomputed independently within each cluster
      for (k in unique(stats::na.omit(res$cluster))) {
        idx = which(res$cluster == k)
        res$rank.in.cluster[idx] = rank(-metric[idx], ties.method = "min", na.last = TRUE)
      }
    }

    DEprot.timecourse.object@results = res
    DEprot.timecourse.object@params$rank.by = rank.by

    return(DEprot.timecourse.object)
  } # END function







#' @title get.timecourse.results
#'
#' @description Retrieves the results table of a \code{DEprot.timecourse} object, optionally restricted to specific clusters and/or to the best-ranked proteins.
#'
#' @param DEprot.timecourse.object An object of class \code{DEprot.timecourse}.
#' @param cluster Numeric value (or vector) indicating the cluster(s) to retrieve. Default: \code{NULL} (all the clusters).
#' @param top.n Numeric value indicating how many best-ranked proteins should be kept (within the cluster when a \code{cluster} is provided, globally otherwise). Default: \code{NULL} (all the proteins).
#' @param trending.only Logic value indicating whether only the proteins showing a significant trend should be kept. Default: \code{TRUE}.
#'
#' @return A data.frame sorted by rank.
#'
#' @importFrom utils head
#'
#' @author Sebastian Gregoricchio
#'
#' @export get.timecourse.results

get.timecourse.results =
  function(DEprot.timecourse.object,
           cluster = NULL,
           top.n = NULL,
           trending.only = TRUE) {

    ### check object
    if (!("DEprot.timecourse" %in% class(DEprot.timecourse.object))) {
      stop("The input must be an object of class 'DEprot.timecourse'.")
    }

    res = DEprot.timecourse.object@results

    if (trending.only == TRUE) {res = res[res$trend.status == "trending",, drop = FALSE]}
    if (!is.null(cluster)) {res = res[which(res$cluster %in% cluster),, drop = FALSE]}

    ## when a cluster is requested, the within-cluster ranking is the relevant one
    rank.col = if (!is.null(cluster)) {"rank.in.cluster"} else {"rank.overall"}
    res = res[order(res[[rank.col]]),, drop = FALSE]

    if (!is.null(top.n)) {res = utils::head(res, top.n)}

    rownames(res) = NULL
    return(res)
  } # END function






# ----------------------------------------------------------------------------------------
###                                     PLOTS                                           ###
# ----------------------------------------------------------------------------------------

#' @title plot.timecourse.protein
#'
#' @description Plots the temporal trend of one or more individual proteins: the measured values (one point per sample), the mean +/- SEM at each timepoint, and the trajectory fitted by \code{analyze.timecourse}.
#'
#' @param DEprot.timecourse.object An object of class \code{DEprot.timecourse}.
#' @param protein.id String (or vector of strings) indicating the proteins to plot. The identifiers must correspond to the full row.names of the counts table.
#' @param values String indicating the quantity displayed on the y-axis. One among: 'counts' (the log2 values as they are), 'log2FC' (log2 fold change relative to \code{reference.time}; the counts being already in log2, this is a simple subtraction), 'zscore' (Z-score of the measured values, the fitted curve being centered/scaled with the same mean and SD so that curve and points remain comparable). Aliases are accepted ('fc', 'foldchange', 'log2', 'z', ...). Default: \code{"counts"}.
#' @param reference.time Numeric value indicating the timepoint used as baseline when \code{values} is 'log2FC'. Default: \code{NULL} (the earliest timepoint).
#' @param color String indicating the color of the line and of the points. Used only when no \code{group.column} was defined in the analyses; with several groups use \code{group.colors} instead. Default: \code{"black"}.
#' @param group.colors Named vector of colors to use for the groups. Default: \code{NULL} (automatic).
#' @param show.points Logic value indicating whether the individual samples should be displayed. Default: \code{TRUE}.
#' @param show.fit Logic value indicating whether the fitted trajectory should be displayed. Default: \code{TRUE}.
#' @param shape.column String indicating a column from the metadata table to use as factor for the shape of the points (e.g., 'replicate'). The column name is used as title of the shape legend. Default: \code{NULL}.
#' @param log.x Logic value indicating whether the x-axis should use a log10(time + 1) scale, useful with log-spaced timepoints. Default: \code{FALSE}.
#' @param show.stats Logic value indicating whether padj, amplitude, shape, cluster and rank should be added as subtitle (single protein only). Default: \code{TRUE}.
#' @param panel.border Logic value indicating whether a border should be drawn around each panel. Default: \code{TRUE}.
#' @param ncol Numeric value indicating the number of columns of the facet grid when several proteins are plotted. Default: \code{NULL} (automatic).
#' @param line.size Numeric value indicating the thickness of the fitted line. Default: \code{0.8}.
#' @param point.size Numeric value indicating the size of the mean points. Default: \code{2}.
#' @param scale.expression Deprecated, kept for backward compatibility: \code{TRUE} is equivalent to \code{values = "zscore"}, \code{FALSE} to \code{values = "counts"}. Default: \code{NULL}.
#'
#' @return A plot of class ggplot2, facetted by protein when several proteins are provided.
#'
#' @seealso \code{\link{analyze.timecourse}}, \code{\link{plot.timecourse.profiles}}
#'
#' @import ggplot2
#' @import ggtext
#' @importFrom ggpubr theme_pubr
#' @importFrom stats aggregate sd setNames
#' @importFrom scales trans_new
#'
#' @author Sebastian Gregoricchio
#'
#' @export plot.timecourse.protein

plot.timecourse.protein =
  function(DEprot.timecourse.object,
           protein.id,
           values = "counts",
           reference.time = NULL,
           color = "black",
           group.colors = NULL,
           show.points = TRUE,
           show.fit = TRUE,
           shape.column = NULL,
           log.x = FALSE,
           show.stats = TRUE,
           panel.border = TRUE,
           ncol = NULL,
           line.size = 0.8,
           point.size = 2,
           scale.expression = NULL) {

    ######################################################################################

    ### check object
    if (!("DEprot.timecourse" %in% class(DEprot.timecourse.object))) {
      stop("The input must be an object of class 'DEprot.timecourse'.")
    }

    tc = DEprot.timecourse.object
    p = tc@params
    mat = tc@counts.used
    meta = tc@tc.metadata


    ### deprecated argument
    if (!is.null(scale.expression)) {
      warning("The option 'scale.expression' is deprecated: please use values = 'zscore' / 'counts' instead.")
      values = ifelse(scale.expression == TRUE, "zscore", "counts")
    }

    values = .tc.check.values(values)


    ### check proteins
    missing.prot = setdiff(protein.id, rownames(mat))

    if (length(missing.prot) > 0) {
      stop(paste0("The following protein(s) are not present in the dataset: ",
                  paste(missing.prot, collapse = ", ")))
    }


    time.raw = suppressWarnings(as.numeric(as.character(meta[[p$time.column]])))
    group.vec = if (is.null(p$group.column)) {rep("all", nrow(meta))} else {as.character(meta[[p$group.column]])}

    ref = .tc.ref.index(time.grid = tc@time.grid, timepoints = tc@timepoints, reference.time = reference.time)



    ######################################################################################
    ### Measured values
    ######################################################################################

    points.tb =
      do.call(rbind,
              lapply(protein.id,
                     function(pr) {
                       v = as.numeric(mat[pr,])

                       if (values == "zscore") {
                         s = stats::sd(v, na.rm = TRUE)
                         v = if (is.na(s) | s == 0) {v - mean(v, na.rm = TRUE)} else {(v - mean(v, na.rm = TRUE)) / s}

                       } else if (values == "log2FC") {
                         ## the baseline is the MEAN measured value at the reference time,
                         ## computed separately within each group
                         for (g in unique(group.vec)) {
                           idx.g = which(group.vec == g)
                           idx.ref = which(group.vec == g & time.raw == ref$time)

                           if (length(idx.ref) == 0) {
                             warning(paste0("No sample available at time ", ref$time, " for the group '", g,
                                            "': the log2FC could not be computed for this group."))
                             next
                           }

                           v[idx.g] = v[idx.g] - mean(as.numeric(mat[pr,idx.ref]), na.rm = TRUE)
                         }
                       }

                       data.frame(prot.id = pr,
                                  column.id = colnames(mat),
                                  time = time.raw,
                                  group = group.vec,
                                  value = v,
                                  shape = if (is.null(shape.column)) {NA_character_} else {as.character(meta[[shape.column]])},
                                  stringsAsFactors = FALSE)
                     }))


    ## mean and standard error of the mean at each timepoint
    mean.tb = stats::aggregate(value ~ prot.id + time + group,
                               data = points.tb,
                               FUN = function(x) {c(avg = mean(x, na.rm = TRUE),
                                                    sem = stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))})
    mean.tb = do.call(data.frame, mean.tb)
    colnames(mean.tb)[(ncol(mean.tb)-1):ncol(mean.tb)] = c("avg", "sem")



    ######################################################################################
    ### Fitted trajectories
    ######################################################################################

    fit.tb =
      do.call(rbind,
              lapply(protein.id,
                     function(pr) {
                       do.call(rbind,
                               lapply(names(tc@fitted.curves),
                                      function(g) {
                                        v = as.numeric(tc@fitted.curves[[g]][pr,])

                                        if (values == "zscore") {
                                          ## scaled with the mean/SD of the MEASURED values, so that
                                          ## the curve stays on the same scale as the points
                                          refv = as.numeric(mat[pr,])
                                          s = stats::sd(refv, na.rm = TRUE)
                                          v = if (is.na(s) | s == 0) {v - mean(refv, na.rm = TRUE)} else {(v - mean(refv, na.rm = TRUE)) / s}

                                        } else if (values == "log2FC") {
                                          ## the curve is shifted by its own value at the reference
                                          ## time, hence it is not forced to pass through the data
                                          v = v - v[ref$index]
                                        }

                                        data.frame(prot.id = pr,
                                                   group = g,
                                                   time = tc@time.grid,
                                                   value = v,
                                                   stringsAsFactors = FALSE)
                                      }))
                     }))



    ######################################################################################
    ### Plot
    ######################################################################################

    plot.tc = ggplot2::ggplot()

    ## a fold change is read against 0
    if (values == "log2FC") {
      plot.tc = plot.tc + ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "gray40", linewidth = 0.3)
    }

    if (show.fit == TRUE) {
      plot.tc =
        plot.tc +
        ggplot2::geom_line(data = fit.tb,
                           ggplot2::aes(x = time, y = value, color = group),
                           linewidth = line.size)
    }

    plot.tc =
      plot.tc +
      ggplot2::geom_errorbar(data = mean.tb,
                             ggplot2::aes(x = time, ymin = avg - sem, ymax = avg + sem, color = group),
                             width = 0,
                             linewidth = 0.4) +
      ggplot2::geom_point(data = mean.tb,
                          ggplot2::aes(x = time, y = avg, color = group),
                          size = point.size)

    if (show.points == TRUE) {
      plot.tc =
        plot.tc +
        ggplot2::geom_point(data = points.tb,
                            ggplot2::aes(x = time, y = value, color = group, shape = shape),
                            size = point.size * 0.6,
                            alpha = 0.45)
    }


    ### colors: 'group.colors' when the groups are defined, a single 'color' otherwise
    if (!is.null(group.colors)) {
      plot.tc = plot.tc + ggplot2::scale_color_manual(values = group.colors)
    } else if (is.null(p$group.column)) {
      ## all the samples share the dummy group 'all', hence a single value is enough
      plot.tc = plot.tc + ggplot2::scale_color_manual(values = stats::setNames(color, "all"))
    }


    ### legends
    if (is.null(p$group.column)) {plot.tc = plot.tc + ggplot2::guides(color = "none")}

    if (is.null(shape.column)) {
      plot.tc = plot.tc + ggplot2::guides(shape = "none")
    } else {
      ## the keys are enlarged, the plotted points being deliberately small
      plot.tc = plot.tc + ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
    }


    ## with log-spaced designs a linear axis crushes all the early timepoints together
    if (log.x == TRUE) {
      plot.tc =
        plot.tc +
        ggplot2::scale_x_continuous(trans = scales::trans_new(name = "log1p",
                                                              transform = function(x) {log10(x + 1)},
                                                              inverse = function(x) {10^x - 1}),
                                    breaks = tc@timepoints)
    } else {
      plot.tc = plot.tc + ggplot2::scale_x_continuous(breaks = tc@timepoints)
    }


    ## statistics displayed as subtitle (only when a single protein is plotted)
    sub.title = NULL

    if (show.stats == TRUE & length(protein.id) == 1) {
      r = tc@results[tc@results$prot.id == protein.id,]

      sub.title = paste0("padj = ", format(r$padj, digits = 3, scientific = TRUE),
                         " | amplitude = ", round(r$amplitude, 2), " log~2~",
                         " | shape: ", r$trend.shape,
                         ifelse(is.na(r$cluster),
                                "",
                                paste0("<br>cluster ", r$cluster, " (rank ", r$rank.in.cluster,
                                       "/", sum(tc@results$cluster == r$cluster, na.rm = TRUE), ")")))
    }


    plot.tc =
      plot.tc +
      ggplot2::labs(x = paste0("time (", p$time.column, ")"),
                    y = .tc.ylab(values = values, ref.time = ref$time, time.column = p$time.column),
                    ## the metadata column name is used as title of the shape legend
                    shape = shape.column,
                    title = if (length(protein.id) == 1) {paste0("**", protein.id, "**")} else {NULL},
                    subtitle = sub.title) +
      ggpubr::theme_pubr(legend = "right") +
      ggplot2::theme(aspect.ratio = 0.8,
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.title.y = ggtext::element_markdown(),
                     plot.title = ggtext::element_markdown(hjust = 0.5),
                     plot.subtitle = ggtext::element_markdown(hjust = 0.5, size = 8),
                     ## naked facet labels, the protein name in bold
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face = "bold"))

    ## the border replaces the axis lines drawn by theme_pubr, which would otherwise
    ## be doubled on the left and bottom sides of each panel
    if (panel.border == TRUE) {
      plot.tc =
        plot.tc +
        ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                       axis.line = ggplot2::element_blank())
    }

    if (length(protein.id) > 1) {
      plot.tc = plot.tc + ggplot2::facet_wrap(~ prot.id, scales = "free_y", ncol = ncol)
    }

    return(plot.tc)
  } # END function






#' @title plot.timecourse.profiles
#'
#' @description Plots the fitted trajectories of the trending proteins, one panel per cluster. The lines are colored by cluster membership (c-means) and the cluster centroid (trend line) is overlaid. With \code{values = "zscore"} this reproduces the classical Mfuzz-like representation.
#'
#' @param DEprot.timecourse.object An object of class \code{DEprot.timecourse}.
#' @param clusters Numeric value (or vector) indicating the cluster(s) to plot. Default: \code{NULL} (all the clusters).
#' @param top.n Numeric value indicating how many best-ranked proteins of each cluster should be displayed. Default: \code{NULL} (all the proteins).
#' @param values String indicating the quantity displayed on the y-axis. One among: 'zscore' (Z-score of each fitted trajectory: shape only, and the scale on which the clustering was actually performed), 'log2FC' (log2 fold change of the fitted curve relative to \code{reference.time}, which preserves the real amplitude of the response), 'counts' (the fitted log2 values as they are: absolute abundances, hence the clusters will overlap vertically). Aliases are accepted ('fc', 'foldchange', 'log2', 'z', ...). Default: \code{"zscore"}.
#' @param reference.time Numeric value indicating the timepoint used as baseline when \code{values} is 'log2FC'. Default: \code{NULL} (the earliest timepoint).
#' @param free.y Logic value indicating whether each facet should have its own y-axis. Default: \code{NULL}, meaning \code{TRUE} when \code{values} is 'counts' and \code{FALSE} otherwise.
#' @param line.color String indicating a single color for the individual protein lines. When provided, it replaces the membership gradient. Default: \code{NULL} (the lines are colored by membership).
#' @param membership.palette Vector of colors used for the membership gradient, ignored when \code{line.color} is provided. Default: \code{c("#3B4CC0", "#7B9FF9", "#EDD1C2", "#E36C55", "#B40426")}.
#' @param centroid.color String indicating the color of the trend line, that is, of the cluster centroid. Default: \code{"black"}.
#' @param centroid.size Numeric value indicating the thickness of the trend line. Default: \code{1}.
#' @param line.alpha Numeric value indicating the transparency of the individual protein lines. Default: \code{0.4}.
#' @param log.x Logic value indicating whether the x-axis should use a log10(time + 1) scale. Default: \code{FALSE}.
#' @param panel.border Logic value indicating whether a border should be drawn around each panel. Default: \code{TRUE}.
#' @param ncol Numeric value indicating the number of columns of the facet grid. Default: \code{NULL} (automatic).
#'
#' @return A plot of class ggplot2.
#'
#' @seealso \code{\link{analyze.timecourse}}, \code{\link{plot.timecourse.protein}}
#'
#' @import ggplot2
#' @import ggtext
#' @importFrom ggpubr theme_pubr
#' @importFrom stats aggregate
#' @importFrom utils head
#' @importFrom scales trans_new
#'
#' @author Sebastian Gregoricchio
#'
#' @export plot.timecourse.profiles

plot.timecourse.profiles =
  function(DEprot.timecourse.object,
           clusters = NULL,
           top.n = NULL,
           values = "zscore",
           reference.time = NULL,
           free.y = NULL,
           line.color = NULL,
           membership.palette = c("#3B4CC0", "#7B9FF9", "#EDD1C2", "#E36C55", "#B40426"),
           centroid.color = "black",
           centroid.size = 1,
           line.alpha = 0.4,
           log.x = FALSE,
           panel.border = TRUE,
           ncol = NULL) {

    ######################################################################################

    ### check object
    if (!("DEprot.timecourse" %in% class(DEprot.timecourse.object))) {
      stop("The input must be an object of class 'DEprot.timecourse'.")
    }

    tc = DEprot.timecourse.object

    if (is.null(tc@clusters)) {
      stop("No clustering is available in this object: nothing to plot.")
    }

    values = .tc.check.values(values)

    ## the absolute abundances span the whole dynamic range: a shared y-axis would flatten
    ## every cluster, hence the free axis is used by default in that case only
    if (is.null(free.y)) {free.y = (values == "counts")}

    ref = .tc.ref.index(time.grid = tc@time.grid, timepoints = tc@timepoints, reference.time = reference.time)



    ######################################################################################
    ### Select the proteins to display
    ######################################################################################

    res = tc@results[tc@results$trend.status == "trending",, drop = FALSE]

    if (!is.null(clusters)) {res = res[which(res$cluster %in% clusters),, drop = FALSE]}

    ## only the best-ranked proteins of each cluster are kept (useful with large datasets)
    if (!is.null(top.n)) {
      res = do.call(rbind,
                    lapply(split(res, res$cluster),
                           function(x) {utils::head(x[order(x$rank.in.cluster),], top.n)}))
    }

    if (nrow(res) == 0) {
      stop("No protein left to plot with the 'clusters'/'top.n' selection provided.")
    }

    ## a hard clustering provides no membership degree: without a fixed color the lines
    ## would all be mapped to NA and drawn in grey by default
    if (is.null(line.color) & all(is.na(res$membership))) {line.color = "grey50"}



    ######################################################################################
    ### Reshape the fitted curves on the requested scale
    ######################################################################################

    group.levels = names(tc@fitted.curves)

    profiles.tb =
      do.call(rbind,
              lapply(group.levels,
                     function(g) {
                       curves = tc@fitted.curves[[g]][res$prot.id,, drop = FALSE]

                       ## the counts being already in log2, the fold change is a subtraction
                       m = switch(values,
                                  "zscore" = .tc.zscore(curves),
                                  "log2FC" = curves - curves[,ref$index],
                                  "counts" = curves)

                       data.frame(prot.id = rep(res$prot.id, times = ncol(m)),
                                  cluster = rep(res$cluster, times = ncol(m)),
                                  membership = rep(res$membership, times = ncol(m)),
                                  group = g,
                                  time = rep(tc@time.grid, each = nrow(m)),
                                  value = as.numeric(m),
                                  stringsAsFactors = FALSE)
                     }))

    ## the centroid is recomputed on the scale that is actually displayed
    centroid.tb = stats::aggregate(value ~ cluster + group + time, data = profiles.tb, FUN = mean)

    ## facet labels carrying the size of each cluster
    n.per.cluster = table(res$cluster)
    profiles.tb$facet = paste0("cluster ", profiles.tb$cluster, " (n = ", n.per.cluster[as.character(profiles.tb$cluster)], ")")
    centroid.tb$facet = paste0("cluster ", centroid.tb$cluster, " (n = ", n.per.cluster[as.character(centroid.tb$cluster)], ")")



    ######################################################################################
    ### Plot
    ######################################################################################

    plot.profiles = ggplot2::ggplot()

    if (values == "log2FC") {
      plot.profiles = plot.profiles + ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "gray40", linewidth = 0.3)
    }

    ## the individual trajectories: either colored by membership, or in a single color
    if (is.null(line.color)) {
      plot.profiles =
        plot.profiles +
        ggplot2::geom_line(data = profiles.tb,
                           ggplot2::aes(x = time, y = value, group = interaction(prot.id, group), color = membership),
                           alpha = line.alpha,
                           linewidth = 0.25) +
        ggplot2::scale_color_gradientn(colors = membership.palette, name = "membership")
    } else {
      plot.profiles =
        plot.profiles +
        ggplot2::geom_line(data = profiles.tb,
                           ggplot2::aes(x = time, y = value, group = interaction(prot.id, group)),
                           color = line.color,
                           alpha = line.alpha,
                           linewidth = 0.25)
    }

    plot.profiles =
      plot.profiles +
      ## trend line = centroid of the cluster, on the scale actually displayed
      ggplot2::geom_line(data = centroid.tb,
                         ggplot2::aes(x = time, y = value, group = group),
                         color = centroid.color,
                         linewidth = centroid.size) +
      ggplot2::facet_wrap(~ facet,
                          ncol = ncol,
                          scales = ifelse(free.y == TRUE, "free_y", "fixed")) +
      ggplot2::labs(x = paste0("time (", tc@params$time.column, ")"),
                    y = .tc.ylab(values = values, ref.time = ref$time, time.column = tc@params$time.column, fitted = TRUE)) +
      ggpubr::theme_pubr(legend = "right") +
      ggplot2::theme(aspect.ratio = 0.8,
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.title.y = ggtext::element_markdown(),
                     ## naked facet labels, the cluster name in bold
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face = "bold"))

    ## the border replaces the axis lines drawn by theme_pubr
    if (panel.border == TRUE) {
      plot.profiles =
        plot.profiles +
        ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                       axis.line = ggplot2::element_blank())
    }

    ## with several groups the trajectories are distinguished by line type
    if (length(group.levels) > 1) {
      plot.profiles = plot.profiles + ggplot2::aes(linetype = group)
    }

    if (log.x == TRUE) {
      plot.profiles =
        plot.profiles +
        ggplot2::scale_x_continuous(trans = scales::trans_new(name = "log1p",
                                                              transform = function(x) {log10(x + 1)},
                                                              inverse = function(x) {10^x - 1}),
                                    breaks = tc@timepoints)
    } else {
      plot.profiles = plot.profiles + ggplot2::scale_x_continuous(breaks = tc@timepoints)
    }

    return(plot.profiles)
  } # END function






#' @title heatmap.timecourse
#'
#' @description Plots a heatmap of the proteins showing a temporal trend, with the rows split by cluster. The values displayed can be the measured mean at each timepoint or the smooth trajectory fitted by \code{analyze.timecourse}.
#'
#' @param DEprot.timecourse.object An object of class \code{DEprot.timecourse}.
#' @param clusters Numeric value (or vector) indicating the cluster(s) to display. Default: \code{NULL} (all the clusters).
#' @param top.n Numeric value indicating how many best-ranked proteins of each cluster should be displayed. Default: \code{NULL} (all the proteins).
#' @param values String indicating the quantity displayed. One among: 'zscore' (row-wise Z-score, the scale on which the clustering was performed), 'log2FC' (log2 fold change relative to \code{reference.time}; the counts being already in log2, this is a simple subtraction), 'counts' (the log2 values as they are). Aliases are accepted ('fc', 'foldchange', 'log2', 'z', ...). Default: \code{"zscore"}.
#' @param reference.time Numeric value indicating the timepoint used as baseline when \code{values} is 'log2FC'. Default: \code{NULL} (the earliest timepoint).
#' @param use.fitted Logic value indicating whether the fitted trajectories should be displayed (smooth heatmap, one column per grid point) instead of the mean measured value at each timepoint (one column per timepoint). Default: \code{FALSE}.
#' @param order.by String indicating how the proteins should be sorted within each cluster. One among: 'rank' (the ranking stored in the object), 'membership' (most prototypical first), 'peak.time' (sorts the proteins by the moment of their maximal change, generating a wave-like pattern), 'amplitude' or 'hclust' (hierarchical clustering of the profiles). Default: \code{"rank"}.
#' @param group.subset String (or vector) indicating which group levels should be displayed when a \code{group.column} was used. Default: \code{NULL} (all the groups).
#' @param palette Vector of colors used when \code{values} is 'counts'. Default: \code{NULL}, meaning that the 'mako' viridis palette is used.
#' @param low.color String indicating the color of the lowest values (only for 'zscore' and 'log2FC'). Default: \code{"#2166AC"}.
#' @param mid.color String indicating the color of the values at zero (only for 'zscore' and 'log2FC'). Default: \code{"white"}.
#' @param high.color String indicating the color of the highest values (only for 'zscore' and 'log2FC'). Default: \code{"firebrick"}.
#' @param na.color String indicating the color used for the missing values. Default: \code{"gray"}.
#' @param color.limits Numeric vector of two elements indicating the limits of the color scale. Values outside the range are squished to the extremes. Default: \code{c(NA, NA)} (automatic and symmetric around zero).
#' @param cell.border.color String indicating the color of the cell borders (ignored when \code{use.fitted} is \code{TRUE}). Default: \code{NA}.
#' @param cell.border.width Numeric value indicating the width of the cell borders. Default: \code{0.5}.
#' @param panel.border Logic value indicating whether a border should be drawn around each facet, meaning around each cluster block. Default: \code{TRUE}.
#' @param panel.border.color String indicating the color of the panel border. Default: \code{"black"}.
#' @param show.protein.names Logic value indicating whether the protein names should be displayed on the y-axis. Default: \code{FALSE}.
#' @param protein.names.pattern String indicating a pattern to be passed to gsub in order to shorten the protein names displayed. Default: \code{NULL} (no changes in the IDs).
#' @param show.cluster.size Logic value indicating whether the number of proteins should be added to the cluster labels. Default: \code{TRUE}.
#' @param title String indicating the title of the plot. Default: \code{NULL} (automatic).
#'
#' @return A plot of class ggplot2.
#'
#' @seealso \code{\link{analyze.timecourse}}, \code{\link{plot.timecourse.profiles}}, \code{\link{timecourse.enrichment}}
#'
#' @import ggplot2
#' @import ggtext
#' @import viridis
#' @importFrom ggpubr theme_pubr
#' @importFrom stats dist hclust
#' @importFrom scales squish
#' @importFrom utils head
#'
#' @author Sebastian Gregoricchio
#'
#' @export heatmap.timecourse

heatmap.timecourse =
  function(DEprot.timecourse.object,
           clusters = NULL,
           top.n = NULL,
           values = "zscore",
           reference.time = NULL,
           use.fitted = FALSE,
           order.by = "rank",
           group.subset = NULL,
           palette = NULL,
           low.color = "#2166AC",
           mid.color = "white",
           high.color = "firebrick",
           na.color = "gray",
           color.limits = c(NA, NA),
           cell.border.color = NA,
           cell.border.width = 0.5,
           panel.border = TRUE,
           panel.border.color = "black",
           show.protein.names = FALSE,
           protein.names.pattern = NULL,
           show.cluster.size = TRUE,
           title = NULL) {

    ######################################################################################

    ### check object
    if (!("DEprot.timecourse" %in% class(DEprot.timecourse.object))) {
      stop("The input must be an object of class 'DEprot.timecourse'.")
    }

    tc = DEprot.timecourse.object

    if (is.null(tc@clusters)) {
      stop("No clustering is available in this object: nothing to plot.")
    }

    values = .tc.check.values(values)

    if (!(tolower(order.by) %in% c("rank", "membership", "peak.time", "amplitude", "hclust"))) {
      stop("The 'order.by' must be one among: 'rank', 'membership', 'peak.time', 'amplitude', 'hclust'.")
    }

    order.by = tolower(order.by)

    ref = .tc.ref.index(time.grid = tc@time.grid, timepoints = tc@timepoints, reference.time = reference.time)



    ######################################################################################
    ### Select the proteins to display
    ######################################################################################

    res = tc@results[tc@results$trend.status == "trending",, drop = FALSE]

    if (!is.null(clusters)) {res = res[which(res$cluster %in% clusters),, drop = FALSE]}

    ## only the best-ranked proteins of each cluster are kept (useful with large datasets)
    if (!is.null(top.n)) {
      res = do.call(rbind,
                    lapply(split(res, res$cluster),
                           function(x) {utils::head(x[order(x$rank.in.cluster),], top.n)}))
    }

    if (nrow(res) == 0) {
      stop("No protein left to plot with the 'clusters'/'top.n' selection provided.")
    }


    ### select the groups to display
    group.levels = names(tc@fitted.curves)

    if (!is.null(group.subset)) {
      if (!all(group.subset %in% group.levels)) {
        stop(paste0("The 'group.subset' provided is not available in this object.\n",
                    "       Available groups: ", paste(group.levels, collapse = ", ")))
      }
      group.levels = group.subset
    }



    ######################################################################################
    ### Build the matrices to display
    ######################################################################################

    ## either the smooth fitted trajectories (many columns) or the mean measured value at
    ## each timepoint (one column per timepoint, the classical representation)
    mat.list =
      lapply(group.levels,
             function(g) {
               if (use.fitted == TRUE) {
                 m = tc@fitted.curves[[g]][res$prot.id,, drop = FALSE]
                 colnames(m) = as.character(tc@time.grid)
                 ref.col = ref$index
               } else {
                 m = tc@observed.means[[g]][res$prot.id,, drop = FALSE]
                 ref.col = which.min(abs(as.numeric(colnames(m)) - ref$time))
               }

               ## the counts being already in log2, the fold change is a subtraction
               m = switch(values,
                          "zscore" = .tc.zscore(m),
                          "log2FC" = m - m[,ref.col],
                          "counts" = m)

               return(m)
             })

    names(mat.list) = group.levels



    ######################################################################################
    ### Define the row order
    ######################################################################################

    ## the proteins are always grouped by cluster; the sorting acts WITHIN each cluster.
    ## 'peak.time' is the most informative one for a time course: it generates the
    ## characteristic wave, since the proteins get sorted by when they respond.
    protein.order =
      unlist(lapply(sort(unique(res$cluster)),
                    function(k) {
                      sub = res[which(res$cluster == k),, drop = FALSE]

                      if (order.by == "hclust") {
                        ## hierarchical clustering of the profiles of the first group only
                        m.k = mat.list[[1]][sub$prot.id,, drop = FALSE]
                        if (nrow(m.k) < 3) {return(sub$prot.id)}
                        hc = stats::hclust(stats::dist(m.k), method = "complete")
                        return(sub$prot.id[hc$order])
                      }

                      ## for all the other metrics a simple sorting is enough
                      metric = switch(order.by,
                                      "rank"       = sub$rank.in.cluster,
                                      "membership" = -sub$membership,
                                      "peak.time"  = sub$peak.time,
                                      "amplitude"  = -sub$amplitude)

                      return(sub$prot.id[order(metric)])
                    }))

    ## the y-axis is reversed so that the first protein appears on top of the heatmap
    protein.order = rev(protein.order)



    ######################################################################################
    ### Reshape into a long table
    ######################################################################################

    heatmap.tb =
      do.call(rbind,
              lapply(group.levels,
                     function(g) {
                       m = mat.list[[g]]

                       data.frame(prot.id = rep(rownames(m), times = ncol(m)),
                                  cluster = rep(res$cluster, times = ncol(m)),
                                  group = g,
                                  time = rep(as.numeric(colnames(m)), each = nrow(m)),
                                  value = as.numeric(m),
                                  stringsAsFactors = FALSE)
                     }))

    heatmap.tb$prot.id = factor(heatmap.tb$prot.id, levels = protein.order)

    ## facet labels, optionally carrying the size of each cluster
    n.per.cluster = table(res$cluster)

    if (show.cluster.size == TRUE) {
      heatmap.tb$facet.row = paste0("cluster ", heatmap.tb$cluster,
                                    "\n(n = ", n.per.cluster[as.character(heatmap.tb$cluster)], ")")
    } else {
      heatmap.tb$facet.row = paste0("cluster ", heatmap.tb$cluster)
    }

    heatmap.tb$facet.row = factor(heatmap.tb$facet.row,
                                  levels = unique(heatmap.tb$facet.row[order(heatmap.tb$cluster)]))

    ## when the timepoints are used as columns they are kept as a discrete axis, so that
    ## the cells have all the same width whatever the real spacing of the design
    if (use.fitted == FALSE) {
      heatmap.tb$time = factor(heatmap.tb$time, levels = as.character(tc@timepoints))
    }



    ######################################################################################
    ### Color scale
    ######################################################################################

    ## for Z-scores and fold changes the scale must be diverging and centered on zero
    if (values %in% c("zscore", "log2FC")) {
      if (all(is.na(color.limits))) {
        max.abs = max(abs(heatmap.tb$value), na.rm = TRUE)
        color.limits = c(-max.abs, max.abs)
      }

      fill.scale = ggplot2::scale_fill_gradient2(low = low.color,
                                                 mid = mid.color,
                                                 high = high.color,
                                                 midpoint = 0,
                                                 limits = color.limits,
                                                 oob = scales::squish,
                                                 na.value = na.color,
                                                 name = .tc.ylab(values = values,
                                                                 ref.time = ref$time,
                                                                 time.column = tc@params$time.column))
    } else {
      if (is.null(palette)) {palette = viridis::viridis(n = 100, option = "mako", direction = -1)}

      fill.scale = ggplot2::scale_fill_gradientn(colors = palette,
                                                 limits = color.limits,
                                                 oob = scales::squish,
                                                 na.value = na.color,
                                                 name = .tc.ylab(values = values,
                                                                 ref.time = ref$time,
                                                                 time.column = tc@params$time.column))
    }



    ######################################################################################
    ### Plot
    ######################################################################################

    plot.heatmap = ggplot2::ggplot(data = heatmap.tb, ggplot2::aes(x = time, y = prot.id, fill = value))

    ## geom_raster is much lighter when the fitted grid is displayed (100 columns)
    if (use.fitted == TRUE) {
      plot.heatmap = plot.heatmap + ggplot2::geom_raster()
    } else {
      plot.heatmap = plot.heatmap + ggplot2::geom_tile(color = cell.border.color, linewidth = cell.border.width)
    }

    ## 'space = "free_y"' makes the height of each cluster block proportional to its size
    if (length(group.levels) > 1) {
      plot.heatmap = plot.heatmap + ggplot2::facet_grid(rows = ggplot2::vars(facet.row),
                                                        cols = ggplot2::vars(group),
                                                        scales = "free_y",
                                                        space = "free_y",
                                                        switch = "y")
    } else {
      plot.heatmap = plot.heatmap + ggplot2::facet_grid(rows = ggplot2::vars(facet.row),
                                                        scales = "free_y",
                                                        space = "free_y",
                                                        switch = "y")
    }

    if (is.null(title)) {
      title = paste0("**Time course**: ", nrow(res), " proteins in ",
                     length(unique(res$cluster)), " clusters")
    }

    plot.heatmap =
      plot.heatmap +
      fill.scale +
      ggplot2::labs(x = paste0("time (", tc@params$time.column, ")"),
                    y = NULL,
                    title = title) +
      ggpubr::theme_pubr(legend = "right") +
      ggplot2::theme(plot.title = ggtext::element_markdown(hjust = 0.5),
                     legend.title = ggtext::element_markdown(),
                     axis.line = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
                     strip.placement = "outside",
                     strip.text.y.left = ggplot2::element_text(angle = 0, size = 8),
                     panel.spacing.y = ggplot2::unit(3, "pt"))

    ## each cluster block gets its own frame
    if (panel.border == TRUE) {
      plot.heatmap = plot.heatmap + ggplot2::theme(panel.border = ggplot2::element_rect(color = panel.border.color,
                                                                                        fill = NA,
                                                                                        linewidth = 0.5))
    }



    ######################################################################################
    ### Axes
    ######################################################################################

    ## a single scale_y_discrete, a second call would override the first one. The tiles
    ## being one unit wide and centered on the ticks, an ADDITIVE expansion of 0.5 stops
    ## the panel exactly on the edge of the first and last cells: 0 would cut them in half
    y.labels = if (!is.null(protein.names.pattern)) {
      function(x) {gsub(protein.names.pattern, "", x)}
    } else {
      ggplot2::waiver()
    }

    plot.heatmap = plot.heatmap + ggplot2::scale_y_discrete(labels = y.labels,
                                                            expand = ggplot2::expansion(add = 0.5))

    ## with thousands of proteins the names are unreadable and are hidden by default
    if (show.protein.names == TRUE) {
      plot.heatmap = plot.heatmap + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))
    } else {
      plot.heatmap = plot.heatmap + ggplot2::theme(axis.text.y = ggplot2::element_blank())
    }

    ## the time axis is moved on top of the heatmap; with the discrete axis the ticks are
    ## meaningless, the labels being already centered on their own column
    if (use.fitted == TRUE) {
      plot.heatmap = plot.heatmap + ggplot2::scale_x_continuous(breaks = tc@timepoints,
                                                                expand = c(0,0),
                                                                position = "top")
    } else {
      plot.heatmap =
        plot.heatmap +
        ggplot2::scale_x_discrete(position = "top",
                                  expand = ggplot2::expansion(add = 0.5)) +
        ggplot2::theme(axis.ticks.x = ggplot2::element_blank())
    }

    return(plot.heatmap)
  } # END function






# ----------------------------------------------------------------------------------------
###                                    ENRICHMENTS                                     ###
# ----------------------------------------------------------------------------------------

#' @title timecourse.enrichment
#'
#' @description Performs an OverRepresentation Analysis (ORA) independently for each cluster of a time-course analysis, using as background (universe) all the proteins that were tested. It returns the combined results and a dotplot in which the size of each dot reflects the enrichment, the color its significance, and the number written inside the dot the count of proteins found in that geneset.
#'
#' @param DEprot.timecourse.object An object of class \code{DEprot.timecourse}.
#' @param TERM2GENE Data.frame containing two columns 'gs_name' (IDs of the gene sets) and 'gene_symbol' (indicating the gene IDs). No default.
#' @param clusters Numeric value (or vector) indicating the cluster(s) to analyze. Default: \code{NULL} (all the clusters).
#' @param top.n Numeric value indicating how many best-ranked proteins of each cluster should be used for the enrichment. Default: \code{NULL} (all the proteins of the cluster).
#' @param universe Character vector indicating the background gene list. Default: \code{NULL}, meaning that all the proteins tested by \code{analyze.timecourse} are used. This is almost always the correct background for a proteomics experiment, since only the quantified proteins could have been detected as trending.
#' @param gsub.pattern.prot.id String indicating a pattern to be passed to gsub and to remove from the prot.id. Default: \code{NULL} (no changes in the IDs).
#' @param min.cluster.size Numeric value indicating the minimal number of proteins required to run the enrichment on a cluster. Default: \code{5}.
#' @param pvalueCutoff Numeric value indicating the adjusted pvalue cutoff on enrichment tests to report. Default: \code{0.05}.
#' @param qvalueCutoff Numeric value indicating the qvalue cutoff on enrichment tests to report as significant. Default: \code{0.05}.
#' @param pAdjustMethod String indicating the method to use for the p-value adjustment. One among "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: \code{"BH"}.
#' @param dotplot.n Numeric value indicating the maximum number of genesets displayed per cluster in the dotplot. Default: \code{5}.
#' @param size.by String indicating the metric mapped on the size of the dots. One among: 'FoldEnrichment' or 'GeneRatio'. Default: \code{"FoldEnrichment"}.
#' @param size.range Numeric vector of two elements indicating the minimal and maximal size of the dots. They must stay large enough for the numbers to fit inside. Default: \code{c(5, 14)}.
#' @param number.size Numeric value indicating the font size of the count displayed inside the dots. Default: \code{2.7}.
#' @param show.numbers Logic value indicating whether the number of proteins should be written inside the dots. Default: \code{TRUE}.
#' @param max.term.length Numeric value indicating the maximal number of characters of the geneset names displayed on the y-axis (longer names are truncated). Default: \code{60}.
#' @param viridis.option String indicating the viridis palette used for the significance color scale. Default: \code{"rocket"}.
#'
#' @return An object of class \code{DEprot.timecourse.enrichment}.
#'
#' @seealso \code{\link{analyze.timecourse}}, \code{\link{heatmap.timecourse}}, \code{\link{geneset.enrichment}}
#'
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#' @import viridis
#' @import clusterProfiler
#' @importFrom methods new
#' @importFrom stats na.omit
#' @importFrom utils head
#'
#' @author Sebastian Gregoricchio
#'
#' @export timecourse.enrichment

timecourse.enrichment =
  function(DEprot.timecourse.object,
           TERM2GENE,
           clusters = NULL,
           top.n = NULL,
           universe = NULL,
           gsub.pattern.prot.id = NULL,
           min.cluster.size = 5,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           dotplot.n = 5,
           size.by = "FoldEnrichment",
           size.range = c(5, 14),
           number.size = 2.7,
           show.numbers = TRUE,
           max.term.length = 60,
           viridis.option = "rocket") {

    ### Internal functions

    ## converts the 'x/y' strings returned by clusterProfiler into a numeric ratio
    parse.ratio =
      function(x) {
        vapply(strsplit(as.character(x), "/"),
               function(v) {as.numeric(v[1]) / as.numeric(v[2])},
               FUN.VALUE = numeric(1))
      }

    ######################################################################################

    ### check object
    if (!("DEprot.timecourse" %in% class(DEprot.timecourse.object))) {
      stop("The input must be an object of class 'DEprot.timecourse'.")
    }

    tc = DEprot.timecourse.object

    if (is.null(tc@clusters)) {
      stop("No clustering is available in this object: there is nothing to test.")
    }

    if (!(tolower(size.by) %in% c("foldenrichment", "generatio"))) {
      stop("The 'size.by' must be either 'FoldEnrichment' or 'GeneRatio'.")
    }

    size.by = ifelse(tolower(size.by) == "foldenrichment", "FoldEnrichment", "GeneRatio")


    ### define the clusters to test
    all.clusters = sort(unique(stats::na.omit(tc@results$cluster)))

    if (is.null(clusters)) {
      clusters = all.clusters
    } else if (!all(clusters %in% all.clusters)) {
      stop(paste0("Some of the 'clusters' requested are not available in this object.\n",
                  "       Available clusters: ", paste(all.clusters, collapse = ", ")))
    }


    ### define the background
    ## the correct universe is the list of the QUANTIFIED proteins: using the whole genome
    ## would make almost every geneset look enriched, since only the detected proteins had
    ## any chance of being called as trending
    if (is.null(universe)) {universe = rownames(tc@counts.used)}
    if (!is.null(gsub.pattern.prot.id)) {universe = gsub(gsub.pattern.prot.id, "", universe)}


    ### parameters used, collected here so that they are available at both the exit points
    enrichment.parameters = list(clusters = clusters,
                                 top.n = top.n,
                                 universe.size = length(universe),
                                 min.cluster.size = min.cluster.size,
                                 pvalueCutoff = pvalueCutoff,
                                 qvalueCutoff = qvalueCutoff,
                                 pAdjustMethod = pAdjustMethod,
                                 size.by = size.by)



    ######################################################################################
    ### Run the ORA cluster by cluster
    ######################################################################################

    enrichment.list = list()
    results.list = list()

    for (k in clusters) {

      ## the proteins of the cluster, optionally restricted to the best-ranked ones
      genes = get.timecourse.results(DEprot.timecourse.object = tc, cluster = k, top.n = top.n)$prot.id

      if (!is.null(gsub.pattern.prot.id)) {genes = gsub(gsub.pattern.prot.id, "", genes)}

      if (length(genes) < min.cluster.size) {
        warning(paste0("Cluster ", k, " contains only ", length(genes),
                       " protein(s): the enrichment has been skipped (see 'min.cluster.size')."))
        next
      }

      enrich = tryCatch(clusterProfiler::enricher(gene = genes,
                                                  universe = universe,
                                                  TERM2GENE = TERM2GENE,
                                                  pvalueCutoff = pvalueCutoff,
                                                  qvalueCutoff = qvalueCutoff,
                                                  pAdjustMethod = pAdjustMethod),
                        error = function(x) {return(NULL)})

      if (is.null(enrich)) {
        warning(paste0("The enrichment failed for the cluster ", k, "."))
        next
      }

      enrichment.list[[paste0("cluster.", k)]] = enrich

      tb = as.data.frame(enrich@result)

      if (nrow(tb) == 0) {next}

      ## the older versions of clusterProfiler do not return the FoldEnrichment column
      if (!("FoldEnrichment" %in% colnames(tb))) {
        tb$FoldEnrichment = parse.ratio(tb$GeneRatio) / parse.ratio(tb$BgRatio)
      }

      tb$GeneRatio.numeric = parse.ratio(tb$GeneRatio)
      tb$cluster = k
      tb$cluster.size = length(genes)

      results.list[[paste0("cluster.", k)]] = tb
    }


    ## early exit: an empty object is still returned, so that the downstream code does not
    ## have to handle a different class when nothing was enriched
    if (length(results.list) == 0) {
      warning("No enrichment could be computed for any of the clusters.")

      return(new(Class = "DEprot.timecourse.enrichment",
                 results = NULL,
                 enrichment.per.cluster = enrichment.list,
                 dotplot = NULL,
                 universe = universe,
                 parameters = enrichment.parameters))
    }

    results = dplyr::bind_rows(results.list)
    rownames(results) = NULL



    ######################################################################################
    ### Dotplot
    ######################################################################################

    ## only the significant genesets are displayed
    signif = dplyr::filter(.data = results, p.adjust <= pvalueCutoff)

    if (nrow(signif) == 0) {
      warning("No geneset passed the significance thresholds: the dotplot is not generated.")
      dotplot.tc = NULL

    } else {

      ## top genesets of EACH cluster, then their union: this keeps the plot readable while
      ## still showing, for every term, in which other clusters it also came up
      top.terms = unique(unlist(lapply(split(signif, signif$cluster),
                                       function(x) {utils::head(x[order(x$p.adjust),]$ID, dotplot.n)})))

      dot.tb = dplyr::filter(.data = signif, ID %in% top.terms)

      ## the y-axis is sorted by the cluster in which each term is the most significant,
      ## so that the dots line up along a diagonal following the temporal order
      term.order = dplyr::summarise(dplyr::group_by(dot.tb, ID),
                                    best.cluster = cluster[which.min(p.adjust)],
                                    best.padj = min(p.adjust),
                                    .groups = "drop")
      term.order = term.order[order(-term.order$best.cluster, term.order$best.padj),]

      dot.tb$ID = factor(dot.tb$ID, levels = term.order$ID)
      dot.tb$log10.padj = -log10(dot.tb$p.adjust)
      dot.tb$size.value = if (size.by == "FoldEnrichment") {dot.tb$FoldEnrichment} else {dot.tb$GeneRatio.numeric}

      ## the numbers must stay readable on both the dark and the light end of the palette
      mid.point = mean(range(dot.tb$log10.padj, na.rm = TRUE))
      dot.tb$label.color = ifelse(dot.tb$log10.padj >= mid.point, "white", "black")


      dotplot.tc =
        ggplot2::ggplot(data = dot.tb,
                        ggplot2::aes(x = factor(cluster), y = ID)) +
        ggplot2::geom_point(ggplot2::aes(size = size.value, fill = log10.padj),
                            shape = 21,
                            color = "grey30",
                            stroke = 0.3)

      ## the count of proteins is written in the middle of each dot
      if (show.numbers == TRUE) {
        dotplot.tc =
          dotplot.tc +
          ggplot2::geom_text(ggplot2::aes(label = Count, color = label.color),
                             size = number.size,
                             fontface = "bold") +
          ggplot2::scale_color_identity()
      }

      dotplot.tc =
        dotplot.tc +
        ggplot2::scale_size_continuous(range = size.range, name = size.by) +
        viridis::scale_fill_viridis(option = viridis.option, direction = -1, begin = 0.15, end = 0.9,
                                    name = "-log~10~(P~adj~)") +
        ggplot2::scale_y_discrete(labels = function(x) {
          ifelse(nchar(x) > max.term.length, paste0(substr(x, 1, max.term.length), "..."), x)
        }) +
        ggplot2::labs(x = "cluster",
                      y = NULL,
                      title = "**Over-representation analyses per cluster**",
                      caption = paste0("The number inside each dot is the count of proteins in the geneset (universe: ",
                                       length(universe), " proteins)")) +
        ggpubr::theme_pubr(legend = "right") +
        ggplot2::theme(axis.line = ggplot2::element_blank(),
                       panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1),
                       plot.title = ggtext::element_markdown(hjust = 0.5),
                       plot.caption = ggplot2::element_text(size = 7, color = "grey30"),
                       legend.title = ggtext::element_markdown(),
                       panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "gray80"),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    }



    ######################################################################################
    ### Build the output object
    ######################################################################################

    DEprot.timecourse.enrichment.object =
      new(Class = "DEprot.timecourse.enrichment",
          results = results,
          enrichment.per.cluster = enrichment.list,
          dotplot = dotplot.tc,
          universe = universe,
          parameters = enrichment.parameters)

    return(DEprot.timecourse.enrichment.object)
  } # END function











# ----------------------------------------------------------------------------------------
###                             TIMECOURSE DATA SIMULATOR                              ###
# ----------------------------------------------------------------------------------------

#' @title .simulate.timecourse
#'
#' @description
#' Generates a small synthetic time-course proteomics dataset, used by the vignette and by the unit tests. A fraction of the proteins follows one of five kinetic archetypes (monotone up, monotone down, early transient, late transient, early drop with rebound), the remaining ones being flat plus noise. A matching geneset annotation is produced as well, each geneset being enriched in the proteins of one archetype, so that the over-representation analyses have something to recover.
#'
#' @param n.proteins Numeric value indicating the number of proteins to simulate. Default: \code{1000}.
#' @param timepoints Numeric vector of the timepoints of the design. Default: \code{c(0, 1, 2, 6, 24)}.
#' @param n.replicates Numeric value indicating the number of replicates per timepoint. Default: \code{3}.
#' @param groups Character vector indicating the group levels to simulate (e.g., treatment and control). When more than one group is provided, only the first one shows the temporal response, the others staying flat. Default: \code{NULL} (a single group).
#' @param fraction.responsive Numeric value indicating the fraction of proteins showing a temporal trend. Default: \code{0.3}.
#' @param amplitude.range Numeric vector of two elements indicating the range of the response amplitudes, in log2 units. Default: \code{c(0.8, 3)}.
#' @param baseline.mean Numeric value indicating the mean log2 intensity of the proteins. Default: \code{22}.
#' @param baseline.sd Numeric value indicating the standard deviation of the baseline log2 intensities. Default: \code{2}.
#' @param noise.sd Numeric value indicating the standard deviation of the technical noise added to each measure. Default: \code{0.25}.
#' @param n.genesets.per.archetype Numeric value indicating how many genesets are generated for each kinetic archetype. Default: \code{2}.
#' @param geneset.size Numeric value indicating the number of members of each geneset. Default: \code{40}.
#' @param geneset.purity Numeric value between 0 and 1 indicating the fraction of the members of a geneset that are drawn from the corresponding archetype, the rest being drawn at random. Default: \code{0.6}.
#' @param na.fraction Numeric value indicating the fraction of measures to be replaced by \code{NA}. Default: \code{0}.
#' @param seed Numeric value indicating the seed used for the simulation. Default: \code{1234}.
#'
#' @return A list with four elements: \code{counts} (matrix of log2 intensities), \code{metadata} (data.frame with the columns column.id, time.hours, replicate, group and combined.id), \code{truth} (data.frame associating each protein to its archetype and amplitude) and \code{TERM2GENE} (data.frame with the columns gs_name and gene_symbol).
#'
#' @importFrom stats rnorm runif
#'
#' @keywords internal

.simulate.timecourse =
  function(n.proteins = 1000,
           timepoints = c(0, 1, 2, 6, 24),
           n.replicates = 3,
           groups = NULL,
           fraction.responsive = 0.3,
           amplitude.range = c(0.8, 3),
           baseline.mean = 22,
           baseline.sd = 2,
           noise.sd = 0.25,
           n.genesets.per.archetype = 2,
           geneset.size = 40,
           geneset.purity = 0.6,
           na.fraction = 0,
           seed = 1234) {

    set.seed(seed)

    if (is.null(groups)) {groups = "all"}


    ######################################################################################
    ### Kinetic archetypes
    ######################################################################################

    ## the shapes are defined on a log-spaced time, rescaled between 0 and 1: this is how
    ## a real degradation/induction experiment behaves, most of the action happening in
    ## the early timepoints even when the design runs up to 24 h
    u = log2(timepoints + 1)
    u = (u - min(u)) / (max(u) - min(u))

    archetypes = list(monotone.up      = u,
                      monotone.down    = -u,
                      early.transient  = sin(pi * u^0.6),
                      late.transient   = sin(pi * u^1.8),
                      early.drop       = -sin(pi * u^0.6))

    ## each shape is normalized so that the amplitude drawn below is the real log2 range
    archetypes = lapply(archetypes, function(s) {s / max(abs(s))})
    archetype.names = names(archetypes)



    ######################################################################################
    ### Assign the proteins to the archetypes
    ######################################################################################

    protein.id = paste0("protein.", seq_len(n.proteins))

    n.responsive = round(n.proteins * fraction.responsive)
    responsive = sample(protein.id, n.responsive)

    truth = data.frame(prot.id = protein.id,
                       archetype = "flat",
                       amplitude = 0,
                       stringsAsFactors = FALSE)

    ## the responsive proteins are split evenly across the five archetypes
    truth$archetype[match(responsive, truth$prot.id)] = rep(archetype.names, length.out = n.responsive)
    truth$amplitude[match(responsive, truth$prot.id)] = stats::runif(n = n.responsive,
                                                                     min = amplitude.range[1],
                                                                     max = amplitude.range[2])



    ######################################################################################
    ### Build the metadata
    ######################################################################################

    metadata = expand.grid(replicate = paste0("rep", seq_len(n.replicates)),
                           time.hours = timepoints,
                           group = groups,
                           stringsAsFactors = FALSE)

    ## 'if' and NOT 'ifelse': the latter returns a result of the same shape as the TEST,
    ## which is a scalar here, so only the first group name would be kept and recycled
    ## over all the samples, generating duplicated column IDs
    group.prefix = if (length(groups) > 1) {paste0(metadata$group, "_")} else {""}

    metadata$column.id = paste0(group.prefix, "T", metadata$time.hours, "_", metadata$replicate)
    metadata$combined.id = paste0(group.prefix, "T", metadata$time.hours)

    metadata = metadata[,c("column.id", "time.hours", "replicate", "group", "combined.id")]
    rownames(metadata) = NULL

    ## safety net: a duplicated ID would break the counts matrix and the DEprot object
    if (any(duplicated(metadata$column.id))) {
      stop("The simulated sample IDs are not unique: please check the 'groups' provided.")
    }



    ######################################################################################
    ### Generate the counts
    ######################################################################################

    baseline = stats::rnorm(n = n.proteins, mean = baseline.mean, sd = baseline.sd)
    names(baseline) = protein.id

    counts = matrix(NA_real_,
                    nrow = n.proteins,
                    ncol = nrow(metadata),
                    dimnames = list(protein.id, metadata$column.id))

    for (j in seq_len(nrow(metadata))) {
      i.time = which(timepoints == metadata$time.hours[j])

      ## the effect of the time is applied only to the first group: the others are used as
      ## controls, which is what makes the group x time interaction testable
      responds = (metadata$group[j] == groups[1])

      shift = vapply(seq_len(n.proteins),
                     function(i) {
                       if (truth$archetype[i] == "flat" | !responds) {return(0)}
                       return(truth$amplitude[i] * archetypes[[truth$archetype[i]]][i.time])
                     },
                     FUN.VALUE = numeric(1))

      counts[,j] = baseline + shift + stats::rnorm(n = n.proteins, mean = 0, sd = noise.sd)
    }

    ## optional missing values, to mimic an unimputed table
    if (na.fraction > 0) {
      n.na = round(length(counts) * na.fraction)
      counts[sample(seq_along(counts), n.na)] = NA
    }



    ######################################################################################
    ### Generate a matching geneset annotation
    ######################################################################################

    TERM2GENE =
      do.call(rbind,
              lapply(archetype.names,
                     function(a) {
                       members.pool = truth$prot.id[truth$archetype == a]
                       others.pool = truth$prot.id[truth$archetype != a]

                       do.call(rbind,
                               lapply(seq_len(n.genesets.per.archetype),
                                      function(s) {
                                        ## a geneset is only partially pure, otherwise the
                                        ## enrichment would be unrealistically clean
                                        n.pure = min(length(members.pool), round(geneset.size * geneset.purity))
                                        n.rest = geneset.size - n.pure

                                        genes = c(sample(members.pool, n.pure),
                                                  sample(others.pool, n.rest))

                                        data.frame(gs_name = paste0(gsub("\\.", "_", toupper(a)), "_SET_", s),
                                                   gene_symbol = genes,
                                                   stringsAsFactors = FALSE)
                                      }))
                     }))

    ## a few genesets unrelated to any archetype, as negative controls
    TERM2GENE =
      rbind(TERM2GENE,
            do.call(rbind,
                    lapply(seq_len(3),
                           function(s) {
                             data.frame(gs_name = paste0("BACKGROUND_SET_", s),
                                        gene_symbol = sample(protein.id, geneset.size),
                                        stringsAsFactors = FALSE)
                           })))

    rownames(TERM2GENE) = NULL



    return(list(counts = counts,
                metadata = metadata,
                truth = truth,
                TERM2GENE = TERM2GENE))
  } # END function
