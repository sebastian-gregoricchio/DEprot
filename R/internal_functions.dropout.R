##########################
### INTERNAL FUNCTIONS ###
###      dropout       ###
##########################

# ----------------------------------------------------------------------------------------

#' @title estimate.dropout
#'
#' @description
#' Estimates one logistic dropout curve per set of columns (a group of replicates or a single sample).
#' The probability that a value is missing is modelled as a function of the protein abundance,
#' P(missing) = plogis(intercept + slope * abundance), where the abundance is the mean of the values
#' measured for that protein across the whole table. A protein entirely absent from a group still
#' contributes to the fit of that group, which is what carries the information about the left tail.
#'
#' @param counts Numeric matrix of counts (proteins x samples), NAs included.
#' @param column.sets Named list of column indexes (or column names), one element per group or per sample.
#' @param min.missing Minimum number of missing values required in a set to attempt the fit. Default: \code{20}.
#'
#' @return Named list, one element per column set, each a numeric vector with \code{intercept} and \code{slope}
#'   (both \code{NA} when the curve could not be estimated).
#'
#' @importFrom stats glm binomial coef
#'
#' @keywords internal

.estimate.dropout =
  function(counts,
           column.sets,
           min.missing = 20) {

    ## protein abundance proxy, shared by all the sets:
    ## the mean of the values actually measured, which is defined only for the proteins detected at least once
    abundance = rowMeans(counts, na.rm = TRUE)
    usable = is.finite(abundance)

    fits = list()

    for (i in 1:length(column.sets)) {
      set.counts = counts[usable, column.sets[[i]], drop = FALSE]
      n.missing = rowSums(is.na(set.counts))
      n.observed = ncol(set.counts) - n.missing

      ## a curve requires missing values, and both fully detected and partially missing proteins
      if (sum(n.missing) < min.missing | length(unique(n.missing)) < 2) {
        fits[[i]] = c(intercept = NA, slope = NA)
        next
      }

      fit = tryCatch(suppressWarnings(stats::glm(cbind(n.missing, n.observed) ~ abundance[usable],
                                                 family = stats::binomial(link = "logit"))),
                     error = function(e) {return(NULL)})

      if (is.null(fit)) {
        fits[[i]] = c(intercept = NA, slope = NA)
      } else {
        fits[[i]] = c(intercept = unname(stats::coef(fit)[1]),
                      slope = unname(stats::coef(fit)[2]))
      }
    }

    names(fits) = names(column.sets)
    return(fits)
  }




# ----------------------------------------------------------------------------------------

#' @title dropout.weights
#'
#' @description
#' Builds the matrix of masking weights used to choose which values are hidden in the test dataset.
#' Each column receives the curve estimated for the set it belongs to, hence samples with different
#' detection sensitivity (an IgG control and an IP, for instance) are masked at different rates and in
#' different regions of the intensity range. The cells that are already missing get a null weight: they
#' cannot be masked because there would be no measured value to compare the imputation with.
#'
#' @param test.counts Numeric matrix of the test dataset (proteins x samples).
#' @param abundance Numeric vector of the abundance proxy of the proteins of \code{test.counts}, in the same order.
#' @param fits Named list of dropout curves, as returned by \code{\link{.estimate.dropout}}.
#' @param column.sets Named list of column indexes, with the same names as \code{fits}.
#' @param flat Logical value, when \code{TRUE} all the measured cells receive the same weight (MCAR masking). Default: \code{FALSE}.
#'
#' @return Numeric matrix with the same dimensions as \code{test.counts}.
#'
#' @keywords internal

.dropout.weights =
  function(test.counts,
           abundance,
           fits,
           column.sets,
           flat = FALSE) {

    weights = matrix(1,
                     nrow = nrow(test.counts),
                     ncol = ncol(test.counts),
                     dimnames = dimnames(test.counts))

    if (isFALSE(flat)) {
      for (i in 1:length(column.sets)) {
        fit = fits[[names(column.sets)[i]]]

        ## a curve that could not be estimated, or that increases with the abundance, is not used:
        ## the corresponding samples keep a flat weight
        if (any(is.na(fit)) | isTRUE(fit["slope"] >= 0)) {
          next
        }

        w = 1 / (1 + exp(-(fit["intercept"] + fit["slope"] * abundance)))
        w[!is.finite(w)] = 0
        weights[, column.sets[[i]]] = w
      }
    }

    ## the values already missing cannot be used as ground truth
    weights[is.na(test.counts)] = 0

    return(weights)
  }




# ----------------------------------------------------------------------------------------

#' @title sample.idx
#'
#' @description
#' Wrapper of \code{sample} that returns an empty vector instead of an error when there is nothing to
#' draw, and that does not interpret a single value as the upper limit of a sequence.
#'
#' @param x Vector of indexes to sample from.
#' @param size Number of elements to draw.
#' @param prob Vector of weights, of the same length as \code{x}. Default: \code{NULL}, uniform.
#'
#' @return Vector of indexes, of length lower than or equal to \code{size}.
#'
#' @keywords internal

.sample.idx =
  function(x,
           size,
           prob = NULL) {

    if (length(x) == 0 | size <= 0) {
      return(c())
    }

    if (!is.null(prob)) {
      prob = prob[1:length(x)]
      ## indexes that cannot be drawn are removed rather than given a null weight,
      ## otherwise 'sample' fails when the requested size exceeds the number of usable elements
      x = x[prob > 0]
      prob = prob[prob > 0]
    }

    size = min(size, length(x))

    if (size == 0) {
      return(c())
    } else if (length(x) == 1) {
      return(x)
    } else {
      return(sample(x = x, size = size, replace = FALSE, prob = prob))
    }
  }
