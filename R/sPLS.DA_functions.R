# ----------------------------------------------------------------------------------------
###                          sPLS-DA INTERNAL FUNCTIONS                                 ###
# ----------------------------------------------------------------------------------------

#' @title .splsda.prepare
#'
#' @description
#' Internal. Extracts the counts, converts them to log2, applies the sample subset and removes the proteins that cannot enter the model. Returns the sample-by-protein matrix expected by mixOmics.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param which.data String indicating which type of counts should be used.
#' @param sample.subset String vector indicating the samples to keep, or \code{NULL}.
#' @param scale.data Logical value indicating whether the data will be scaled by mixOmics.
#'
#' @return A list with the elements: \code{mat} (proteins x samples, log2), \code{X} (samples x proteins), \code{metadata} and \code{data.used}.
#'
#' @import dplyr
#' @importFrom stats var
#'
#' @keywords internal

.splsda.prepare =
  function(DEprot.object,
           which.data,
           sample.subset = NULL,
           scale.data = TRUE) {

    ### Extract the requested counts table (aliases handled by the shared helper)
    counts = .tc.get.counts(DEprot.object = DEprot.object, which.data = which.data)
    mat = counts$mat


    ## sPLS-DA works on the intensities directly: on a linear scale the few most abundant
    ## proteins would dominate both the deflation and the variable selection
    if (!is.numeric(DEprot.object@log.base)) {
      message("The log.base is not numeric, linear counts are assumed. Counts matrix will be converted to log2(score+1) values to analyze the data.")
      mat = log2(mat + 1)
    } else if (as.numeric(DEprot.object@log.base) != 2) {
      message("The log.base is not 2, counts will be converted to log2 values to analyze the data.")
      mat = log2(DEprot.object@log.base^mat)
    }


    ### Subset the samples
    if (!is.null(sample.subset)) {
      missing.samples = setdiff(sample.subset, colnames(mat))

      if (length(missing.samples) > 0) {
        stop("The following sample(s) indicated in the 'sample.subset' are not available in the counts table:\n",
             "       ", paste(missing.samples, collapse = ", "))
      }

      mat = mat[, colnames(mat) %in% sample.subset, drop = FALSE]
    }

    meta = DEprot.object@metadata
    meta = meta[match(colnames(mat), meta$column.id), , drop = FALSE]
    rownames(meta) = NULL

    if (ncol(mat) < 4) {
      stop("Only ", ncol(mat), " sample(s) available: a sPLS-DA cannot be cross-validated on such a design.")
    }


    ### Clean the matrix
    mat[is.nan(mat)] = NA

    ## Unlike a PCA, which can fall back on NIPALS, both the tuning and the performance
    ## estimation resample the samples: any missing value would be silently handled in a
    ## different way at every fold, and the error rates would not be comparable
    if (TRUE %in% is.na(mat)) {
      n.na = sum(is.na(mat))
      stop("The counts table contains ", n.na, " missing value(s), which sPLS-DA cannot handle.\n",
           "       Use 'which.data = \"imputed\"' or impute the data with `impute.counts()`.")
    }

    ## A protein with the same value in every sample carries no information and cannot be
    ## rescaled to unit variance. This happens easily on a subset of samples, so these
    ## proteins are dropped rather than making the whole analysis fail
    row.variance = apply(X = mat, MARGIN = 1, FUN = function(x) {stats::var(x, na.rm = TRUE)})
    constant.rows = which(is.na(row.variance) | row.variance == 0)

    if (length(constant.rows) > 0) {
      warning("Data contain ", length(constant.rows), " row(s) (proteins) with a null variance across the samples analyzed. ",
              "These rows will be removed to perform the sPLS-DA.")
      mat = mat[-constant.rows, , drop = FALSE]
    }

    if (nrow(mat) < 2) {
      stop("Only ", nrow(mat), " protein(s) are usable: nothing can be modeled.")
    }


    return(list(mat = mat,
                X = t(mat),
                metadata = meta,
                data.used = counts$data.used))
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .splsda.groups
#'
#' @description
#' Internal. Builds and validates the class vector used as response of the sPLS-DA.
#'
#' @param metadata data.frame of the metadata of the samples analyzed.
#' @param group.column String indicating the metadata column holding the classes.
#' @param reference.group String indicating the class used to orient the components, or \code{NULL}.
#'
#' @return A list with the elements \code{Y} (factor) and \code{reference.group}.
#'
#' @keywords internal

.splsda.groups =
  function(metadata,
           group.column,
           reference.group = NULL) {

    if (!(group.column %in% colnames(metadata))) {
      stop("The 'group.column' ('", group.column, "') is not present in the metadata table.\n",
           "       Available columns: ", paste(colnames(metadata), collapse = ", "))
    }

    Y = metadata[, group.column]

    if (TRUE %in% is.na(Y)) {
      stop("The 'group.column' ('", group.column, "') contains missing values: every sample must be assigned to a class.")
    }

    Y = droplevels(factor(Y))

    if (nlevels(Y) < 2) {
      stop("The 'group.column' ('", group.column, "') defines a single class: there is nothing to discriminate.")
    }

    group.sizes = table(Y)

    if (min(group.sizes) < 2) {
      stop("The class(es) ", paste(names(group.sizes)[group.sizes < 2], collapse = ", "),
           " contain a single sample: a discriminant model cannot be cross-validated.")
    }

    if (min(group.sizes) < 3) {
      warning("The smallest class contains ", min(group.sizes), " samples only. ",
              "The error rates estimated on such a design are very unstable and should be read as indicative.")
    }


    ### Reference group used to orient the components
    if (is.null(reference.group)) {
      reference.group = levels(Y)[1]
    } else if (!(reference.group %in% levels(Y))) {
      stop("The 'reference.group' ('", reference.group, "') is not one of the classes of '", group.column, "'.\n",
           "       Available classes: ", paste(levels(Y), collapse = ", "))
    }


    return(list(Y = Y,
                reference.group = reference.group))
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .splsda.validation
#'
#' @description
#' Internal. Adapts the cross-validation scheme to the size of the smallest class. Asking for more folds than the samples available in a class produces empty folds, and mixOmics would either fail or silently return an over-optimistic error rate.
#'
#' @param Y Factor of the classes.
#' @param validation String indicating the validation scheme ('Mfold' or 'loo').
#' @param folds Number of folds.
#' @param nrepeat Number of repeats of the cross-validation.
#'
#' @return A list with the elements \code{validation}, \code{folds} and \code{nrepeat}.
#'
#' @keywords internal

.splsda.validation =
  function(Y,
           validation = "Mfold",
           folds = 5,
           nrepeat = 10) {

    min.size = min(table(Y))

    if (tolower(validation) %in% c("loo", "leave-one-out", "leaveoneout")) {
      ## leave-one-out is deterministic: repeating it would return the same partition
      return(list(validation = "loo",
                  folds = length(Y),
                  nrepeat = 1))
    }

    if (min.size < 3) {
      warning("The smallest class contains ", min.size, " samples: the cross-validation has been switched to leave-one-out.")
      return(list(validation = "loo",
                  folds = length(Y),
                  nrepeat = 1))
    }

    if (folds > min.size) {
      warning("The number of 'folds' (", folds, ") is larger than the smallest class (", min.size,
              " samples): 'folds' has been set to ", min.size, ".")
      folds = min.size
    }

    return(list(validation = "Mfold",
                folds = round(folds, 0),
                nrepeat = round(nrepeat, 0)))
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .splsda.parallel
#'
#' @description
#' Internal. Calls a mixOmics function requesting several cores, falling back to a serial run when the parallel back-end is not available. The name of the argument controlling the parallelization changed along the mixOmics versions, so a failure of the parallel call is not by itself informative.
#'
#' @param FUN The mixOmics function to call.
#' @param args Named list of arguments passed to \code{FUN}.
#' @param n.cores Number of cores requested.
#'
#' @return The output of \code{FUN}.
#'
#' @keywords internal

.splsda.parallel =
  function(FUN,
           args,
           n.cores = 1) {

    result = NULL

    if (n.cores > 1) {
      result = tryCatch(do.call(FUN, c(args, list(cpus = n.cores))),
                        error = function(e) {return(NULL)})

      if (is.null(result)) {
        message("The parallel back-end is not available for this version of mixOmics: the computation will run on a single core.")
      }
    }

    if (is.null(result)) {
      result = do.call(FUN, args)
    }

    return(result)
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .splsda.sign
#'
#' @description
#' Internal. Computes the sign to apply to each component so that the reference class always sits on the positive side.
#'
#' @param variates Matrix of the sample scores (samples x components).
#' @param Y Factor of the classes.
#' @param reference.group String indicating the class placed on the positive side.
#'
#' @return A numeric vector of 1 and -1, one value per component.
#'
#' @keywords internal

.splsda.sign =
  function(variates,
           Y,
           reference.group) {

    ## mixOmics returns an arbitrary orientation for each component. Left as it is, the sign of
    ## the loadings (and therefore the direction of any NES computed from them) would change from
    ## one run to the next, which makes the results impossible to compare
    group.means = apply(X = variates,
                        MARGIN = 2,
                        FUN = function(x) {mean(x[Y == reference.group], na.rm = TRUE)})

    signs = ifelse(group.means < 0, yes = -1, no = 1)
    signs[is.na(signs)] = 1

    return(as.numeric(signs))
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .splsda.expl.var
#'
#' @description
#' Internal. Extracts the proportion of variance explained by each component, for both the protein block (X) and the class block (Y).
#'
#' @param model A mixOmics model of class \code{mixo_splsda} or \code{mixo_plsda}.
#' @param ncomp Number of components.
#'
#' @return A data.frame with one row per component.
#'
#' @keywords internal

.splsda.expl.var =
  function(model,
           ncomp) {

    ## the slot holding the explained variance was renamed along the mixOmics 6.x versions
    if (!is.null(model$prop_expl_var)) {
      ev = model$prop_expl_var
    } else if (!is.null(model$explained_variance)) {
      ev = model$explained_variance
    } else {
      ev = list(X = rep(NA_real_, ncomp), Y = rep(NA_real_, ncomp))
    }

    var.x = as.numeric(ev$X)[seq_len(ncomp)]
    var.y = as.numeric(ev$Y)[seq_len(ncomp)]

    importance =
      data.frame(comp = factor(seq_len(ncomp), levels = seq_len(ncomp)),
                 Proportion.of.Variance = var.x,
                 Cumulative.Proportion = cumsum(ifelse(is.na(var.x), 0, var.x)),
                 Percentage.of.Variance = var.x * 100,
                 Proportion.of.Variance.Y = var.y,
                 Cumulative.Proportion.Y = cumsum(ifelse(is.na(var.y), 0, var.y)),
                 Percentage.of.Variance.Y = var.y * 100)

    return(importance)
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .splsda.contribution
#'
#' @description
#' Internal. Assigns each protein to the class showing its highest median intensity. This is the same rule applied by \code{mixOmics::plotLoadings} and it is what makes a loading readable: a coefficient alone says how much a protein weighs on a component, not in which group it is high.
#'
#' @param mat Counts matrix used for the model (proteins x samples, log2).
#' @param Y Factor of the classes, aligned to the columns of \code{mat}.
#' @param proteins Character vector of the proteins to annotate.
#'
#' @return A character vector of class names, one per protein.
#'
#' @importFrom stats median
#'
#' @keywords internal

.splsda.contribution =
  function(mat,
           Y,
           proteins) {

    if (length(proteins) == 0) {return(character(0))}

    medians = vapply(X = levels(Y),
                     FUN = function(g) {apply(X = mat[proteins, Y == g, drop = FALSE],
                                              MARGIN = 1,
                                              FUN = stats::median, na.rm = TRUE)},
                     FUN.VALUE = numeric(length(proteins)))

    ## vapply drops the matrix shape when a single protein is annotated
    medians = matrix(medians,
                     nrow = length(proteins),
                     dimnames = list(proteins, levels(Y)))

    return(colnames(medians)[max.col(medians, ties.method = "first")])
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .splsda.loadings
#'
#' @description
#' Internal. Reshapes the loading matrix of a mixOmics model into a long table, applying the sign correction and annotating the class in which each protein is the highest.
#'
#' @param model A mixOmics model of class \code{mixo_splsda} or \code{mixo_plsda}.
#' @param ncomp Number of components.
#' @param signs Numeric vector of the sign corrections, one per component.
#' @param mat Counts matrix used for the model (proteins x samples, log2).
#' @param Y Factor of the classes.
#' @param sparse Logical value indicating whether the model is sparse (only the selected proteins are annotated).
#'
#' @return A data.frame with the columns \code{prot.id}, \code{component}, \code{loading}, \code{abs.loading}, \code{selected} and \code{contrib.group}.
#'
#' @keywords internal

.splsda.loadings =
  function(model,
           ncomp,
           signs,
           mat,
           Y,
           sparse = TRUE) {

    loadings.mat = as.matrix(model$loadings$X)

    tb.list = list()

    for (k in seq_len(ncomp)) {
      loading = as.numeric(loadings.mat[, k]) * signs[k]

      tb = data.frame(prot.id = rownames(loadings.mat),
                      component = k,
                      loading = loading,
                      abs.loading = abs(loading),
                      selected = loading != 0,
                      stringsAsFactors = FALSE)

      ## for the dense model every protein has a non-null coefficient, so there is no
      ## selection to report and annotating all of them would be pointlessly slow
      if (sparse == TRUE) {
        tb$contrib.group = NA_character_
        selected.proteins = tb$prot.id[tb$selected]
        tb$contrib.group[tb$selected] = .splsda.contribution(mat = mat, Y = Y, proteins = selected.proteins)
      } else {
        tb$selected = NA
        tb$contrib.group = NA_character_
      }

      tb.list[[k]] = tb[order(-tb$abs.loading), ]
    }

    loadings = do.call(rbind, tb.list)
    rownames(loadings) = NULL

    return(loadings)
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .splsda.stability
#'
#' @description
#' Internal. Collects the selection frequency of each protein across the cross-validation folds. A protein selected in every fold is a robust marker; one selected in a third of them is an artifact of the particular split.
#'
#' @param perf.result Output of \code{mixOmics::perf}.
#' @param ncomp Number of components.
#'
#' @return A data.frame with the columns \code{prot.id}, \code{component} and \code{frequency}, or \code{NULL}.
#'
#' @import dplyr
#'
#' @keywords internal

.splsda.stability =
  function(perf.result,
           ncomp) {

    stable = tryCatch(perf.result$features$stable, error = function(e) {return(NULL)})

    if (is.null(stable)) {return(NULL)}

    ## depending on the version and on 'nrepeat', mixOmics returns either one table per
    ## component, or one list per repeat each holding one table per component
    if (is.list(stable[[1]])) {
      per.component = lapply(X = seq_len(ncomp),
                             FUN = function(k) {lapply(stable, function(x) {x[[k]]})})
    } else {
      per.component = lapply(X = seq_len(ncomp), FUN = function(k) {list(stable[[k]])})
    }


    tb.list = list()

    for (k in seq_len(ncomp)) {
      tables = per.component[[k]]
      tables = tables[!vapply(tables, is.null, FUN.VALUE = logical(1))]

      if (length(tables) == 0) {next}

      long = do.call(rbind,
                     lapply(X = tables,
                            FUN = function(x) {data.frame(prot.id = names(x),
                                                          frequency = as.numeric(x),
                                                          stringsAsFactors = FALSE)}))

      ## averaged over the repeats, so that the frequency stays a proportion of folds
      long = dplyr::summarise(dplyr::group_by(long, prot.id),
                              frequency = mean(frequency, na.rm = TRUE),
                              .groups = "drop")

      long = as.data.frame(long)
      long$component = k

      tb.list[[k]] = long[, c("prot.id", "component", "frequency")]
    }

    if (length(tb.list) == 0) {return(NULL)}

    stability = do.call(rbind, tb.list)
    stability = stability[order(stability$component, -stability$frequency), ]
    rownames(stability) = NULL

    return(stability)
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .splsda.auc
#'
#' @description
#' Internal. Reshapes the cross-validated AUC returned by \code{mixOmics::perf} into a long table.
#'
#' @param perf.result Output of \code{mixOmics::perf}.
#'
#' @return A data.frame with the columns \code{component}, \code{comparison}, \code{AUC} and \code{AUC.sd}, or \code{NULL}.
#'
#' @keywords internal

.splsda.auc =
  function(perf.result) {

    auc = tryCatch(perf.result$auc, error = function(e) {return(NULL)})

    if (is.null(auc)) {return(NULL)}

    tb.list = list()

    for (k in seq_along(auc)) {
      m = tryCatch(as.data.frame(auc[[k]]), error = function(e) {return(NULL)})

      if (is.null(m)) {next}
      if (nrow(m) == 0) {next}

      sd.column = grep("sd", colnames(m), ignore.case = TRUE)
      auc.column = setdiff(grep("auc", colnames(m), ignore.case = TRUE), sd.column)

      if (length(auc.column) == 0) {next}

      tb.list[[k]] =
        data.frame(component = k,
                   comparison = rownames(m),
                   AUC = as.numeric(m[, auc.column[1]]),
                   AUC.sd = if (length(sd.column) > 0) {as.numeric(m[, sd.column[1]])} else {NA_real_},
                   stringsAsFactors = FALSE)
    }

    if (length(tb.list) == 0) {return(NULL)}

    auc.tb = do.call(rbind, tb.list)
    rownames(auc.tb) = NULL

    return(auc.tb)
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .splsda.error.table
#'
#' @description
#' Internal. Reshapes the classification error rates returned by \code{mixOmics::perf} into a long table.
#'
#' @param perf.result Output of \code{mixOmics::perf}.
#'
#' @return A data.frame with the columns \code{component}, \code{distance}, \code{metric} and \code{error.rate}, or \code{NULL}.
#'
#' @keywords internal

.splsda.error.table =
  function(perf.result) {

    error.rate = tryCatch(perf.result$error.rate, error = function(e) {return(NULL)})

    if (is.null(error.rate)) {return(NULL)}

    melt.matrix =
      function(m, metric) {
        if (is.null(m)) {return(NULL)}
        m = as.matrix(m)
        data.frame(component = rep(seq_len(nrow(m)), times = ncol(m)),
                   distance = rep(colnames(m), each = nrow(m)),
                   metric = metric,
                   error.rate = as.numeric(m),
                   stringsAsFactors = FALSE)
      }

    tb = rbind(melt.matrix(error.rate$overall, "Overall error rate"),
               melt.matrix(error.rate$BER, "Balanced error rate"))

    if (is.null(tb)) {return(NULL)}

    rownames(tb) = NULL

    return(tb)
  } # END function




# ----------------------------------------------------------------------------------------
###                                    MAIN FUNCTION                                    ###
# ----------------------------------------------------------------------------------------


#' @title perform.sPLSDA
#'
#' @description Performs a sparse Partial Least Squares Discriminant Analysis (sPLS-DA) on a \code{DEprot} object. Differently from a PCA, which looks for the directions carrying the most variance whatever their origin, sPLS-DA looks for the directions that separate a set of classes given by the user, and selects at the same time the small group of proteins responsible for that separation.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param group.column String indicating the column of the metadata table holding the classes to discriminate (e.g., \code{"condition"}). No default.
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'randomized', 'random', 'imputed', 'imp'. Default: \code{"imputed"}.
#' @param sample.subset String vector indicating the column names (samples) to keep in the counts table (the 'column.id' in the metadata table). Default: \code{NULL} (no subsetting).
#' @param ncomp Integer number indicating the number of components to compute. Values lower than 3 are raised to 3, since the third component is what tells whether the separation seen on the first two is the whole story. Default: \code{3}.
#' @param keepX Number of proteins to retain on each component. A single value is recycled over all the components, a vector must have one value per component. When \code{NULL} (default) the value is tuned by cross-validation over \code{test.keepX}.
#' @param test.keepX Numeric vector of the values of \code{keepX} tested during the tuning. Used only when \code{keepX = NULL}. Default: \code{NULL} (\code{c(5, 10, 15, 20, 30, 50, 100)}, capped to the number of proteins available).
#' @param scale.data Logical value indicating whether each protein should be scaled to unit variance. Default: \code{TRUE}.
#' @param reference.group String indicating the class placed on the positive side of every component. Default: \code{NULL} (the first level of \code{group.column}).
#' @param validate Logical value indicating whether the performance of the model should be estimated by cross-validation. This is a second resampling run, roughly as expensive as the tuning, so it can be skipped while exploring. Default: \code{TRUE}.
#' @param validation String indicating the cross-validation scheme. One among: 'Mfold', 'loo'. Default: \code{"Mfold"}.
#' @param folds Number of folds of the cross-validation. Automatically capped to the size of the smallest class. Default: \code{5}.
#' @param nrepeat Number of times the cross-validation is repeated. Ignored for leave-one-out. Default: \code{10}.
#' @param dist String indicating the prediction distance used during the tuning. One among: 'max.dist', 'centroids.dist', 'mahalanobis.dist'. All three are reported by the validation. Default: \code{"max.dist"}.
#' @param n.cores Number of cores used for the tuning and the validation. Default: \code{1}.
#' @param seed Numeric value used to seed the resampling, so that the tuning and the error rates are reproducible. Default: \code{1234}.
#'
#' @details
#' Both the tuning of \code{keepX} and the estimation of the performances happen inside this function: neither of them changes as long as the model does not change, so keeping them apart would only mean carrying around three objects that must be kept in sync by hand.
#'
#' Alongside the sparse model, a plain (non-sparse) PLS-DA is fitted on the same data and the same number of components. It costs nothing, since no resampling is involved, and it gives a complete loading vector over every protein of the matrix. That vector is what \link{get.sPLSDA.ranking} returns for a GSEA: the sparse loadings are zero for all the proteins that were not selected, which leaves a ranking made mostly of ties and an enrichment score that depends on the order in which those ties happen to be sorted.
#'
#' The orientation of the components returned by mixOmics is arbitrary. Here the sign of each component is fixed so that \code{reference.group} always sits on the positive side, which also means that a protein with a positive loading is a protein higher in that group. The mixOmics object stored in the \code{splsda} slot keeps its original orientation, so that anything computed on it directly stays consistent with the mixOmics documentation.
#'
#' @return A \code{DEprot.sPLSDA} object.
#'
#' @seealso \link{plot.sPLSDA.scatter}, \link{plot.sPLSDA.loadings}, \link{get.sPLSDA.proteins}, \link{sPLSDA.enrichment}, \link{perform.PCA}
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom mixOmics splsda plsda perf vip tune.splsda
#' @importFrom methods new
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' # Model with a fixed number of proteins per component (no tuning)
#' splsda <- perform.sPLSDA(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                          group.column = "condition",
#'                          keepX = c(10, 10, 5),
#'                          validate = FALSE)
#'
#' \donttest{
#' # Tuned model, cross-validated (both steps resample the samples)
#' splsda.tuned <- perform.sPLSDA(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                                group.column = "condition",
#'                                test.keepX = c(5, 10, 20),
#'                                folds = 3,
#'                                nrepeat = 1)
#' }
#'
#' @export perform.sPLSDA

perform.sPLSDA =
  function(DEprot.object,
           group.column,
           which.data = "imputed",
           sample.subset = NULL,
           ncomp = 3,
           keepX = NULL,
           test.keepX = NULL,
           scale.data = TRUE,
           reference.group = NULL,
           validate = TRUE,
           validation = "Mfold",
           folds = 5,
           nrepeat = 10,
           dist = "max.dist",
           n.cores = 1,
           seed = 1234) {

    ### check object
    if (!methods::is(DEprot.object, "DEprot")) {
      stop("The input must be an object of class 'DEprot'.")
    }


    ### Prepare the counts and the classes
    prepared = .splsda.prepare(DEprot.object = DEprot.object,
                               which.data = which.data,
                               sample.subset = sample.subset,
                               scale.data = scale.data)

    mat = prepared$mat
    X = prepared$X

    groups = .splsda.groups(metadata = prepared$metadata,
                            group.column = group.column,
                            reference.group = reference.group)

    Y = groups$Y
    reference.group = groups$reference.group


    ### Number of components
    ## A two-group design does not stop at one component: the following ones keep being
    ## extracted from the deflated matrix, and they are exactly what shows whether the
    ## separation is driven by one protein module or by several
    ncomp = round(ncomp, 0)

    if (ncomp < 3) {
      message("The 'ncomp' has been raised to 3: the third component is what tells whether the separation seen on the first two is the whole story.")
      ncomp = 3
    }

    max.ncomp = min(c(nrow(X) - 1, ncol(X)))

    if (ncomp > max.ncomp) {
      warning("The 'ncomp' requested (", ncomp, ") exceeds what this design can support: it has been set to ", max.ncomp, ".")
      ncomp = max.ncomp
    }

    if (ncomp < 2) {
      stop("This design supports ", ncomp, " component(s) only: a sPLS-DA cannot be computed.")
    }


    ### Cross-validation scheme, adapted to the smallest class
    cv = .splsda.validation(Y = Y,
                            validation = validation,
                            folds = folds,
                            nrepeat = nrepeat)



    ######################################################################################
    ### Tuning of keepX
    ######################################################################################

    tuning = NULL

    if (is.null(keepX)) {

      if (is.null(test.keepX)) {
        test.keepX = c(5, 10, 15, 20, 30, 50, 100)
      }

      test.keepX = sort(unique(pmin(round(test.keepX, 0), nrow(mat))))
      test.keepX = test.keepX[test.keepX > 0]

      if (length(test.keepX) < 2) {
        stop("The 'test.keepX' grid must contain at least two distinct values not exceeding the number of proteins (", nrow(mat), ").")
      }

      message("Tuning 'keepX' by cross-validation over ", length(test.keepX), " values...")

      set.seed(seed)

      tune.result =
        .splsda.parallel(FUN = mixOmics::tune.splsda,
                         args = list(X = X,
                                     Y = Y,
                                     ncomp = ncomp,
                                     test.keepX = test.keepX,
                                     validation = cv$validation,
                                     folds = cv$folds,
                                     nrepeat = cv$nrepeat,
                                     dist = dist,
                                     scale = scale.data,
                                     progressBar = FALSE),
                         n.cores = n.cores)

      keepX = as.numeric(tune.result$choice.keepX)[seq_len(ncomp)]

      ## the tuning grid is kept in a long format, ready to be plotted
      error.matrix = as.matrix(tune.result$error.rate)

      tuning.tb =
        data.frame(keepX = as.numeric(rep(rownames(error.matrix), times = ncol(error.matrix))),
                   component = rep(seq_len(ncol(error.matrix)), each = nrow(error.matrix)),
                   error.rate = as.numeric(error.matrix),
                   stringsAsFactors = FALSE)

      error.sd = tryCatch(as.numeric(as.matrix(tune.result$error.rate.sd)), error = function(e) {return(NULL)})

      if (!is.null(error.sd) & length(error.sd) == nrow(tuning.tb)) {
        tuning.tb$error.rate.sd = error.sd
      } else {
        tuning.tb$error.rate.sd = NA_real_
      }

      tuning = list(test.keepX = test.keepX,
                    choice.keepX = keepX,
                    error.rate = tuning.tb,
                    dist = dist,
                    tune.object = tune.result)

    } else {
      ### keepX provided by the user
      keepX = round(keepX, 0)

      if (length(keepX) == 1) {
        keepX = rep(keepX, ncomp)
      } else if (length(keepX) != ncomp) {
        stop("The 'keepX' must be a single value or a vector of ", ncomp, " values (one per component), but ",
             length(keepX), " were provided.")
      }

      if (TRUE %in% (keepX > nrow(mat))) {
        warning("Some 'keepX' values exceed the number of proteins available (", nrow(mat), "): they have been capped.")
        keepX = pmin(keepX, nrow(mat))
      }

      if (TRUE %in% (keepX < 1)) {
        stop("The 'keepX' values must be positive.")
      }
    }



    ######################################################################################
    ### Sparse and dense models
    ######################################################################################

    set.seed(seed)

    model = mixOmics::splsda(X = X,
                             Y = Y,
                             ncomp = ncomp,
                             keepX = keepX,
                             scale = scale.data)

    ## companion non-sparse model, fitted on the same data: it is not resampled, so it costs
    ## nothing, and it provides the complete loading vector needed to rank all the proteins
    dense.model = mixOmics::plsda(X = X,
                                  Y = Y,
                                  ncomp = ncomp,
                                  scale = scale.data)


    ### Fix the orientation of the components
    signs = .splsda.sign(variates = model$variates$X, Y = Y, reference.group = reference.group)
    dense.signs = .splsda.sign(variates = dense.model$variates$X, Y = Y, reference.group = reference.group)


    ### Sample scores
    variates = model$variates$X[, seq_len(ncomp), drop = FALSE]
    variates = sweep(x = variates, MARGIN = 2, STATS = signs, FUN = "*")
    variates = as.data.frame(variates)
    colnames(variates) = paste0("comp", seq_len(ncomp))
    variates$column.id = rownames(model$variates$X)

    components = dplyr::left_join(variates, prepared$metadata, by = "column.id")


    ### Loadings
    loadings = .splsda.loadings(model = model,
                                ncomp = ncomp,
                                signs = signs,
                                mat = mat,
                                Y = Y,
                                sparse = TRUE)

    full.loadings = .splsda.loadings(model = dense.model,
                                     ncomp = ncomp,
                                     signs = dense.signs,
                                     mat = mat,
                                     Y = Y,
                                     sparse = FALSE)

    full.loadings = full.loadings[, c("prot.id", "component", "loading", "abs.loading")]

    ## VIP scores of the dense model: an unsigned measure of how much each protein
    ## contributes, usable as an alternative ranking
    vip.scores = tryCatch(mixOmics::vip(dense.model), error = function(e) {return(NULL)})

    if (!is.null(vip.scores)) {
      vip.tb = data.frame(prot.id = rep(rownames(vip.scores), times = ncomp),
                          component = rep(seq_len(ncomp), each = nrow(vip.scores)),
                          vip = as.numeric(vip.scores[, seq_len(ncomp), drop = FALSE]),
                          stringsAsFactors = FALSE)

      full.loadings = dplyr::left_join(full.loadings, vip.tb, by = c("prot.id", "component"))
    } else {
      full.loadings$vip = NA_real_
    }


    ### Selected proteins, one vector per component
    selected.proteins =
      lapply(X = seq_len(ncomp),
             FUN = function(k) {loadings$prot.id[loadings$component == k & loadings$selected]})

    names(selected.proteins) = paste0("comp", seq_len(ncomp))


    ### Explained variance
    importance = .splsda.expl.var(model = model, ncomp = ncomp)



    ######################################################################################
    ### Performance
    ######################################################################################

    performance = NULL

    if (validate == TRUE) {

      message("Estimating the classification performances by cross-validation...")

      set.seed(seed)

      perf.result =
        tryCatch(.splsda.parallel(FUN = mixOmics::perf,
                                  args = list(object = model,
                                              validation = cv$validation,
                                              folds = cv$folds,
                                              nrepeat = cv$nrepeat,
                                              auc = TRUE,
                                              progressBar = FALSE),
                                  n.cores = n.cores),
                 error = function(e) {return(NULL)})

      if (is.null(perf.result)) {
        warning("The cross-validation of the performances failed: the 'performance' slot is left empty.")
      } else {
        performance = list(error.rate = .splsda.error.table(perf.result),
                           error.rate.class = tryCatch(perf.result$error.rate.class, error = function(e) {return(NULL)}),
                           choice.ncomp = tryCatch(perf.result$choice.ncomp, error = function(e) {return(NULL)}),
                           stability = .splsda.stability(perf.result = perf.result, ncomp = ncomp),
                           auc = .splsda.auc(perf.result),
                           perf.object = perf.result)
      }
    }



    ######################################################################################
    ### Ready-made plots
    ######################################################################################

    parameters = list(group.column = group.column,
                      groups = levels(Y),
                      reference.group = reference.group,
                      ncomp = ncomp,
                      keepX = keepX,
                      tuned = !is.null(tuning),
                      scale.data = scale.data,
                      validation = cv$validation,
                      folds = cv$folds,
                      nrepeat = cv$nrepeat,
                      dist = dist,
                      n.proteins = nrow(mat),
                      n.samples = ncol(mat),
                      seed = seed,
                      signs = signs)


    ### Build the object before the plots, so that the plotting functions can read it
    DEprot.sPLSDA.object =
      new(Class = "DEprot.sPLSDA",
          sPLSDA.metadata = prepared$metadata,
          sample.subset = paste(prepared$metadata$column.id, collapse = ", "),
          data.used = prepared$data.used,
          group.column = group.column,
          groups = Y,
          splsda = model,
          plsda = dense.model,
          components = components,
          loadings = loadings,
          full.loadings = full.loadings,
          selected.proteins = selected.proteins,
          importance = importance,
          counts.used = mat,
          tuning = tuning,
          performance = performance,
          cumulative.plot = NULL,
          scatter.plot.123 = NULL,
          tuning.plot = NULL,
          performance.plot = NULL,
          parameters = parameters)


    DEprot.sPLSDA.object@cumulative.plot =
      tryCatch(plot.sPLSDA.cumulative(DEprot.sPLSDA.object = DEprot.sPLSDA.object),
               error = function(e) {return(NULL)})

    DEprot.sPLSDA.object@scatter.plot.123 =
      tryCatch(plot.sPLSDA.scatter.123(DEprot.sPLSDA.object = DEprot.sPLSDA.object,
                                       color.column = group.column),
               error = function(e) {return(NULL)})

    if (!is.null(tuning)) {
      DEprot.sPLSDA.object@tuning.plot =
        tryCatch(plot.sPLSDA.tuning(DEprot.sPLSDA.object = DEprot.sPLSDA.object),
                 error = function(e) {return(NULL)})
    }

    if (!is.null(performance)) {
      DEprot.sPLSDA.object@performance.plot =
        tryCatch(plot.sPLSDA.performance(DEprot.sPLSDA.object = DEprot.sPLSDA.object),
                 error = function(e) {return(NULL)})
    }


    return(DEprot.sPLSDA.object)
  } # END function




# ----------------------------------------------------------------------------------------
###                                      GETTERS                                        ###
# ----------------------------------------------------------------------------------------


#' @title get.sPLSDA.results
#'
#' @description Simplifies the access to the loadings of a sPLS-DA. When the original object carries a protein annotation table (\code{protein.info} slot, see \link{load.counts2} and \link{add.protein.info}), the annotation columns can be appended.
#'
#' @param DEprot.sPLSDA.object An object of class \code{DEprot.sPLSDA}.
#' @param component Numeric value (or vector) indicating the component(s) to return. Default: \code{NULL} (all the components).
#' @param selected.only Logical value indicating whether only the proteins retained by the sparse model should be returned. When \code{FALSE}, the complete loadings of the companion non-sparse model are returned instead. Default: \code{TRUE}.
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses} from which the protein annotation is taken. Default: \code{NULL} (no annotation appended).
#' @param protein.info.columns String vector indicating which columns of the \code{protein.info} slot should be appended. Two keywords are available: \code{"none"} and \code{"all"}. Default: \code{"none"}.
#' @param protein.info.prefix String to prepend to the names of the appended annotation columns. Default: \code{NULL}.
#'
#' @return A data.frame.
#'
#' @seealso \link{perform.sPLSDA}, \link{get.sPLSDA.proteins}, \link{get.sPLSDA.ranking}
#'
#' @import dplyr
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' splsda <- perform.sPLSDA(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                          group.column = "condition",
#'                          keepX = 5,
#'                          validate = FALSE)
#'
#' # Proteins selected on the first component
#' get.sPLSDA.results(DEprot.sPLSDA.object = splsda, component = 1)
#'
#' @export get.sPLSDA.results

get.sPLSDA.results =
  function(DEprot.sPLSDA.object,
           component = NULL,
           selected.only = TRUE,
           DEprot.object = NULL,
           protein.info.columns = "none",
           protein.info.prefix = NULL) {

    ### check object
    if (!("DEprot.sPLSDA" %in% class(DEprot.sPLSDA.object))) {
      stop("The input must be an object of class 'DEprot.sPLSDA'.")
    }

    if (selected.only == TRUE) {
      data = DEprot.sPLSDA.object@loadings
      data = data[data$selected %in% TRUE, , drop = FALSE]
    } else {
      data = DEprot.sPLSDA.object@full.loadings
    }


    ### Component subset
    if (!is.null(component)) {
      available = sort(unique(DEprot.sPLSDA.object@loadings$component))
      missing.components = setdiff(component, available)

      if (length(missing.components) > 0) {
        stop("The component(s) ", paste(missing.components, collapse = ", "), " are not available.\n",
             "       Available components: ", paste(available, collapse = ", "))
      }

      data = data[data$component %in% component, , drop = FALSE]
    }

    data = data[order(data$component, -data$abs.loading), , drop = FALSE]
    rownames(data) = NULL


    ### Append the protein annotation, when requested
    data = .append.protein.info(data = data,
                                DEprot.object = DEprot.object,
                                protein.info.columns = protein.info.columns,
                                protein.info.prefix = protein.info.prefix)

    return(data)
  } # END function




# ----------------------------------------------------------------------------------------


#' @title get.sPLSDA.proteins
#'
#' @description Extracts the proteins retained by the sparse model. The output is a plain character vector, directly usable as \code{protein.subset} in \link{heatmap.counts} or as gene list of an over-representation analysis.
#'
#' @param DEprot.sPLSDA.object An object of class \code{DEprot.sPLSDA}.
#' @param component Numeric value (or vector) indicating the component(s) to use. Default: \code{NULL} (all the components).
#' @param mode String indicating how the components are combined. One among: 'union', 'intersect', 'list' (one vector per component). Default: \code{"union"}.
#' @param direction String indicating which side of the component to keep. One among: 'both', 'positive' (proteins higher in the reference group), 'negative'. Default: \code{"both"}.
#' @param top.n Numeric value indicating how many best-ranked proteins of each component should be returned. Default: \code{NULL} (all the selected proteins).
#'
#' @return A character vector, or a named list when \code{mode = "list"}.
#'
#' @seealso \link{perform.sPLSDA}, \link{get.sPLSDA.results}, \link{sPLSDA.enrichment}
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' splsda <- perform.sPLSDA(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                          group.column = "condition",
#'                          keepX = 5,
#'                          validate = FALSE)
#'
#' # All the selected proteins
#' get.sPLSDA.proteins(DEprot.sPLSDA.object = splsda)
#'
#' # Top 3 proteins of the first component, higher in the reference group
#' get.sPLSDA.proteins(DEprot.sPLSDA.object = splsda,
#'                     component = 1,
#'                     direction = "positive",
#'                     top.n = 3)
#'
#' @export get.sPLSDA.proteins

get.sPLSDA.proteins =
  function(DEprot.sPLSDA.object,
           component = NULL,
           mode = "union",
           direction = "both",
           top.n = NULL) {

    ### check object
    if (!("DEprot.sPLSDA" %in% class(DEprot.sPLSDA.object))) {
      stop("The input must be an object of class 'DEprot.sPLSDA'.")
    }

    if (!(tolower(mode) %in% c("union", "intersect", "list"))) {
      stop("The 'mode' must be one among: 'union', 'intersect', 'list'.")
    }

    if (!(tolower(direction) %in% c("both", "positive", "negative", "pos", "neg"))) {
      stop("The 'direction' must be one among: 'both', 'positive', 'negative'.")
    }


    loadings = DEprot.sPLSDA.object@loadings
    loadings = loadings[loadings$selected %in% TRUE, , drop = FALSE]

    available = sort(unique(DEprot.sPLSDA.object@loadings$component))

    if (is.null(component)) {
      component = available
    } else {
      missing.components = setdiff(component, available)

      if (length(missing.components) > 0) {
        stop("The component(s) ", paste(missing.components, collapse = ", "), " are not available.\n",
             "       Available components: ", paste(available, collapse = ", "))
      }
    }


    protein.list =
      lapply(X = component,
             FUN = function(k) {
               tb = loadings[loadings$component == k, , drop = FALSE]

               if (tolower(direction) %in% c("positive", "pos")) {
                 tb = tb[tb$loading > 0, , drop = FALSE]
               } else if (tolower(direction) %in% c("negative", "neg")) {
                 tb = tb[tb$loading < 0, , drop = FALSE]
               }

               tb = tb[order(-tb$abs.loading), , drop = FALSE]

               if (!is.null(top.n)) {tb = utils::head(tb, round(top.n, 0))}

               return(tb$prot.id)
             })

    names(protein.list) = paste0("comp", component)


    if (tolower(mode) == "list") {
      return(protein.list)
    } else if (tolower(mode) == "intersect") {
      return(Reduce(intersect, protein.list))
    } else {
      return(unique(unlist(protein.list, use.names = FALSE)))
    }
  } # END function




# ----------------------------------------------------------------------------------------


#' @title get.sPLSDA.ranking
#'
#' @description Builds the protein ranking used as input of a GSEA. By default the ranking comes from the companion non-sparse model stored in the \code{plsda} slot, and therefore covers every protein of the matrix.
#'
#' @param DEprot.sPLSDA.object An object of class \code{DEprot.sPLSDA}.
#' @param component Numeric value indicating the component to use. Default: \code{1}.
#' @param metric String indicating the score used for the ranking. One among: 'loading' (signed coefficient of the non-sparse model), 'vip' (unsigned VIP score of the non-sparse model), 'sparse' (signed coefficient of the sparse model). Default: \code{"loading"}.
#'
#' @details
#' The sparse loadings are a poor ranking for a GSEA and \code{metric = "sparse"} exists only for the cases in which one really wants them. Three reasons: every protein that was not selected has a coefficient of exactly zero, which leaves a handful of ranked entries and thousands of ties whose order the running-sum statistic depends on; the sparsity penalty keeps one member of a pair of correlated proteins and drops the other, which is exactly what happens inside a pathway; and a set-level aggregation over a vector that is mostly zeros is not testing what it looks like it is testing.
#'
#' A positive score means a protein higher in the \code{reference.group} chosen at the \link{perform.sPLSDA} call. With more than two classes, a component separates a mixture of classes rather than a single contrast, and the sign should be read as "towards the reference group" rather than as a fold change.
#'
#' @return A named numeric vector, sorted in decreasing order.
#'
#' @seealso \link{perform.sPLSDA}, \link{sPLSDA.enrichment}
#'
#' @import dplyr
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' splsda <- perform.sPLSDA(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                          group.column = "condition",
#'                          keepX = 5,
#'                          validate = FALSE)
#'
#' ranking <- get.sPLSDA.ranking(DEprot.sPLSDA.object = splsda, component = 1)
#'
#' @export get.sPLSDA.ranking

get.sPLSDA.ranking =
  function(DEprot.sPLSDA.object,
           component = 1,
           metric = "loading") {

    ### check object
    if (!("DEprot.sPLSDA" %in% class(DEprot.sPLSDA.object))) {
      stop("The input must be an object of class 'DEprot.sPLSDA'.")
    }

    if (!(tolower(metric) %in% c("loading", "loadings", "vip", "sparse"))) {
      stop("The 'metric' must be one among: 'loading', 'vip', 'sparse'.")
    }

    available = sort(unique(DEprot.sPLSDA.object@loadings$component))

    if (length(component) != 1) {
      stop("The 'component' must be a single numeric value.")
    }

    if (!(component %in% available)) {
      stop("The component ", component, " is not available.\n",
           "       Available components: ", paste(available, collapse = ", "))
    }


    ### Extract the scores
    if (tolower(metric) == "sparse") {
      warning("The sparse loadings are zero for every protein that was not selected: the ranking is mostly made of ties ",
              "and the enrichment scores depend on the order in which those ties are sorted. ",
              "Use 'metric = \"loading\"' unless the ties are what you are after.")

      tb = DEprot.sPLSDA.object@loadings
      tb = tb[tb$component == component, c("prot.id", "loading")]
      colnames(tb)[2] = "score"

    } else if (tolower(metric) == "vip") {
      tb = DEprot.sPLSDA.object@full.loadings
      tb = tb[tb$component == component, c("prot.id", "vip")]
      colnames(tb)[2] = "score"

      if (TRUE %in% is.na(tb$score)) {
        stop("The VIP scores are not available in this object: use 'metric = \"loading\"'.")
      }

    } else {
      tb = DEprot.sPLSDA.object@full.loadings
      tb = tb[tb$component == component, c("prot.id", "loading")]
      colnames(tb)[2] = "score"
    }


    ### Collapse the duplicated protein IDs, as done for the differential rankings
    tb = dplyr::summarise(dplyr::group_by(tb, prot.id),
                          score = mean(score, na.rm = TRUE),
                          .groups = "drop")

    tb = tb[tb$prot.id != "", , drop = FALSE]
    tb = tb[order(-tb$score), , drop = FALSE]

    ranking = tb$score
    names(ranking) = tb$prot.id

    return(ranking)
  } # END function




# ----------------------------------------------------------------------------------------
###                                    ENRICHMENT                                       ###
# ----------------------------------------------------------------------------------------


#' @title sPLSDA.enrichment
#'
#' @description Performs a Gene Set Enrichment Analysis (GSEA) or an OverRepresentation Analysis (ORA) on the results of a sPLS-DA. The GSEA is run on the complete ranking of the companion non-sparse model, the ORA on the proteins retained by the sparse one, using as background all the proteins that entered the model.
#'
#' @param DEprot.sPLSDA.object An object of class \code{DEprot.sPLSDA}.
#' @param TERM2GENE Data.frame containing two columns 'gs_name' (IDs of the gene sets) and 'gene_symbol' (indicating the gene IDs). No default.
#' @param enrichment.type String indicating the type of analyses to perform. One among: GSEA, ORA. No default.
#' @param component Numeric value indicating the component to analyze. Default: \code{1}.
#' @param gsea.rank.method String indicating the score used to rank the proteins for the GSEA. One among: 'loading', 'vip', 'sparse' (see \link{get.sPLSDA.ranking}). Default: \code{"loading"}.
#' @param direction String indicating which side of the component is tested by the ORA. One among: 'both', 'positive', 'negative'. Ignored for the GSEA. Default: \code{"both"}.
#' @param top.n Numeric value indicating how many best-ranked selected proteins should be used for the ORA. Default: \code{NULL} (all the selected proteins).
#' @param universe Character vector indicating the background gene list of the ORA. Default: \code{NULL}, meaning all the proteins that entered the model.
#' @param gsub.pattern.prot.id String indicating a pattern to be passed to gsub and to remove from the prot.id. Default: \code{NULL} (no changes in the IDs).
#' @param pvalueCutoff Numeric value indicating the adjusted pvalue cutoff on enrichment tests to report. Default: \code{0.05}.
#' @param qvalueCutoff Numeric value indicating the qvalue cutoff on enrichment tests to report as significant (only for ORA). Default: \code{0.05}.
#' @param pAdjustMethod String indicating the method to use for the p-value adjustment. One among "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: \code{"BH"}.
#' @param dotplot.n Numeric value indicating the maximum number of categories to plot in the dotplot. Default: \code{10}.
#'
#' @details
#' The object returned is a \code{DEprot.enrichResult}, the same class produced by \link{geneset.enrichment}, so the results of a sPLS-DA and those of a differential analysis can be handled by the same downstream functions (\link{NES.plot}, \link{plot.GSEA}, \link{simplify.enrichment}, \link{combine.enrichments}, \link{divergent.enrichment}).
#'
#' The background of the ORA is the list of the proteins that entered the model, which is the only sensible universe here: a protein that was never quantified had no chance of being selected.
#'
#' @return An object of class \code{DEprot.enrichResult}.
#'
#' @seealso \link{perform.sPLSDA}, \link{get.sPLSDA.ranking}, \link{geneset.enrichment}
#'
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#' @import viridis
#' @import clusterProfiler
#' @importFrom methods new
#' @importFrom utils head
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' splsda <- perform.sPLSDA(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                          group.column = "condition",
#'                          keepX = 5,
#'                          validate = FALSE)
#'
#' # Over-representation of the proteins selected by the sparse model
#' ora <- sPLSDA.enrichment(DEprot.sPLSDA.object = splsda,
#'                          TERM2GENE = DEprot::test.toolbox$geneset,
#'                          enrichment.type = "ORA",
#'                          component = 1,
#'                          pvalueCutoff = 1,
#'                          qvalueCutoff = 1)
#'
#' # GSEA on the complete ranking of the non-sparse model
#' gsea <- sPLSDA.enrichment(DEprot.sPLSDA.object = splsda,
#'                           TERM2GENE = DEprot::test.toolbox$geneset,
#'                           enrichment.type = "GSEA",
#'                           component = 1,
#'                           pvalueCutoff = 1)
#'
#' @export sPLSDA.enrichment

sPLSDA.enrichment =
  function(DEprot.sPLSDA.object,
           TERM2GENE,
           enrichment.type,
           component = 1,
           gsea.rank.method = "loading",
           direction = "both",
           top.n = NULL,
           universe = NULL,
           gsub.pattern.prot.id = NULL,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           dotplot.n = 10) {

    ### check object
    if (!("DEprot.sPLSDA" %in% class(DEprot.sPLSDA.object))) {
      stop("The input must be an object of class 'DEprot.sPLSDA'.")
    }

    if (!(toupper(enrichment.type) %in% c("ORA", "GSEA", "OVER"))) {
      stop("The 'enrichment.type' must be one among: 'ORA', 'GSEA'.")
    }

    ### Check TERM2GENE
    if (!("data.frame" %in% class(TERM2GENE))) {
      stop("The 'TERM2GENE' must be a 2-column data.frame. One of this columns must be named 'gs_name' (gs: gene set)")
    } else if (!("gs_name" %in% colnames(TERM2GENE)) | ncol(TERM2GENE) > 2) {
      stop("The 'TERM2GENE' must be a 2-column data.frame. One of this columns must be named 'gs_name' (gs: gene set)")
    }


    parameters = DEprot.sPLSDA.object@parameters


    ### Labels of the two sides of the component, used by the downstream plotting functions
    ## with more than two classes a component separates the reference group from a mixture,
    ## which is what the label says
    if (length(parameters$groups) == 2) {
      var.2 = setdiff(parameters$groups, parameters$reference.group)
    } else {
      var.2 = "other groups"
    }

    contrast.info = list(metadata.column = parameters$group.column,
                         var.1 = parameters$reference.group,
                         var.2 = var.2,
                         group.1 = DEprot.sPLSDA.object@sPLSDA.metadata$column.id[DEprot.sPLSDA.object@groups == parameters$reference.group],
                         group.2 = DEprot.sPLSDA.object@sPLSDA.metadata$column.id[DEprot.sPLSDA.object@groups != parameters$reference.group],
                         component = component)

    plot.title = paste0(parameters$group.column, " (component ", component, "): **",
                        contrast.info$var.1, "** *vs* **", contrast.info$var.2, "**")



    ######################################################################################
    ### GSEA
    ######################################################################################

    if (toupper(enrichment.type) == "GSEA") {

      rank.method = tolower(gsea.rank.method)

      gene.list = get.sPLSDA.ranking(DEprot.sPLSDA.object = DEprot.sPLSDA.object,
                                     component = component,
                                     metric = rank.method)

      if (!is.null(gsub.pattern.prot.id)) {
        names(gene.list) = gsub(gsub.pattern.prot.id, "", names(gene.list))
        gene.list = gene.list[names(gene.list) != ""]
      }

      enrichment.discovery = tryCatch(clusterProfiler::GSEA(geneList = gene.list,
                                                            pvalueCutoff = pvalueCutoff,
                                                            pAdjustMethod = pAdjustMethod,
                                                            TERM2GENE = TERM2GENE),
                                      error = function(e) {return(NULL)})

      dotplot_fold.enrichment = NULL
      used.universe = NULL


      ######################################################################################
      ### ORA
      ######################################################################################

    } else {

      rank.method = NA

      genes = get.sPLSDA.proteins(DEprot.sPLSDA.object = DEprot.sPLSDA.object,
                                  component = component,
                                  mode = "union",
                                  direction = direction,
                                  top.n = top.n)

      if (length(genes) == 0) {
        stop("No protein was selected on the component ", component,
             " with 'direction = \"", direction, "\"': there is nothing to test.")
      }

      ## the correct background is the list of proteins that entered the model: using the
      ## whole proteome would make almost every geneset look enriched
      if (is.null(universe)) {universe = rownames(DEprot.sPLSDA.object@counts.used)}

      if (!is.null(gsub.pattern.prot.id)) {
        genes = gsub(gsub.pattern.prot.id, "", genes)
        universe = gsub(gsub.pattern.prot.id, "", universe)
      }

      enrichment.discovery = tryCatch(clusterProfiler::enricher(gene = genes,
                                                                universe = universe,
                                                                TERM2GENE = TERM2GENE,
                                                                pvalueCutoff = pvalueCutoff,
                                                                qvalueCutoff = qvalueCutoff,
                                                                pAdjustMethod = pAdjustMethod),
                                      error = function(e) {return(NULL)})

      used.universe = universe

      dotplot_fold.enrichment =
        tryCatch(clusterProfiler::dotplot(enrichment.discovery,
                                          x = "FoldEnrichment",
                                          showCategory = dotplot.n) +
                   ggtitle(plot.title) +
                   viridis::scale_fill_viridis(option = "mako", direction = -1, begin = 0.3) +
                   theme(plot.title = ggtext::element_markdown(hjust = 0.5),
                         axis.ticks.y = element_blank()),
                 error = function(e) {return(NULL)})
    }


    if (is.null(enrichment.discovery)) {
      warning("The enrichment failed: no result could be computed on the component ", component, ".")
    }



    ######################################################################################
    ### Plots shared by the two modes
    ######################################################################################

    pathway.network.clusters = tryCatch(aPEAR::findPathClusters(enrichment.discovery@result),
                                        error = function(e) {return(NULL)})

    if (toupper(enrichment.type) == "GSEA") {
      pathway.network.plot = tryCatch(aPEAR::enrichmentNetwork(enrichment.discovery@result,
                                                               drawEllipses = TRUE,
                                                               colorBy = "NES"),
                                      error = function(e) {return(NULL)})
    } else {
      pathway.network.plot = tryCatch(aPEAR::enrichmentNetwork(enrichment.discovery@result,
                                                               colorBy = "pvalue",
                                                               colorType = "pval",
                                                               pCutoff = -5,
                                                               nodeSize = "Count"),
                                      error = function(e) {return(NULL)})
    }

    protein.network = tryCatch(clusterProfiler::cnetplot(enrichment.discovery),
                               error = function(e) {return(NULL)})

    dotplot_gene.ratio =
      tryCatch(clusterProfiler::dotplot(enrichment.discovery,
                                        x = "GeneRatio",
                                        showCategory = dotplot.n) +
                 ggtitle(plot.title) +
                 viridis::scale_fill_viridis(option = "rocket", direction = -1, begin = 0.3) +
                 theme(plot.title = ggtext::element_markdown(hjust = 0.5),
                       axis.ticks.y = element_blank()),
               error = function(e) {return(NULL)})



    ######################################################################################
    ### Build the object
    ######################################################################################

    DEprot.enrichResult.object =
      new(Class = "DEprot.enrichResult",
          enrichment.discovery = enrichment.discovery,
          protein.network = protein.network,
          pathway.network = list(clusters = pathway.network.clusters,
                                 plot = pathway.network.plot),
          NES.plot = NULL,
          dotplot_gene.ratio = dotplot_gene.ratio,
          dotplot_fold.enrichment = dotplot_fold.enrichment,
          parameters = list(enrichment.type = toupper(enrichment.type),
                            contrast = contrast.info,
                            source = "sPLS-DA",
                            component = component,
                            direction = direction,
                            top.n = top.n,
                            universe = used.universe,
                            diff.status.category = NULL,
                            gsub.pattern.prot.id = gsub.pattern.prot.id,
                            gsea.rank.method = rank.method,
                            pvalueCutoff = pvalueCutoff,
                            qvalueCutoff = qvalueCutoff,
                            pAdjustMethod = pAdjustMethod,
                            dotplot.n = dotplot.n),
          affinity.propagation = NULL)


    ## the NES barplot is built by the exported function, which reads the labels of the two
    ## sides from the object itself
    if (toupper(enrichment.type) == "GSEA") {
      DEprot.enrichResult.object@NES.plot =
        tryCatch(NES.plot(enrichResult = DEprot.enrichResult.object, title = plot.title),
                 error = function(e) {return(NULL)})
    }

    return(DEprot.enrichResult.object)
  } # END function
