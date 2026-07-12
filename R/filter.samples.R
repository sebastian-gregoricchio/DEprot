#' @title filter.samples
#'
#' @description Removes (or keeps) a subset of samples from an object of class
#' \code{DEprot} or \code{DEprot.analyses} and \strong{fully re-derives every
#' processing step} that had been applied to the original object. Because
#' normalization (MBQN / HarmonizR), bottom-distribution randomization and,
#' most importantly, imputation (e.g. \code{missForest}) depend on the whole
#' set of samples/sub-groups, dropping samples and simply sub-setting the
#' pre-computed matrices would produce values that are inconsistent with what
#' those algorithms would return on the reduced data. This function therefore
#' goes back to the lowest count level that can be faithfully re-derived
#' (\code{raw} > \code{normalized} > \code{randomized} > \code{imputed}),
#' rebuilds a fresh object with the retained samples and re-runs, in order, the
#' same steps as in the original object, recovering the parameters directly from
#' the original object's slots and re-using the existing DEprot functions
#' (\link{normalize.counts}, \link{harmonize.batches},
#' \link{randomize.missing.values}, \link{impute.counts} and, for
#' \code{DEprot.analyses}, \link{diff.analyses} / \link{diff.analyses.limma} /
#' \link{diff.analyses.prolfqua}).
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param samples Character vector with the identifiers of the samples to remove
#'   or to keep. Identifiers must correspond to values of the \code{column.id}
#'   column of the object metadata (i.e. the column names of the counts tables).
#' @param mode String, one among \code{"remove"} and \code{"keep"}. If
#'   \code{"remove"} (default) the samples listed in \code{samples} are dropped;
#'   if \code{"keep"} only the samples listed in \code{samples} are retained.
#' @param batch.column String indicating the metadata column that stores the
#'   batch identifiers, used only when the original object was batch-corrected
#'   with \link{harmonize.batches}. This is normally recovered from the object,
#'   so this argument acts as an override; if it is not stored (legacy objects)
#'   and this argument is \code{NULL} (default) the function looks for a column
#'   literally called \code{"batch"}. Default: \code{NULL}.
#' @param diff.method String defining which differential-analyses engine should
#'   be used to recompute the contrasts of a \code{DEprot.analyses} object. One
#'   among \code{"auto"}, \code{"t.test"}, \code{"limma"} and \code{"prolfqua"}.
#'   With \code{"auto"} (default) the engine is read directly from the object
#'   (the \code{stat.test} entry of \code{differential.analyses.params} records
#'   whether \code{limma}, \code{prolfqua} or the t/Wilcoxon engine was used);
#'   for legacy objects that do not store it, \code{prolfqua} is still detected
#'   and the function otherwise falls back to \code{diff.analyses} (t/Wilcoxon).
#'   Set this argument explicitly only to override the stored engine. Default:
#'   \code{"auto"}.
#' @param stat.test String passed to \link{diff.analyses} (t/Wilcoxon engine
#'   only), one among \code{"t.test"} and \code{"wilcoxon"}. Normally recovered
#'   from the object, so this argument acts as an override. If \code{NULL}
#'   (default) the stored test is used, falling back to \code{"t.test"}. Default:
#'   \code{NULL}.
#' @param replicate.column String indicating the metadata column with the
#'   replicate identifiers. It is recovered from the object for every engine, so
#'   this argument only acts as an override (e.g. if the stored column was
#'   renamed). Default: \code{NULL}.
#' @param min.samples.per.group Numeric value indicating the minimum number of
#'   retained samples that each group of a contrast must still contain for the
#'   contrast to be recomputed (\code{DEprot.analyses} only). Contrasts in which
#'   at least one group falls below this threshold are skipped with a warning.
#'   Default: \code{2}.
#' @param verbose Logical value indicating whether progress messages should be
#'   printed. Imputation can be slow, hence a message is printed before each
#'   step. Default: \code{TRUE}.
#'
#' @details
#' \strong{Parameter recovery.} The parameters of each step are collected from
#' the corresponding slot of the original object:
#' \itemize{
#'   \item \emph{Normalization} (\code{normalization.method} slot): a
#'         \code{data.frame} whose \code{"package"} row is either \code{"MBQN"}
#'         or \code{"HarmonizR"}. For MBQN, \code{balancing.function} and
#'         \code{NRI.RI.ratio.threshold} are recovered. For HarmonizR, the
#'         \code{batch.column}, \code{algorithm}, \code{ComBat.mode},
#'         \code{block} and \code{cores} are recovered.
#'   \item \emph{Randomization} (\code{randomization.method} slot): all the
#'         parameters (\code{group.column}, \code{percentage.missing},
#'         \code{tail.percentage}, \code{which.data} and \code{seed}) are
#'         recovered.
#'   \item \emph{Imputation} (\code{imputation.method} slot): \code{method},
#'         \code{seed}, \code{data.used} (the counts that were imputed) and the
#'         method-specific parameters are recovered -- including, for
#'         \code{missForest}, the \code{parallel.mode} and
#'         \code{variable.wise.OOBerror}, and for \code{SVD}/\code{BPCA}/\code{PPCA}
#'         the number of principal components tested.
#'   \item \emph{Differential analyses} (\code{differential.analyses.params} and
#'         \code{contrasts} slots): the contrasts, fold-change thresholds,
#'         p-value/FDR threshold, adjustment method and counts used are recovered.
#'         The engine is identified from the stored \code{stat.test}
#'         (\code{"limma"}, \code{"prolfqua"} or the t/Wilcoxon test name); the
#'         \code{replicate.column} is recovered for every engine, together with
#'         the t/Wilcoxon \code{paired.test}, the limma \code{rep.model} /
#'         \code{fitting.method}, and the prolfqua \code{strategy} /
#'         \code{moderate.variance}.
#' }
#' The seeds saved during the original randomization/imputation are re-used so
#' the re-computation is as reproducible as possible on the reduced data.
#'
#' \strong{HarmonizR fallback.} When the original object was harmonized, the
#' re-harmonization is attempted with the same \code{algorithm}/\code{ComBat.mode}.
#' If, after removing samples, a batch becomes too small for ComBat (or ComBat
#' fails for any reason), the function automatically falls back to the
#' \code{"limma"} algorithm.
#'
#' @return An object of the same class as the input (\code{DEprot} or
#'   \code{DEprot.analyses}) containing only the retained samples, with all the
#'   count tables, boxplots and (when applicable) differential results
#'   recomputed from scratch.
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' \dontrun{
#' # Remove two samples from an imputed DEprot object
#' dpo.filtered <- filter.samples(DEprot.object = dpo.imputed,
#'                                samples = c("sample_3", "sample_7"),
#'                                mode = "remove")
#'
#' # Keep only a subset of samples from a DEprot.analyses object
#' dpa.filtered <- filter.samples(DEprot.object = dpo.analyses,
#'                                samples = c("ctrl_1", "ctrl_2", "treat_1", "treat_2"),
#'                                mode = "keep",
#'                                replicate.column = "replicate")
#'
#' # Object that had been batch-corrected: provide the batch column
#' dpo.filtered <- filter.samples(DEprot.object = dpo.harmonized,
#'                                samples = "outlier_batchB",
#'                                mode = "remove",
#'                                batch.column = "batch")
#' }
#'
#' @import dplyr
#' @importFrom methods slot slotNames is
#'
#' @export filter.samples

filter.samples =
  function(DEprot.object,
           samples,
           mode = "remove",
           batch.column = NULL,
           diff.method = "auto",
           stat.test = NULL,
           replicate.column = NULL,
           min.samples.per.group = 2,
           verbose = TRUE) {

    ############################################################################
    ### small local helpers
    ############################################################################
    `%||%` = function(a, b) {if (is.null(a) || (length(a) == 1 && is.na(a))) b else a}

    say = function(...) {if (isTRUE(verbose)) {message(...)}}

    # run an expression silencing only *messages* (warnings are always kept)
    # when verbose = FALSE
    quietly = function(expr) {if (isTRUE(verbose)) {force(expr)} else {suppressMessages(expr)}}




    ############################################################################
    ### Check input object
    ############################################################################
    obj.class = class(DEprot.object)

    if (!any(c("DEprot", "DEprot.analyses") %in% obj.class)) {
      stop("The input must be an object of class 'DEprot' or 'DEprot.analyses'.", call. = FALSE)
    }
    is.analyses = "DEprot.analyses" %in% obj.class

    mode = tolower(mode[1])
    if (!(mode %in% c("remove", "keep"))) {
      stop("'mode' must be one among 'remove' and 'keep'.", call. = FALSE)
    }

    if (missing(samples) || is.null(samples) || length(samples) == 0) {
      stop("Provide 'samples': a character vector of sample IDs (values of the 'column.id' metadata column) to remove or keep.", call. = FALSE)
    }
    samples = as.character(samples)

    diff.method = tolower(diff.method[1])
    if (!(diff.method %in% c("auto", "t.test", "ttest", "wilcoxon", "limma", "prolfqua"))) {
      stop("'diff.method' must be one among 'auto', 't.test', 'limma' and 'prolfqua'.", call. = FALSE)
    }





    ############################################################################
    ### Check that the requested samples exist and define the retained ones
    ############################################################################
    meta = DEprot.object@metadata
    if (is.null(meta) || !("column.id" %in% colnames(meta))) {
      stop("The metadata table of the object must contain a 'column.id' column.", call. = FALSE)
    }
    all.samples = as.character(meta$column.id)

    missing.samples = setdiff(samples, all.samples)
    if (length(missing.samples) > 0) {
      stop(paste0("The following sample(s) are not present in the 'column.id' column of the metadata table:\n  ",
                  paste(missing.samples, collapse = ", "),
                  "\n\nAvailable samples:\n  ",
                  paste(all.samples, collapse = ", ")),
           call. = FALSE)
    }

    if (mode == "remove") {
      keep.samples = setdiff(all.samples, samples)
    } else { # keep
      keep.samples = intersect(all.samples, samples)
    }
    # preserve the original column order
    keep.samples = all.samples[all.samples %in% keep.samples]

    if (length(keep.samples) == 0) {
      stop("The requested operation would remove all the samples. Nothing to return.", call. = FALSE)
    }
    if (length(keep.samples) < 2) {
      warning(paste0("Only ", length(keep.samples),
                     " sample would be retained: normalization/imputation and differential analyses require at least 2 samples and will most likely fail."),
              call. = FALSE, immediate. = TRUE)
    }

    n.removed = length(all.samples) - length(keep.samples)
    say("[filter.samples] Retaining ", length(keep.samples), "/", length(all.samples),
        " samples (", n.removed, " removed).")




    ############################################################################
    ### Collect the state / parameters of the original object
    ############################################################################
    was.normalized = isTRUE(DEprot.object@normalized)
    was.randomized = isTRUE(DEprot.object@randomized)
    was.imputed = isTRUE(DEprot.object@imputed)

    norm.method = DEprot.object@normalization.method
    random.params = DEprot.object@randomization.method
    imp.params = DEprot.object@imputation.method

    log.base = DEprot.object@log.base %||% 2


    # availability of each count level (non-null and non-empty)
    non_empty = function(x) {!is.null(x) && length(dim(x)) == 2 && nrow(x) > 0 && ncol(x) > 0}

    has.raw = non_empty(DEprot.object@raw.counts)
    has.norm = non_empty(DEprot.object@norm.counts)
    has.random = non_empty(DEprot.object@random.counts)
    has.imputed = non_empty(DEprot.object@imputed.counts)


    # helper to read a value from a normalization.method data.frame (param/value)
    norm.param = function(name) {
      if (!is.data.frame(norm.method)) {return(NULL)}
      idx = which(as.character(norm.method$param) == name)
      if (length(idx) == 0) {return(NULL)}
      as.character(norm.method$value[idx[1]])
    }

    # Is the normalization step something we can re-derive from raw counts?
    # Both MBQN and HarmonizR store a data.frame in 'normalization.method'; they
    # are told apart by the "package" row ("MBQN" vs "HarmonizR").
    norm.package = norm.param("package")
    norm.is.mbqn = identical(norm.package, "MBQN")
    norm.is.harmonizr = identical(norm.package, "HarmonizR")
    norm.rederivable = norm.is.mbqn || norm.is.harmonizr




    ############################################################################
    ### Decide the starting level and rebuild a fresh (raw) DEprot object
    ############################################################################
    # We restart from the lowest level that we can faithfully re-derive.
    if (has.raw && (!was.normalized || norm.rederivable)) {
      start.level = "raw"
    } else if (has.norm) {
      start.level = "normalized"
    } else if (has.random) {
      start.level = "randomized"
    } else if (has.imputed) {
      start.level = "imputed"
    } else {
      stop("The object does not contain any usable count table (raw/normalized/randomized/imputed).", call. = FALSE)
    }

    if (start.level != "raw") {
      if (was.normalized && !norm.rederivable && start.level == "normalized") {
        warning("No 'raw.counts' available (or normalization method not recoverable): the normalization step cannot be re-derived, ",
                "the existing normalized counts of the retained samples are used as starting point.",
                call. = FALSE, immediate. = TRUE)
      } else {
        warning(paste0("The lowest count table that can be used is '", start.level,
                       "': steps below this level cannot be re-derived and the corresponding counts are simply sub-set. ",
                       "For a fully consistent re-computation, provide an object that still contains the 'raw' counts."),
                call. = FALSE, immediate. = TRUE)
      }
    }

    level.slot = c(raw = "raw.counts", normalized = "norm.counts",
                   randomized = "random.counts", imputed = "imputed.counts")[[start.level]]
    level.type = c(raw = "raw", normalized = "normalized",
                   randomized = "randomized", imputed = "imputed")[[start.level]]

    start.counts = methods::slot(DEprot.object, level.slot)
    start.counts = start.counts[, keep.samples, drop = FALSE]

    meta.sub = meta[match(keep.samples, as.character(meta$column.id)), , drop = FALSE]
    rownames(meta.sub) = NULL

    say("[filter.samples] Rebuilding the object from the '", start.level, "' counts...")
    dpo = quietly(DEprot::load.counts2(counts = start.counts,
                                       metadata = meta.sub,
                                       data.type = level.type,
                                       log.base = log.base,
                                       column.id = "column.id"))





    ############################################################################
    ### Re-normalize / re-harmonize (only when starting from raw counts)
    ############################################################################
    if (start.level == "raw" && was.normalized) {

      if (norm.is.harmonizr) {
        ## batch harmonization (HarmonizR) -----------------------------------------------
        # recover the parameters directly from the normalization.method
        # data.frame (package / batch.column / algorithm / ComBat.mode / block /
        # cores).
        algorithm = norm.param("algorithm") %||% "ComBat"
        combat.mode = suppressWarnings(as.numeric(norm.param("ComBat.mode")))
        if (length(combat.mode) == 0 || is.na(combat.mode)) {combat.mode = 1}
        block.stored = norm.param("block")
        hblock = if (is.null(block.stored) || tolower(block.stored) %in% c("null", "na", "")) {NULL} else {block.stored}
        hcores = suppressWarnings(as.numeric(norm.param("cores")))
        if (length(hcores) == 0 || is.na(hcores)) {hcores = 1}
        stored.bcol = norm.param("batch.column")

        # batch column: the 'batch.column' argument (if given) overrides the
        # value stored in the object; otherwise use the stored one, then a column
        # literally called 'batch'.
        bcol = batch.column %||% stored.bcol
        if (is.null(bcol)) {
          if ("batch" %in% colnames(meta.sub)) {
            bcol = "batch"
            say("[filter.samples] batch column not available: using the 'batch' metadata column for harmonization.")
          } else {
            stop("The original object was batch-corrected with 'harmonize.batches' but the batch column could not be determined ",
                 "(not stored and no column called 'batch' was found in the metadata).\n",
                 "Please provide the batch column via the 'batch.column' argument.", call. = FALSE)
          }
        } else if (!(bcol %in% colnames(meta.sub))) {
          stop(paste0("The 'batch.column' ('", bcol, "') is not present in the metadata of the object."), call. = FALSE)
        }

        # proactively switch to limma if a batch became too small for ComBat
        if (tolower(algorithm) == "combat") {
          batch.sizes = table(as.character(meta.sub[[bcol]]))
          if (any(batch.sizes < 2)) {
            warning(paste0("After sample removal at least one batch contains fewer than 2 samples: ComBat cannot be applied. ",
                           "Falling back to the 'limma' algorithm for the batch correction."),
                    call. = FALSE, immediate. = TRUE)
            algorithm = "limma"
          }
        }

        say("[filter.samples] Re-harmonizing batches (HarmonizR, algorithm = '", algorithm,
            if (tolower(algorithm) == "combat") {paste0("', ComBat.mode = ", combat.mode)} else {"'"}, ")...")

        dpo = tryCatch(
          quietly(DEprot::harmonize.batches(DEprot.object = dpo,
                                            batch.column = bcol,
                                            algorithm = algorithm,
                                            ComBat.mode = combat.mode,
                                            block = hblock,
                                            cores = hcores,
                                            verbose = verbose)),
          error = function(e) {
            if (tolower(algorithm) == "limma") {stop(e)}
            warning(paste0("Re-harmonization with algorithm '", algorithm, "' failed (",
                           conditionMessage(e), "). Falling back to the 'limma' algorithm."),
                    call. = FALSE, immediate. = TRUE)
            quietly(DEprot::harmonize.batches(DEprot.object = dpo,
                                              batch.column = bcol,
                                              algorithm = "limma",
                                              block = hblock,
                                              cores = hcores,
                                              verbose = verbose))
          })

      } else if (norm.is.mbqn) {
        ## MBQN normalization ------------------------------------------------------------
        # recover the parameters from the normalization.method data.frame
        get.norm.param = function(param.name) {
          idx = which(as.character(norm.method$param) == param.name)
          if (length(idx) == 0) {return(NULL)}
          as.character(norm.method$value[idx[1]])
        }
        bal.fun = get.norm.param("function")
        if (is.null(bal.fun) || is.na(bal.fun) || bal.fun == "NA") {bal.fun = "median"}
        nri.ri  = suppressWarnings(as.numeric(get.norm.param("NRI/RI ratio threshold")))
        if (length(nri.ri) == 0 || is.na(nri.ri)) {nri.ri = 0.5}

        say("[filter.samples] Re-normalizing (MBQN, balancing.function = '", bal.fun,
            "', NRI.RI.ratio.threshold = ", nri.ri, ")...")
        dpo = quietly(DEprot::normalize.counts(DEprot.object = dpo,
                                               balancing.function = bal.fun,
                                               NRI.RI.ratio.threshold = nri.ri))
      }
    }




    ############################################################################
    ### Re-randomize the missing values
    ############################################################################
    if (was.randomized) {
      if (!is.list(random.params)) {
        stop("The object is flagged as randomized but the 'randomization.method' slot does not contain the expected list of parameters.", call. = FALSE)
      }

      group.col = random.params$group.column
      if (is.null(group.col) || !(group.col %in% colnames(meta.sub))) {
        stop(paste0("The randomization group column ('", group.col %||% "NA",
                    "') stored in the object is not present in the metadata table."), call. = FALSE)
      }

      rand.which = random.params$data.used %||% "normalized"
      # make sure the requested level is available in the rebuilt object;
      # otherwise fall back to the closest available lower level
      rand.avail = c(raw = non_empty(dpo@raw.counts),
                     normalized = non_empty(dpo@norm.counts))
      if (rand.which %in% c("raw") && !isTRUE(rand.avail[["raw"]])) {rand.which = "normalized"}
      if (rand.which %in% c("normalized", "norm") && !isTRUE(rand.avail[["normalized"]])) {rand.which = "raw"}

      say("[filter.samples] Re-randomizing missing values (group.column = '", group.col,
          "', which.data = '", rand.which, "')...")
      dpo = DEprot::randomize.missing.values(DEprot.object = dpo,
                                             group.column = group.col,
                                             percentage.missing = random.params$percentage.missing %||% 100,
                                             tail.percentage = random.params$tail.percentage %||% 3,
                                             which.data = rand.which,
                                             seed = random.params$seed %||% floor(stats::runif(1, 0, 50000)),
                                             verbose = verbose)
    }


    ############################################################################
    ### Re-impute
    ############################################################################
    if (was.imputed) {
      if (!is.list(imp.params) || is.null(imp.params$method)) {
        stop("The object is flagged as imputed but the 'imputation.method' slot does not contain the expected list of parameters.", call. = FALSE)
      }

      imp.method = imp.params$method

      # which.data used for the imputation is now stored ('data.used'); for legacy
      # objects it is inferred from the levels present in the ORIGINAL object (the
      # same order impute.counts would have used by default).
      imp.which = imp.params$data.used %||%
        (if (was.randomized) {"randomized"} else if (was.normalized) {"normalized"} else {"raw"})

      # base arguments common to every method
      imp.args = list(DEprot.object = dpo,
                      method        = imp.method,
                      which.data    = imp.which,
                      seed          = imp.params$seed,
                      verbose       = verbose)

      # method-specific arguments (recovered from the imputation list)
      m = tolower(imp.method)
      if (m == "missforest") {
        if (!is.null(imp.params$max.iterations)) {imp.args$missForest.max.iterations = imp.params$max.iterations}
        if (!is.null(imp.params$cores))          {imp.args$missForest.cores          = imp.params$cores}
        # 'parallel.mode' holds the user-requested mode; 'parallelization.mode' is
        # forced to 'none' when cores <= 1, so it is only used as a legacy fallback.
        pmode = imp.params$parallel.mode %||% imp.params$parallelization.mode
        if (!is.null(pmode) && tolower(pmode) %in% c("variables", "forests")) {
          imp.args$missForest.parallel.mode = tolower(pmode)
        }
        if (!is.null(imp.params$variable.wise.OOBerror)) {
          imp.args$missForest.variable.wise.OOBerror = as.logical(imp.params$variable.wise.OOBerror)
        }
      } else if (m == "knn") {
        if (!is.null(imp.params$n.nearest.neighbours)) {imp.args$kNN.n.nearest.neighbours = imp.params$n.nearest.neighbours}
      } else if (m == "lls") {
        if (!is.null(imp.params$cluster.size)) {imp.args$LLS.k = imp.params$cluster.size}
      } else if (m == "regimpute") {
        if (!is.null(imp.params$parameters$fillmethod))       {imp.args$RegImpute.fillmethod     = imp.params$parameters$fillmethod}
        if (!is.null(imp.params$parameters$maxiter_RegImpute)) {imp.args$RegImpute.max.iterations = imp.params$parameters$maxiter_RegImpute}
      } else if (m %in% c("svd", "bpca", "ppca")) {
        if (!is.null(imp.params$PCs.tested)) {imp.args$pcaMethods.nPCs.to.test = imp.params$PCs.tested}
      }
      # tkNN / corkNN: method + seed are sufficient (their neighbour count is
      # derived internally from the data dimensions).

      say("[filter.samples] Re-imputing missing values (method = '", imp.method,
          "', which.data = '", imp.which, "'). This step can take some time...")
      dpo = do.call(DEprot::impute.counts, imp.args)
      say("[filter.samples] Imputation completed.")
    }





    ############################################################################
    ### Re-run the differential analyses (DEprot.analyses only)
    ############################################################################
    if (is.analyses) {

      contrasts.info = DEprot.object@contrasts
      diff.params    = DEprot.object@differential.analyses.params

      if (is.null(contrasts.info) || length(contrasts.info) == 0) {
        warning("The object is of class 'DEprot.analyses' but does not contain any contrast: the differential analyses are not recomputed and a 'DEprot' object is returned.",
                call. = FALSE, immediate. = TRUE)
        return(dpo)
      }

      # rebuild the contrast.list (metadata.column, var.1, var.2) and check that
      # every group still has enough retained samples
      contrast.list = list()
      any.paired = FALSE
      dropped.contrasts = character(0)

      for (i in seq_along(contrasts.info)) {
        ci = contrasts.info[[i]]
        mcol = ci$metadata.column
        v1   = ci$var.1
        v2   = ci$var.2
        any.paired = any.paired || isTRUE(ci$paired.test)

        if (is.null(mcol) || !(mcol %in% colnames(meta.sub))) {
          dropped.contrasts = c(dropped.contrasts, paste0(mcol %||% "NA", ": ", v1, " vs ", v2))
          next
        }
        n1 = sum(as.character(meta.sub[[mcol]]) == v1)
        n2 = sum(as.character(meta.sub[[mcol]]) == v2)
        if (n1 < min.samples.per.group || n2 < min.samples.per.group) {
          dropped.contrasts = c(dropped.contrasts, paste0(mcol, ": ", v1, " (n=", n1, ") vs ", v2, " (n=", n2, ")"))
          next
        }
        contrast.list[[length(contrast.list) + 1]] = c(mcol, v1, v2)
      }

      if (length(dropped.contrasts) > 0) {
        warning(paste0("The following contrast(s) are skipped because at least one group has fewer than ",
                       min.samples.per.group, " retained sample(s):\n  ",
                       paste(dropped.contrasts, collapse = "\n  ")),
                call. = FALSE, immediate. = TRUE)
      }
      if (length(contrast.list) == 0) {
        warning("None of the original contrasts can be recomputed with the retained samples: a 'DEprot' object (without differential analyses) is returned.",
                call. = FALSE, immediate. = TRUE)
        return(dpo)
      }

      # recover the shared differential parameters
      diff.params = if (is.list(diff.params)) {diff.params} else {list()}
      linear.FC.th = diff.params$linear.FC.th %||% 2
      linear.FC.unresp.range = diff.params$linear.FC.unresp.range %||% c(1/1.1, 1.1)
      padj.th = diff.params$padj.th %||% 0.05
      padj.method = diff.params$padj.method %||% "BH"
      diff.which = diff.params$counts.used %||% (if (was.imputed) {"imputed"} else if (was.randomized) {"randomized"} else if (was.normalized) {"normalized"} else {"raw"})

      # the engine is now recorded in 'stat.test' ("limma", "prolfqua", or the
      # actual t/Wilcoxon test name); legacy objects fall back to detection.
      stored.stat = diff.params$stat.test

      # choose the engine
      engine = diff.method
      if (engine == "ttest") {engine = "t.test"}
      if (engine == "wilcoxon") {engine = "t.test"; if (is.null(stat.test)) {stat.test = "wilcoxon"}}
      if (engine == "auto") {
        if (!is.null(stored.stat)) {
          s = tolower(as.character(stored.stat))
          engine = if (s == "limma") {"limma"} else if (s == "prolfqua") {"prolfqua"} else {"t.test"}
        } else {
          # legacy objects: prolfqua is detectable, t.test vs limma is not
          is.prolfqua = identical(as.character(padj.method), "effective FDR") ||
            !is.null(diff.params$strategy) ||
            !is.null(diff.params$moderate.variance)
          if (is.prolfqua) {
            engine = "prolfqua"
          } else {
            engine = "t.test"
            warning("The differential-analyses engine is not stored in this (legacy) object: assuming the t/Wilcoxon engine ('diff.analyses'). ",
                    "If the original analysis used limma, re-run with diff.method = 'limma'.",
                    call. = FALSE, immediate. = TRUE)
          }
        }
      }

      # engine-specific parameter recovery
      #   t/Wilcoxon : 'stat.test' (the test name) and 'paired.test' are stored
      stat.test.value = stat.test %||%
        (if (!is.null(stored.stat) && !(tolower(as.character(stored.stat)) %in% c("limma", "prolfqua"))) {as.character(stored.stat)} else {NULL}) %||%
        "t.test"
      paired.test.value = if (!is.null(diff.params$paired.test)) {isTRUE(diff.params$paired.test)} else {isTRUE(any.paired)}

      # limma: 'rep.model' (include.rep.model), 'replicate.column' and 'fitting.method' are stored
      include.rep.model.value = if (!is.null(diff.params$rep.model)) {isTRUE(diff.params$rep.model)} else {isTRUE(any.paired)}
      fitting.method.value    = diff.params$fitting.method %||% "ls"

      # prolfqua: 'strategy' (strategy.id) and 'moderate.variance' are stored
      prolfqua.strategy = diff.params$strategy$strategy.id %||% "lm"
      moderate.value    = diff.params$moderate.variance %||% FALSE

      # replicate column: every engine now stores it, so the 'replicate.column'
      # argument only acts as an override and the stored value is used otherwise.
      rep.col = replicate.column %||% diff.params$replicate.column

      # is a replicate column actually required by the chosen configuration?
      needs.rep = (engine == "t.test"   && isTRUE(paired.test.value)) ||
        (engine == "limma"    && isTRUE(include.rep.model.value)) ||
        (engine == "prolfqua" && tolower(prolfqua.strategy) %in% c("lmer", "logistf"))
      if (isTRUE(needs.rep) && is.null(rep.col)) {
        stop("The original differential analyses require the replicate column (paired test or mixed-effects model) but it is not available. ",
             "Please provide it via the 'replicate.column' argument.", call. = FALSE)
      }
      if (!is.null(rep.col) && !(rep.col %in% colnames(meta.sub))) {
        stop(paste0("The 'replicate.column' ('", rep.col, "') is not present in the metadata table."), call. = FALSE)
      }

      say("[filter.samples] Re-running differential analyses (engine = '", engine,
          "', ", length(contrast.list), " contrast(s), which.data = '", diff.which, "')...")

      if (engine == "prolfqua") {
        dpo = DEprot::diff.analyses.prolfqua(DEprot.object = dpo,
                                             contrast.list = contrast.list,
                                             replicate.column = rep.col,
                                             linear.FC.th = linear.FC.th,
                                             linear.FC.unresp.range = linear.FC.unresp.range,
                                             FDR.th = padj.th,
                                             strategy = prolfqua.strategy,
                                             moderate.variance = moderate.value,
                                             which.data = diff.which,
                                             overwrite.analyses = TRUE)

      } else if (engine == "limma") {
        dpo = DEprot::diff.analyses.limma(DEprot.object = dpo,
                                          contrast.list = contrast.list,
                                          include.rep.model = include.rep.model.value,
                                          replicate.column = rep.col,
                                          linear.FC.th = linear.FC.th,
                                          linear.FC.unresp.range = linear.FC.unresp.range,
                                          padj.th = padj.th,
                                          padj.method = padj.method,
                                          fitting.method = fitting.method.value,
                                          which.data = diff.which,
                                          overwrite.analyses = TRUE)

      } else { # t.test / wilcoxon engine
        dpo = DEprot::diff.analyses(DEprot.object = dpo,
                                    contrast.list = contrast.list,
                                    replicate.column = rep.col,
                                    linear.FC.th = linear.FC.th,
                                    linear.FC.unresp.range = linear.FC.unresp.range,
                                    padj.th = padj.th,
                                    padj.method = padj.method,
                                    stat.test = stat.test.value,
                                    paired.test = paired.test.value,
                                    which.data = diff.which,
                                    overwrite.analyses = TRUE)
      }
    }

    say("[filter.samples] Done.")
    return(dpo)

  } # END function
