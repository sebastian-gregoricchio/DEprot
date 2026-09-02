#' @title .collect.assays
#'
#' @description Internal helper collecting the count matrices of a \code{DEprot} object and aligning
#'   them on the proteins and samples of the primary assay, as required by the Bioconductor classes.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param counts.type String, \code{"auto"} or one among \code{"raw"}, \code{"normalized"},
#'   \code{"randomized"}, \code{"imputed"}.
#' @param assays String \code{"all"} or a character vector of assay names.
#' @param verbose Logical.
#'
#' @return A list with the elements \code{assays} (named list of matrices, primary one first) and
#'   \code{reference} (string).
#'
#' @keywords internal

.collect.assays =
  function(DEprot.object,
           counts.type = "auto",
           assays = "all",
           verbose = TRUE) {

    ## slot corresponding to each assay name, in pipeline order
    slot.map = c(raw = "raw.counts",
                 normalized = "norm.counts",
                 randomized = "random.counts",
                 imputed = "imputed.counts")

    available = names(slot.map)[sapply(slot.map, function(x) {!is.null(methods::slot(DEprot.object, x))})]

    if (length(available) == 0) {
      stop("The object does not contain any count matrix.", call. = FALSE)
    }


    ### ------------------------------------------------------------------- assays to keep
    if (length(assays) == 1 && tolower(as.character(assays)[1]) == "all") {
      keep = available
    } else {
      assays = tolower(as.character(assays))

      if (!all(assays %in% names(slot.map))) {
        stop("Unknown assay(s): ", paste(setdiff(assays, names(slot.map)), collapse = ", "),
             ".\nPossible values: ", paste(names(slot.map), collapse = ", "), ", or 'all'.",
             call. = FALSE)
      }

      keep = intersect(available, assays)

      if (length(keep) == 0) {
        stop("None of the assays requested is available in this object.\nAvailable counts: ",
             paste(available, collapse = ", "), ".", call. = FALSE)
      }
    }


    ### ------------------------------------------------------------------- primary assay
    ## 'available' and 'keep' follow the pipeline order, hence the last one is the most processed
    if (counts.type == "auto") {
      reference = keep[length(keep)]
    } else {
      if (!(counts.type %in% available)) {
        stop("The '", counts.type, "' counts are not available in this object.\nAvailable counts: ",
             paste(available, collapse = ", "), ".", call. = FALSE)
      }
      reference = counts.type
      keep = intersect(available, unique(c(reference, keep)))
    }

    ## the primary assay must come first: it is the one returned by assay() without arguments
    keep = c(reference, setdiff(keep, reference))


    ### ------------------------------------------------------------------- alignment
    ref.mat = as.matrix(methods::slot(DEprot.object, unname(slot.map[reference])))
    features = rownames(ref.mat)
    samples  = colnames(ref.mat)

    assay.list = list()
    n.filled = c()

    for (i in keep) {
      mat = as.matrix(methods::slot(DEprot.object, unname(slot.map[i])))

      if (identical(rownames(mat), features) && identical(colnames(mat), samples)) {
        assay.list[[i]] = mat
        next
      }

      aligned = mat[match(features, rownames(mat)), match(samples, colnames(mat)), drop = FALSE]
      dimnames(aligned) = list(features, samples)
      assay.list[[i]] = aligned

      missing.features = sum(!(features %in% rownames(mat)))
      if (missing.features > 0) {n.filled[i] = missing.features}
    }

    if (isTRUE(verbose) && length(n.filled) > 0) {
      message("Assays re-indexed on the primary one; proteins filled with NA: ",
              paste0(names(n.filled), " (", n.filled, ")", collapse = ", "), ".")
    }

    return(list(assays = assay.list, reference = reference))
  }




#' @title .build.row.annotation
#'
#' @description Internal helper assembling the row annotation of the exported object from the
#'   \code{protein.info} slot and, when available, from the differential results.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param features Character vector of protein IDs (rownames of the primary assay).
#' @param add.protein.info Logical.
#' @param add.results Logical.
#' @param contrast.subset Numeric vector or NULL.
#' @param verbose Logical.
#'
#' @return A data.frame with \code{length(features)} rows.
#'
#' @keywords internal

.build.row.annotation =
  function(DEprot.object,
           features,
           add.protein.info = TRUE,
           add.results = TRUE,
           contrast.subset = NULL,
           verbose = TRUE) {

    row.data = data.frame(prot.id = features, stringsAsFactors = FALSE)
    rownames(row.data) = features


    ### ------------------------------------------------------------------- protein.info
    ## `.hasSlot` (and not `slotNames`) so that objects serialized before the introduction of the
    ## slot are handled gracefully instead of raising an error
    if (isTRUE(add.protein.info) && methods::.hasSlot(DEprot.object, "protein.info")) {
      info = methods::slot(DEprot.object, "protein.info")

      if (!is.null(info)) {
        info = as.data.frame(info)
        info = info[match(features, rownames(info)), , drop = FALSE]
        rownames(info) = features
        colnames(info) = make.unique(c(colnames(row.data), colnames(info)))[-seq_along(colnames(row.data))]
        row.data = cbind(row.data, info)
      }
    }


    ### ------------------------------------------------------------------- differential results
    if (isTRUE(add.results) && (methods::is(DEprot.object, "DEprot.analyses"))) {
      result.list = DEprot.object@analyses.result.list

      if (!is.null(contrast.subset)) {
        if (!all(contrast.subset %in% 1:length(result.list))) {
          stop("Not all the contrasts indicated in the subset are present in the 'analyses.result.list' of the object provided.",
               call. = FALSE)
        }
        contrasts = contrast.subset
      } else {
        contrasts = 1:length(result.list)
      }

      for (i in contrasts) {
        res = as.data.frame(result.list[[i]]$results)

        if (is.null(res) || nrow(res) == 0) {next}

        rownames(res) = as.character(res$prot.id)
        res = res[match(features, rownames(res)), setdiff(colnames(res), "prot.id"), drop = FALSE]
        rownames(res) = features

        ## prefix with the contrast name to avoid collisions between contrasts
        colnames(res) = paste(make.names(names(result.list)[i]), colnames(res), sep = ".")
        row.data = cbind(row.data, res)
      }

      if (isTRUE(verbose) && length(contrasts) > 0) {
        message("Differential results added to the row annotation for ", length(contrasts),
                " contrast(s).")
      }
    }

    return(row.data)
  }




#' @title .build.col.annotation
#'
#' @description Internal helper assembling the column annotation (sample metadata) of the exported
#'   object, reordered to match the columns of the primary assay.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param samples Character vector of sample names (colnames of the primary assay).
#'
#' @return A data.frame with \code{length(samples)} rows.
#'
#' @keywords internal

.build.col.annotation =
  function(DEprot.object,
           samples) {

    meta = as.data.frame(DEprot.object@metadata)

    if ("column.id" %in% colnames(meta)) {
      meta = meta[match(samples, as.character(meta$column.id)), , drop = FALSE]
    } else {
      meta = meta[match(samples, rownames(meta)), , drop = FALSE]
    }

    rownames(meta) = samples

    return(meta)
  }




#' @title .build.export.metadata
#'
#' @description Internal helper collecting the experimental parameters that have no equivalent in the
#'   destination classes, stored in the metadata list of the exported object.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param keep.object Logical, whether the original object should be stored as well.
#'
#' @return A named list.
#'
#' @keywords internal

.build.export.metadata =
  function(DEprot.object,
           keep.object = FALSE) {

    obj.metadata =
      list(exported.by = "DEprot",
           DEprot.version = as.character(tryCatch(utils::packageVersion("DEprot"), error = function(e) {NA})),
           export.date = as.character(Sys.Date()),
           log.base = DEprot.object@log.base,
           log.transformed = DEprot.object@log.transformed,
           normalized = DEprot.object@normalized,
           normalization.method = DEprot.object@normalization.method,
           randomized = DEprot.object@randomized,
           randomization.method = DEprot.object@randomization.method,
           imputed = DEprot.object@imputed,
           imputation.method = DEprot.object@imputation.method,
           contrasts = DEprot.object@contrasts,
           differential.analyses.params = DEprot.object@differential.analyses.params)

    if (isTRUE(keep.object)) {
      obj.metadata$DEprot.object = DEprot.object
    }

    return(obj.metadata)
  }




#' @title .build.summarized.experiment
#'
#' @description Internal helper building a \code{SummarizedExperiment} from the collected assays and
#'   annotations.
#'
#' @param assay.list Named list of matrices sharing the same dimnames.
#' @param row.data A data.frame.
#' @param col.data A data.frame.
#' @param obj.metadata A named list.
#' @param install.missing One among \code{"ask"}, \code{"always"}, \code{"never"}.
#'
#' @return A \code{SummarizedExperiment} object.
#'
#' @keywords internal

.build.summarized.experiment =
  function(assay.list,
           row.data,
           col.data,
           obj.metadata,
           install.missing = "ask") {

    .require.package(package = "SummarizedExperiment",
                     repo = "Bioconductor",
                     install.missing = install.missing,
                     reason = "export a DEprot object as a SummarizedExperiment")

    se = SummarizedExperiment::SummarizedExperiment(assays = assay.list,
                                                    rowData = row.data,
                                                    colData = col.data,
                                                    metadata = obj.metadata)

    return(se)
  }




#' @title .build.qfeatures
#'
#' @description Internal helper building a \code{QFeatures} object containing a single, protein-level
#'   assay set.
#'
#' @param assay.list Named list of matrices sharing the same dimnames.
#' @param row.data A data.frame.
#' @param col.data A data.frame.
#' @param obj.metadata A named list.
#' @param assay.name String, name of the assay set.
#' @param install.missing One among \code{"ask"}, \code{"always"}, \code{"never"}.
#'
#' @return A \code{QFeatures} object.
#'
#' @keywords internal

.build.qfeatures =
  function(assay.list,
           row.data,
           col.data,
           obj.metadata,
           assay.name = "proteins",
           install.missing = "ask") {

    .require.package(package = "QFeatures",
                     repo = "Bioconductor",
                     install.missing = install.missing,
                     reason = "export a DEprot object as a QFeatures object")

    se = .build.summarized.experiment(assay.list = assay.list,
                                      row.data = row.data,
                                      col.data = col.data,
                                      obj.metadata = obj.metadata,
                                      install.missing = install.missing)

    experiment.list = list(se)
    names(experiment.list) = as.character(assay.name)[1]

    qf = QFeatures::QFeatures(experiments = experiment.list,
                              colData = col.data,
                              metadata = obj.metadata)

    return(qf)
  }




#' @title .build.msnset
#'
#' @description Internal helper building an \code{MSnSet}. Since this class holds a single expression
#'   matrix, only the primary assay is exported.
#'
#' @param assay.list Named list of matrices sharing the same dimnames.
#' @param row.data A data.frame.
#' @param col.data A data.frame.
#' @param obj.metadata A named list.
#' @param install.missing One among \code{"ask"}, \code{"always"}, \code{"never"}.
#' @param verbose Logical.
#'
#' @return An \code{MSnSet} object.
#'
#' @keywords internal

.build.msnset =
  function(assay.list,
           row.data,
           col.data,
           obj.metadata,
           install.missing = "ask",
           verbose = TRUE) {

    .require.package(package = "MSnbase",
                     repo = "Bioconductor",
                     install.missing = install.missing,
                     reason = "export a DEprot object as an MSnSet")

    if (isTRUE(verbose) && length(assay.list) > 1) {
      message("An MSnSet holds a single matrix: only the '", names(assay.list)[1],
              "' counts have been exported.")
    }

    ms = MSnbase::MSnSet(exprs = assay.list[[1]],
                         fData = row.data,
                         pData = col.data)

    ## the parameters that MSnSet cannot store natively are kept in the experiment description
    if ("other" %in% methods::slotNames(ms@experimentData)) {
      ms@experimentData@other = obj.metadata
    }

    return(ms)
  }




#' @title .build.elist
#'
#' @description Internal helper building a \code{limma::EList}. Since this class holds a single
#'   expression matrix, only the primary assay is exported.
#'
#' @param assay.list Named list of matrices sharing the same dimnames.
#' @param row.data A data.frame.
#' @param col.data A data.frame.
#' @param verbose Logical.
#'
#' @return An \code{EList} object.
#'
#' @keywords internal

.build.elist =
  function(assay.list,
           row.data,
           col.data,
           verbose = TRUE) {

    if (isTRUE(verbose) && length(assay.list) > 1) {
      message("An EList holds a single matrix: only the '", names(assay.list)[1],
              "' counts have been exported.")
    }

    el = methods::new("EList",
                      list(E = assay.list[[1]],
                           genes = row.data,
                           targets = col.data))

    return(el)
  }
