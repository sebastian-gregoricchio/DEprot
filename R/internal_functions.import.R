#' @title .require.package
#'
#' @description Internal helper checking for an optional (\code{Suggests}) dependency and,
#'   if absent, offering to install it. It never installs silently: in non-interactive
#'   sessions it always fails with an informative message reporting the install command.
#'
#' @param package String, name of the package.
#' @param repo String indicating where to get the package from: \code{"CRAN"},
#'   \code{"Bioconductor"}, or a GitHub slug such as \code{"tvpham/iq"}.
#' @param install.missing One among \code{"ask"} (default), \code{"always"}, \code{"never"}.
#' @param reason Optional string completing the sentence "<package> is required to <reason>".
#'
#' @return Invisibly \code{TRUE} if the package is available after the call, otherwise it stops.
#'
#' @keywords internal

.require.package =
  function(package,
           repo = "CRAN",
           install.missing = "ask",
           reason = NULL) {

    ### already available -> nothing to do
    if (requireNamespace(package, quietly = TRUE)) {
      return(invisible(TRUE))
    }

    install.missing = match.arg(tolower(as.character(install.missing)[1]),
                                choices = c("ask", "always", "never"))

    ### define the install command (shown to the user and used if accepted)
    if (tolower(repo) == "cran") {
      cmd.string = paste0("install.packages('", package, "')")
      installer  = function() {utils::install.packages(package)}

    } else if (tolower(repo) %in% c("bioc", "bioconductor")) {
      cmd.string = paste0("BiocManager::install('", package, "')")
      installer  = function() {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          utils::install.packages("BiocManager")
        }
        BiocManager::install(package, update = FALSE, ask = FALSE)
      }

    } else {
      ## anything else is interpreted as a GitHub slug ("user/repo")
      cmd.string = paste0("remotes::install_github('", repo, "')")
      installer  = function() {
        if (!requireNamespace("remotes", quietly = TRUE)) {
          utils::install.packages("remotes")
        }
        remotes::install_github(repo, upgrade = "never", build_vignettes = FALSE)
      }
    }

    base.msg = paste0("The package '", package, "' is required",
                      if (is.null(reason)) "" else paste0(" to ", reason), ".")

    ### never install without an explicit, interactive confirmation
    if (install.missing == "never" || !interactive()) {
      stop(base.msg, "\nInstall it with:  ", cmd.string,
           "\n(or re-run in an interactive session with `install.missing = \"always\"`).",
           call. = FALSE)
    }

    if (install.missing == "ask") {
      message(base.msg)
      answer = readline(paste0("Do you want to install '", package, "' now? [yes/no] "))
      if (!(tolower(trimws(answer)) %in% c("y", "yes", "yeah", "yep", "yo", "ok", "1"))) {
        stop("Installation declined. ", base.msg,
             "\nYou can install it later with:  ", cmd.string, call. = FALSE)
      }
    }

    ### install and re-check
    message("Installing '", package, "' ...")
    try(installer(), silent = FALSE)

    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Installation of '", package, "' failed. Please install it manually:  ",
           cmd.string, call. = FALSE)
    }

    return(invisible(TRUE))
  }




#' @title .read.table.any
#'
#' @description Internal reader dispatching on the file extension: delimited text is read with
#'   \code{data.table::fread}, parquet with \code{nanoparquet} (or \code{arrow} if already installed).
#'
#' @param file String, path to the file.
#' @param install.missing One among \code{"ask"}, \code{"always"}, \code{"never"}.
#'
#' @return A \code{data.frame}.
#'
#' @keywords internal

.read.table.any =
  function(file,
           install.missing = "ask") {

    if (!file.exists(file)) {
      stop("File not found: ", file, call. = FALSE)
    }

    if (grepl("\\.parquet$", file, ignore.case = TRUE)) {
      if (requireNamespace("nanoparquet", quietly = TRUE)) {
        return(as.data.frame(nanoparquet::read_parquet(file)))
      } else if (requireNamespace("arrow", quietly = TRUE)) {
        return(as.data.frame(arrow::read_parquet(file)))
      } else {
        .require.package(package = "nanoparquet",
                         repo = "CRAN",
                         install.missing = install.missing,
                         reason = "read parquet reports (e.g. the DIA-NN >= 1.9 'report.parquet')")
        return(as.data.frame(nanoparquet::read_parquet(file)))
      }
    }

    return(data.table::fread(file, data.table = FALSE, check.names = FALSE))
  }




#' @title .clean.run.names
#'
#' @description Internal helper converting raw-file paths into clean sample names: removes
#'   Spectronaut-style \code{"[1] "} prefixes, directory paths and raw-file extensions.
#'
#' @param x Character vector.
#'
#' @return A character vector.
#'
#' @keywords internal

.clean.run.names =
  function(x) {
    x = gsub("^\\[[0-9]+\\]\\s*", "", x)                 # Spectronaut "[1] sample"
    x = gsub("\\\\", "/", x)                             # windows separators
    x = basename(x)
    x = sub("\\.(raw|mzML|mzXML|d|wiff|dia|htrms|timsTOF)$", "", x, ignore.case = TRUE)
    return(trimws(x))
  }




#' @title .pick.column
#'
#' @description Internal helper selecting a column: returns \code{user} when the user provided one
#'   (checking that it exists), otherwise the first available among \code{candidates}.
#'
#' @param user String provided by the user, or \code{"auto"}.
#' @param candidates Character vector of candidate column names, in order of preference.
#' @param columns Character vector of the column names actually available.
#' @param what String used in the error message.
#'
#' @return A string, or \code{NA_character_} if nothing matches.
#'
#' @keywords internal

.pick.column =
  function(user,
           candidates,
           columns,
           what = "column") {

    if (!identical(tolower(as.character(user)[1]), "auto")) {
      if (!(user %in% columns)) {
        stop("The ", what, " '", user, "' is not present in the file.\nAvailable columns: ",
             paste(columns, collapse = ", "), call. = FALSE)
      }
      return(user)
    }

    hit = intersect(candidates, columns)
    if (length(hit) == 0) {
      return(NA_character_)
    }
    return(hit[1])
  }




#' @title .long.to.matrix
#'
#' @description Internal helper pivoting a long table (id / sample / quantity) into a
#'   proteins-by-samples numeric matrix. Duplicated id-sample pairs are aggregated with \code{FUN}.
#'   Zeros are converted to \code{NA}.
#'
#' @param df A data.frame.
#' @param id.col String, column holding the protein/feature ID.
#' @param sample.col String, column holding the run/sample ID.
#' @param quantity.col String, column holding the intensity.
#' @param FUN Function used to aggregate duplicated entries. Default: sum ignoring NAs.
#' @param zero.to.na Logical, whether zeros must be converted to \code{NA}. \code{TRUE} (default) is
#'   correct for linear intensities, where a zero means "not measured". Set it to \code{FALSE} for
#'   already log-transformed values, where zero is a legitimate abundance.
#'
#' @return A numeric matrix.
#'
#' @keywords internal

.long.to.matrix =
  function(df,
           id.col,
           sample.col,
           quantity.col,
           FUN = function(x) {sum(x, na.rm = TRUE)},
           zero.to.na = TRUE) {

    df = data.frame(id = as.character(df[[id.col]]),
                    sample = as.character(df[[sample.col]]),
                    value = suppressWarnings(as.numeric(df[[quantity.col]])),
                    stringsAsFactors = FALSE)

    if (isTRUE(zero.to.na)) {
      df$value[which(df$value == 0)] = NA
    }

    df = df[!is.na(df$id) & df$id != "" & !is.na(df$sample) & is.finite(df$value), , drop = FALSE]

    if (nrow(df) == 0) {
      stop("No finite values found in the column '", quantity.col, "'.", call. = FALSE)
    }

    ## data.table is much faster than aggregate() on multi-million-row reports
    dt = data.table::as.data.table(df)
    agg = as.data.frame(dt[, list(value = FUN(value)), by = list(id, sample)])

    wide = reshape2::dcast(agg, id ~ sample, value.var = "value")
    mat = as.matrix(wide[, -1, drop = FALSE])
    rownames(mat) = wide[, 1]

    if (isTRUE(zero.to.na)) {
      mat[which(mat == 0)] = NA
    }

    return(mat)
  }




#' @title .wide.to.matrix
#'
#' @description Internal helper extracting a proteins-by-samples matrix from an already wide report
#'   (DIA-NN matrices, MaxQuant proteinGroups, FragPipe combined_protein, Spectronaut pivot, ...).
#'
#' @param df A data.frame.
#' @param id.col String, column holding the protein/feature ID.
#' @param quant.cols Named character vector: values are the intensity columns in \code{df},
#'   names are the sample names to use in the output.
#'
#' @return A numeric matrix.
#'
#' @keywords internal

.wide.to.matrix =
  function(df,
           id.col,
           quant.cols) {

    if (nrow(df) == 0) {
      stop("No entries left in the table: the contaminant/decoy filters removed everything.",
           call. = FALSE)
    }

    ids = trimws(as.character(df[[id.col]]))

    keep = !is.na(ids) & ids != ""
    if (any(!keep)) {
      warning(sum(!keep), " row(s) with a missing '", id.col, "' have been discarded.", call. = FALSE)
      df = df[keep, , drop = FALSE]
      ids = ids[keep]
    }

    mat = vapply(X = unname(quant.cols),
                 FUN = function(i) {suppressWarnings(as.numeric(as.character(df[[i]])))},
                 FUN.VALUE = numeric(nrow(df)))

    mat = matrix(mat, nrow = nrow(df), ncol = length(quant.cols))
    colnames(mat) = names(quant.cols)
    rownames(mat) = make.unique(ids, sep = ".")
    mat[which(mat == 0)] = NA

    return(mat)
  }




#' @title .maxlfq.rollup
#'
#' @description Internal helper performing the MaxLFQ protein summarization of a long-format
#'   report through the optional package \code{iq}.
#'
#' @param tb A data.frame in long format.
#' @param primary.id String, protein-group column.
#' @param secondary.id Character vector of the columns defining a precursor/fragment.
#' @param sample.col String, run/sample column.
#' @param quantity.col String, intensity column.
#' @param median.normalization Logical, whether \code{iq} should median-normalize. Default \code{FALSE}.
#' @param install.missing One among \code{"ask"}, \code{"always"}, \code{"never"}.
#' @param verbose Logical.
#'
#' @return A numeric matrix of log2 estimates (proteins x samples).
#'
#' @keywords internal

.maxlfq.rollup =
  function(tb,
           primary.id,
           secondary.id,
           sample.col,
           quantity.col,
           median.normalization = FALSE,
           install.missing = "ask",
           verbose = TRUE) {

    .require.package(package = "iq",
                     repo = "CRAN",
                     install.missing = install.missing,
                     reason = paste0("perform the MaxLFQ summarization (`summarization = \"maxlfq\"`).\n",
                                     "Alternatively, use `summarization = \"none\"` to take the protein-level ",
                                     "quantities already present in the report, or \"sum\"/\"median\""))

    secondary.id = intersect(secondary.id, colnames(tb))
    if (length(secondary.id) == 0) {
      stop("None of the expected precursor/fragment-level columns was found: MaxLFQ summarization ",
           "is not possible on this file. Use `summarization = \"none\"` instead.", call. = FALSE)
    }

    if (is.na(quantity.col) || !(quantity.col %in% colnames(tb))) {
      stop("The intensity column required for the MaxLFQ summarization was not found in the report.",
           call. = FALSE)
    }

    keep.cols = unique(c(primary.id, secondary.id, sample.col, quantity.col))
    sub = tb[, keep.cols, drop = FALSE]
    sub[[quantity.col]] = suppressWarnings(as.numeric(sub[[quantity.col]]))
    sub = sub[!is.na(sub[[quantity.col]]) & sub[[quantity.col]] > 0, , drop = FALSE]

    if (isTRUE(verbose)) {
      message("Running MaxLFQ (iq) on ", format(nrow(sub), big.mark = ","),
              " rows using '", quantity.col, "' [", paste(secondary.id, collapse = " + "), "] ...")
    }

    norm.data = iq::preprocess(sub,
                               primary_id = primary.id,
                               secondary_id = secondary.id,
                               sample_id = sample.col,
                               intensity_col = quantity.col,
                               median_normalization = median.normalization,
                               log2_intensity_cutoff = NULL,
                               pdf_out = NULL)

    ## fast implementation when available, documented three-step pipeline as fallback
    res = tryCatch(expr = iq::fast_MaxLFQ(norm.data),
                   error = function(e) {
                     iq::create_protein_table(iq::create_protein_list(norm.data), method = "maxLFQ")
                   })

    mat = as.matrix(res$estimate)   # log2 space

    return(mat)
  }




#' @title .finalize.import
#'
#' @description Internal helper shared by \link{import.external} and \link{import.msstats}: cleans the
#'   sample names, reconciles the counts with the metadata, resolves the \code{data.type} and builds
#'   the \code{DEprot} object through \link{load.counts2}.
#'
#' @param imported List returned by one of the internal readers, with the elements \code{counts},
#'   \code{log.base}, \code{data.type}, \code{normalization.method}, \code{imputation.method},
#'   \code{quantity}, \code{feature} and (optionally) \code{annotation}.
#' @param metadata A data.frame, a path to a table, or \code{NULL}.
#' @param column.id String indicating the metadata column matching the sample names.
#' @param data.type String, \code{"auto"} or one of the \code{load.counts2} types.
#' @param log.base Numeric or \code{NULL}: when provided it overrides the value proposed by the reader.
#' @param clean.sample.names Logical.
#' @param subset.to.metadata Logical.
#' @param verbose Logical.
#'
#' @return A \code{DEprot} object.
#'
#' @keywords internal

.finalize.import =
  function(imported,
           metadata = NULL,
           column.id = "column.id",
           data.type = "auto",
           log.base = NULL,
           clean.sample.names = TRUE,
           subset.to.metadata = TRUE,
           verbose = TRUE) {

    cnt = imported$counts

    if (is.null(cnt) || nrow(cnt) == 0 || ncol(cnt) == 0) {
      stop("No quantitative values could be extracted. Check the import parameters.", call. = FALSE)
    }


    ### ------------------------------------------------------------------- sample names
    if (isTRUE(clean.sample.names)) {
      new.names = .clean.run.names(colnames(cnt))

      ## keep the original names if the cleaning would collapse distinct samples
      if (!any(duplicated(new.names))) {
        ## the annotation table must follow the renaming
        if (!is.null(imported$annotation)) {
          idx = match(imported$annotation[[column.id]], colnames(cnt))
          imported$annotation[[column.id]][!is.na(idx)] = new.names[idx[!is.na(idx)]]
        }
        colnames(cnt) = new.names
      } else if (isTRUE(verbose)) {
        message("Sample names were left untouched: cleaning them would generate duplicates.")
      }
    }

    if (any(duplicated(colnames(cnt)))) {
      warning("Duplicated sample names: `make.unique` has been applied.", call. = FALSE)
      colnames(cnt) = make.unique(colnames(cnt), sep = ".")
    }


    ### ------------------------------------------------------------------- metadata
    if (is.null(metadata)) {
      if (!is.null(imported$annotation)) {
        ## annotation recovered from the imported object (conditions, replicates, mixtures, ...)
        metadata = as.data.frame(imported$annotation)
        metadata = metadata[match(colnames(cnt), metadata[[column.id]]), , drop = FALSE]
        rownames(metadata) = NULL
        if (isTRUE(verbose)) {
          message("No 'metadata' provided: it has been reconstructed from the imported object (",
                  paste(setdiff(colnames(metadata), column.id), collapse = ", "), ").")
        }
      } else {
        metadata = data.frame(colnames(cnt), colnames(cnt), stringsAsFactors = FALSE)
        colnames(metadata) = c(column.id, "sample.id")
        metadata = metadata[, unique(colnames(metadata)), drop = FALSE]
        warning("No 'metadata' provided: a minimal table has been generated from the sample names.\n",
                "Add the experimental groups before performing any differential analysis.",
                call. = FALSE)
      }

    } else if (is.character(metadata)) {
      metadata = data.table::fread(metadata, data.table = FALSE)
    }

    metadata = as.data.frame(metadata)

    if (!(column.id %in% colnames(metadata))) {
      stop("The column '", column.id, "' is not present in the metadata table.\nAvailable columns: ",
           paste(colnames(metadata), collapse = ", "), call. = FALSE)
    }

    meta.ids = as.character(metadata[, column.id])
    only.meta   = setdiff(meta.ids, colnames(cnt))
    only.counts = setdiff(colnames(cnt), meta.ids)

    if (length(only.meta) > 0) {
      stop("Some samples listed in the metadata are absent from the imported data:\n  ",
           paste(only.meta, collapse = ", "),
           "\n\nSample names extracted from the data:\n  ",
           paste(colnames(cnt), collapse = ", "),
           "\n\nUse `clean.sample.names = FALSE` if the metadata contains the full raw-file names.",
           call. = FALSE)
    }

    if (length(only.counts) > 0) {
      if (isTRUE(subset.to.metadata)) {
        if (isTRUE(verbose)) {
          message(length(only.counts), " sample(s) present in the data but absent from the metadata ",
                  "have been discarded: ", paste(only.counts, collapse = ", "), ".")
        }
      } else {
        stop("Some samples of the imported data are absent from the metadata:\n  ",
             paste(only.counts, collapse = ", "),
             "\nSet `subset.to.metadata = TRUE` to discard them.", call. = FALSE)
      }
    }

    cnt = cnt[, meta.ids, drop = FALSE]

    ## drop features that became all-NA after the sample subsetting
    keep = rowSums(!is.na(cnt)) > 0
    if (any(!keep)) {
      cnt = cnt[keep, , drop = FALSE]
    }


    ### ------------------------------------------------------------------- data.type and methods
    if (identical(tolower(as.character(data.type)[1]), "auto")) {
      data.type = imported$data.type
    }

    is.normalized = tolower(data.type) %in% c("n", "nor", "norm", "normalized", "normalised")
    is.imputed    = tolower(data.type) %in% c("i", "im", "imp", "imputed")

    normalization.method = if (is.normalized) imported$normalization.method else NA
    imputation.method    = if (is.imputed) {
      if (is.null(imported$imputation.method)) NA else imported$imputation.method
    } else {
      NA
    }

    if (isTRUE(verbose)) {
      message("Imported ", format(nrow(cnt), big.mark = ","), " ", imported$feature, " x ",
              ncol(cnt), " samples (quantity: '", imported$quantity,
              "'; loaded as '", data.type, "').")
    }


    ### ------------------------------------------------------------------- build the object
    DEprot.object =
      load.counts2(counts = cnt,
                   metadata = metadata,
                   data.type = data.type,
                   log.base = if (is.null(log.base)) imported$log.base else log.base,
                   normalization.method = normalization.method,
                   randomization.method = NA,
                   imputation.method = imputation.method,
                   column.id = column.id)

    return(DEprot.object)
  }
