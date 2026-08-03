#' @title import.external
#'
#' @description Imports the output of the most common proteomics search/quantification tools
#'   (DIA-NN, Spectronaut, FragPipe, MaxQuant, Proteome Discoverer) and converts it into a
#'   \code{DEprot} object. The function only reads, filters and reshapes the report: the object
#'   itself is then built by \link{load.counts2}, so log2 standardization, removal of all-NA rows,
#'   boxplot generation and slot filling behave exactly as in a manual load.
#'
#' @param file String, path to the report file. Delimited text (\code{.tsv}, \code{.csv}, \code{.txt})
#'   is read with \code{data.table::fread}; \code{.parquet} requires the optional package
#'   \code{nanoparquet} (or \code{arrow}, if already installed).
#' @param metadata A data.frame (or a path to a table) containing at least the column indicated by
#'   \code{column.id}, whose values must match the sample names extracted from the report.
#'   If \code{NULL} (default) a minimal metadata table is generated from the sample names, with a warning.
#' @param source String indicating the tool that generated the file. One among: \code{"auto"} (default),
#'   \code{"diann"}, \code{"diann.matrix"}, \code{"spectronaut"}, \code{"spectronaut.pivot"},
#'   \code{"fragpipe"}, \code{"maxquant"}, \code{"proteome.discoverer"}, \code{"generic"}.
#' @param quantity String indicating the intensity column (long reports) or the column suffix/prefix
#'   (wide reports). Default \code{"auto"} uses the tool-specific recommended one, see 'Details'.
#' @param protein.id String indicating the column to use as protein identifier (row names of the
#'   count matrix). Default \code{"auto"}.
#' @param sample.id String indicating the run/sample column of long reports. Default \code{"auto"}.
#' @param summarization String defining how precursor/fragment-level long reports are rolled up to the
#'   protein level: \code{"maxlfq"} (default, requires the optional package \code{iq}), \code{"sum"},
#'   \code{"median"}, or \code{"none"} (use a protein-level column already present in the report,
#'   e.g. \code{PG.MaxLFQ} or \code{PG.Quantity}). Ignored for wide reports.
#' @param q.value Numeric, run-specific precursor q-value cutoff applied to long reports.
#'   Default \code{0.01}; use \code{NULL} to skip.
#' @param pg.q.value Numeric, protein-group q-value cutoff applied to long reports.
#'   Default \code{0.01}; use \code{NULL} to skip.
#' @param remove.contaminants Logical, whether to drop contaminant/decoy entries. Default \code{TRUE}.
#' @param contaminant.pattern String, regular expression identifying contaminants/decoys.
#'   Default \code{NULL} uses a tool-aware pattern.
#' @param clean.sample.names Logical, whether to strip paths and raw-file extensions from the sample
#'   names. Default \code{TRUE}.
#' @param subset.to.metadata Logical, whether to keep only the samples listed in the metadata.
#'   Default \code{TRUE}.
#' @param column.id String indicating the metadata column matching the sample names.
#'   Default \code{"column.id"}.
#' @param data.type String, one among \code{"auto"} (default), \code{"raw"}, \code{"normalized"},
#'   \code{"randomized"}, \code{"imputed"}. When \code{"auto"} it is inferred from the quantity used.
#' @param id.column String, used only with \code{source = "generic"}: the column holding the protein IDs.
#' @param sample.columns Character vector, used only with \code{source = "generic"}: the intensity
#'   columns. If named, the names are used as sample names.
#' @param install.missing String defining the behaviour when an optional package is missing:
#'   \code{"ask"} (default, prompts in interactive sessions), \code{"always"}, \code{"never"}.
#' @param verbose Logical, whether to print progress messages. Default \code{TRUE}.
#'
#' @details
#' No new hard dependency is introduced: all delimited reports are parsed with packages already
#' imported by DEprot. Two optional (\code{Suggests}) packages are used, and only when needed:
#' \describe{
#'   \item{\code{iq} (CRAN)}{MaxLFQ summarization of long precursor-level reports
#'     (\code{summarization = "maxlfq"}).}
#'   \item{\code{nanoparquet} (CRAN)}{reading \code{.parquet} reports (DIA-NN >= 1.9);
#'     \code{arrow} is used instead when already installed.}
#' }
#' If a required optional package is missing, the user is prompted to install it
#' (see \code{install.missing}). Nothing is ever installed in a non-interactive session.
#'
#' Default quantity column per source:
#' \tabular{lll}{
#'   \strong{source} \tab \strong{default quantity} \tab \strong{inferred data.type} \cr
#'   \code{diann} \tab \code{Precursor.Normalised} (maxlfq) / \code{PG.MaxLFQ} (none) \tab normalized \cr
#'   \code{diann.matrix} \tab all non-annotation columns \tab normalized \cr
#'   \code{spectronaut} \tab \code{F.PeakArea} (maxlfq) / \code{PG.Quantity} (none) \tab normalized \cr
#'   \code{spectronaut.pivot} \tab \code{PG.Quantity} \tab normalized \cr
#'   \code{fragpipe} \tab \code{"MaxLFQ Intensity"} \tab normalized \cr
#'   \code{maxquant} \tab \code{"LFQ intensity"} \tab normalized \cr
#'   \code{proteome.discoverer} \tab \code{"Abundances (Normalized)"} \tab normalized \cr
#' }
#'
#' Intensities are kept on a linear scale (zeros converted to \code{NA}) except for
#' \code{summarization = "maxlfq"}, which returns log2 estimates from \code{iq}. The corresponding
#' \code{log.base} is passed to \link{load.counts2} automatically.
#'
#' @return A \code{DEprot} object (S4 vector).
#'
#' @import dplyr
#' @import methods
#' @importFrom data.table fread as.data.table
#' @importFrom reshape2 dcast
#' @importFrom stats median
#' @importFrom utils head install.packages
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' \dontrun{
#' # DIA-NN protein-group matrix (no optional dependency required)
#' dpo <- import.external(file = "report.pg_matrix.tsv",
#'                        metadata = sample.config,
#'                        source = "diann.matrix")
#'
#' # DIA-NN main report, MaxLFQ roll-up from precursors (uses `iq`)
#' dpo <- import.external(file = "report.parquet",
#'                        metadata = sample.config,
#'                        source = "diann",
#'                        summarization = "maxlfq")
#'
#' # FragPipe
#' dpo <- import.external(file = "combined_protein.tsv",
#'                        metadata = sample.config,
#'                        source = "fragpipe",
#'                        quantity = "MaxLFQ Intensity")
#'
#' # any other wide table
#' dpo <- import.external(file = "my_table.tsv",
#'                        metadata = sample.config,
#'                        source = "generic",
#'                        id.column = "Accession",
#'                        sample.columns = c(WT_1 = "int_1", WT_2 = "int_2"))
#' }
#'
#' @seealso \link{load.counts2}
#'
#' @export import.external

import.external =
  function(file,
           metadata = NULL,
           source = "auto",
           quantity = "auto",
           protein.id = "auto",
           sample.id = "auto",
           summarization = "maxlfq",
           q.value = 0.01,
           pg.q.value = 0.01,
           remove.contaminants = TRUE,
           contaminant.pattern = NULL,
           clean.sample.names = TRUE,
           subset.to.metadata = TRUE,
           column.id = "column.id",
           data.type = "auto",
           id.column = NULL,
           sample.columns = NULL,
           install.missing = "ask",
           verbose = TRUE) {

    ### ------------------------------------------------------------------- checks
    if (missing(file)) {
      stop("Provide 'file': the path to the report generated by the external tool.", call. = FALSE)
    }

    source = match.arg(tolower(as.character(source)[1]),
                       choices = c("auto", "diann", "diann.matrix", "spectronaut",
                                   "spectronaut.pivot", "fragpipe", "maxquant",
                                   "proteome.discoverer", "generic"))

    summarization = match.arg(tolower(as.character(summarization)[1]),
                              choices = c("maxlfq", "sum", "median", "none"))

    install.missing = match.arg(tolower(as.character(install.missing)[1]),
                                choices = c("ask", "always", "never"))


    ### ------------------------------------------------------------------- read
    if (isTRUE(verbose)) {message("Reading '", basename(file), "' ...")}
    tb = .read.table.any(file = file, install.missing = install.missing)

    if (source == "auto") {
      source = .detect.source(tb)
      if (isTRUE(verbose)) {message("Detected source: '", source, "'.")}
    }


    ### ------------------------------------------------------------------- dispatch
    imported =
      switch(source,
             "diann" = .read.diann(tb = tb, quantity = quantity, protein.id = protein.id,
                                   sample.id = sample.id, summarization = summarization,
                                   q.value = q.value, pg.q.value = pg.q.value,
                                   remove.contaminants = remove.contaminants,
                                   contaminant.pattern = contaminant.pattern,
                                   install.missing = install.missing, verbose = verbose),

             "diann.matrix" = .read.diann.matrix(tb = tb, protein.id = protein.id,
                                                 remove.contaminants = remove.contaminants,
                                                 contaminant.pattern = contaminant.pattern),

             "spectronaut" = .read.spectronaut(tb = tb, quantity = quantity, protein.id = protein.id,
                                               sample.id = sample.id, summarization = summarization,
                                               q.value = q.value, pg.q.value = pg.q.value,
                                               remove.contaminants = remove.contaminants,
                                               contaminant.pattern = contaminant.pattern,
                                               install.missing = install.missing, verbose = verbose),

             "spectronaut.pivot" = .read.spectronaut.pivot(tb = tb, quantity = quantity,
                                                           protein.id = protein.id,
                                                           remove.contaminants = remove.contaminants,
                                                           contaminant.pattern = contaminant.pattern),

             "fragpipe" = .read.fragpipe(tb = tb, quantity = quantity, protein.id = protein.id,
                                         remove.contaminants = remove.contaminants,
                                         contaminant.pattern = contaminant.pattern),

             "maxquant" = .read.maxquant(tb = tb, quantity = quantity, protein.id = protein.id,
                                         remove.contaminants = remove.contaminants,
                                         contaminant.pattern = contaminant.pattern,
                                         verbose = verbose),

             "proteome.discoverer" = .read.proteome.discoverer(tb = tb, quantity = quantity,
                                                               protein.id = protein.id,
                                                               remove.contaminants = remove.contaminants,
                                                               contaminant.pattern = contaminant.pattern),

             "generic" = .read.generic(tb = tb, id.column = id.column, sample.columns = sample.columns))

    ### ------------------------------------------------------------------- build the object
    DEprot.object =
      .finalize.import(imported = imported,
                       metadata = metadata,
                       column.id = column.id,
                       data.type = data.type,
                       clean.sample.names = clean.sample.names,
                       subset.to.metadata = subset.to.metadata,
                       verbose = verbose)

    return(DEprot.object)
  } # END function




# =====================================================================================
# Internal source-specific readers.
# Each returns: list(counts, log.base, data.type, normalization.method, quantity, feature)
# =====================================================================================


#' @title .detect.source
#' @description Internal helper guessing the tool of origin from the column names of a report.
#' @param tb A data.frame.
#' @return A string.
#' @keywords internal

.detect.source =
  function(tb) {

    cn = colnames(tb)

    ## DIA-NN: long reports always carry a run column, matrices never do
    if ("Protein.Group" %in% cn || "Protein.Ids" %in% cn) {
      if (any(c("Run", "File.Name") %in% cn)) {
        return("diann")
      } else {
        return("diann.matrix")
      }
    }

    ## Spectronaut
    if (any(grepl("^(R|PG|EG|FG|F|PEP)\\.", cn))) {
      if (any(c("R.FileName", "R.Raw.File.Name", "R.Label", "R.Condition") %in% cn)) {
        return("spectronaut")
      } else {
        return("spectronaut.pivot")
      }
    }

    ## MaxQuant
    if (any(c("Majority protein IDs", "Protein IDs") %in% cn)) {
      return("maxquant")
    }

    ## FragPipe
    if (all(c("Protein ID", "Entry Name") %in% cn) || "Protein Probability" %in% cn) {
      return("fragpipe")
    }

    ## Proteome Discoverer
    if (any(c("Accession", "Master.Protein.Accessions") %in% cn) && any(grepl("^Abundance", cn))) {
      return("proteome.discoverer")
    }

    stop("Impossible to identify automatically the tool that generated this file.\n",
         "Set the 'source' argument explicitly, or use `source = \"generic\"` together with ",
         "'id.column' and 'sample.columns'.\nFirst columns found: ",
         paste(utils::head(cn, 10), collapse = ", "), call. = FALSE)
  }




#' @title .filter.long.report
#' @description Internal helper applying q-value and contaminant filters to a long report.
#' @param tb A data.frame.
#' @param id.col String, protein ID column.
#' @param precursor.q.col String (or NA), column of the precursor q-values.
#' @param protein.q.col String (or NA), column of the protein-group q-values.
#' @param q.value Numeric or NULL.
#' @param pg.q.value Numeric or NULL.
#' @param remove.contaminants Logical.
#' @param contaminant.pattern String or NULL.
#' @param verbose Logical.
#' @return A data.frame.
#' @keywords internal

.filter.long.report =
  function(tb,
           id.col,
           precursor.q.col = NA,
           protein.q.col = NA,
           q.value = 0.01,
           pg.q.value = 0.01,
           remove.contaminants = TRUE,
           contaminant.pattern = NULL,
           verbose = TRUE) {

    n0 = nrow(tb)

    if (!is.null(q.value) && !is.na(precursor.q.col) && precursor.q.col %in% colnames(tb)) {
      qv = suppressWarnings(as.numeric(tb[[precursor.q.col]]))
      tb = tb[!is.na(qv) & qv <= q.value, , drop = FALSE]
    }

    if (!is.null(pg.q.value) && !is.na(protein.q.col) && protein.q.col %in% colnames(tb)) {
      qv = suppressWarnings(as.numeric(tb[[protein.q.col]]))
      tb = tb[!is.na(qv) & qv <= pg.q.value, , drop = FALSE]
    }

    if (isTRUE(remove.contaminants)) {
      pattern = if (is.null(contaminant.pattern)) {
        "(^|;)(CON__|REV__|rev_|contam_|Cont_)|DECOY"
      } else {
        contaminant.pattern
      }
      tb = tb[!grepl(pattern, as.character(tb[[id.col]])), , drop = FALSE]
    }

    if (isTRUE(verbose) && nrow(tb) != n0) {
      message("Filtering removed ", format(n0 - nrow(tb), big.mark = ","), " of ",
              format(n0, big.mark = ","), " rows.")
    }

    if (nrow(tb) == 0) {
      stop("All rows were removed by the q-value/contaminant filters. Relax the thresholds.",
           call. = FALSE)
    }

    return(tb)
  }




#' @title .read.diann
#' @description Internal reader for the DIA-NN main report (long format, .tsv or .parquet).
#' @param tb A data.frame.
#' @param quantity String or \code{"auto"}.
#' @param protein.id String or \code{"auto"}.
#' @param sample.id String or \code{"auto"}.
#' @param summarization One among \code{"maxlfq"}, \code{"sum"}, \code{"median"}, \code{"none"}.
#' @param q.value Numeric or NULL.
#' @param pg.q.value Numeric or NULL.
#' @param remove.contaminants Logical.
#' @param contaminant.pattern String or NULL.
#' @param install.missing One among \code{"ask"}, \code{"always"}, \code{"never"}.
#' @param verbose Logical.
#' @return A list.
#' @keywords internal

.read.diann =
  function(tb, quantity, protein.id, sample.id, summarization,
           q.value, pg.q.value, remove.contaminants, contaminant.pattern,
           install.missing, verbose) {

    cn = colnames(tb)

    id.col  = .pick.column(protein.id, c("Protein.Group", "Protein.Ids", "Genes"), cn, "protein.id")
    run.col = .pick.column(sample.id, c("Run", "File.Name"), cn, "sample.id")

    if (is.na(id.col) || is.na(run.col)) {
      stop("Columns 'Protein.Group'/'Run' not found: this does not look like a DIA-NN main report.",
           call. = FALSE)
    }

    tb = .filter.long.report(tb = tb, id.col = id.col,
                             precursor.q.col = .pick.column("auto", c("Q.Value", "Global.Q.Value"), cn),
                             protein.q.col = .pick.column("auto", c("PG.Q.Value", "Global.PG.Q.Value", "Lib.PG.Q.Value"), cn),
                             q.value = q.value, pg.q.value = pg.q.value,
                             remove.contaminants = remove.contaminants,
                             contaminant.pattern = contaminant.pattern, verbose = verbose)

    if (summarization == "maxlfq") {
      quant.col = .pick.column(quantity, c("Precursor.Normalised", "Precursor.Quantity", "Ms1.Area"), cn, "quantity")
      mat = .maxlfq.rollup(tb = tb, primary.id = id.col,
                           secondary.id = c("Precursor.Id", "Modified.Sequence", "Precursor.Charge"),
                           sample.col = run.col, quantity.col = quant.col,
                           install.missing = install.missing, verbose = verbose)
      log.base = 2
      norm.method = paste0("DIA-NN + MaxLFQ (iq) on '", quant.col, "'")
      data.type = if (grepl("Normalised", quant.col)) "normalized" else "raw"

    } else if (summarization == "none") {
      quant.col = .pick.column(quantity, c("PG.MaxLFQ", "PG.Quantity", "Genes.MaxLFQ", "Genes.Quantity"), cn, "quantity")
      if (is.na(quant.col)) {
        stop("No protein-level quantity column found (e.g. 'PG.MaxLFQ'). ",
             "Use `summarization = \"maxlfq\"` to compute it from the precursors.", call. = FALSE)
      }
      ## protein-level values are constant within Protein.Group x Run: deduplicate first
      mat = .long.to.matrix(df = unique(tb[, c(id.col, run.col, quant.col)]),
                            id.col = id.col, sample.col = run.col, quantity.col = quant.col,
                            FUN = function(x) {mean(x, na.rm = TRUE)})
      log.base = 1
      norm.method = paste0("DIA-NN '", quant.col, "'")
      data.type = "normalized"

    } else {
      quant.col = .pick.column(quantity, c("Precursor.Normalised", "Precursor.Quantity"), cn, "quantity")
      FUN = if (summarization == "sum") {
        function(x) {sum(x, na.rm = TRUE)}
      } else {
        function(x) {stats::median(x, na.rm = TRUE)}
      }
      mat = .long.to.matrix(tb, id.col, run.col, quant.col, FUN = FUN)
      log.base = 1
      norm.method = paste0("DIA-NN '", quant.col, "' (", summarization, " roll-up)")
      data.type = if (grepl("Normalised", quant.col)) "normalized" else "raw"
    }

    return(list(counts = mat,
                log.base = log.base,
                data.type = data.type,
                normalization.method = norm.method,
                quantity = quant.col,
                feature = if (grepl("Genes", id.col)) "genes" else "protein groups"))
  }




#' @title .read.diann.matrix
#' @description Internal reader for the DIA-NN \code{*_matrix.tsv} files
#'   (pg_matrix, gg_matrix, pr_matrix, unique_genes_matrix).
#' @param tb A data.frame.
#' @param protein.id String or \code{"auto"}.
#' @param remove.contaminants Logical.
#' @param contaminant.pattern String or NULL.
#' @return A list.
#' @keywords internal

.read.diann.matrix =
  function(tb, protein.id, remove.contaminants, contaminant.pattern) {

    cn = colnames(tb)

    annotation.cols = c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                        "First.Protein.Description", "Precursor.Id", "Modified.Sequence",
                        "Stripped.Sequence", "Precursor.Charge", "Proteotypic")

    id.col = .pick.column(protein.id, c("Protein.Group", "Genes", "Precursor.Id", "Protein.Ids"), cn, "protein.id")
    if (is.na(id.col)) {
      stop("No identifier column found in this DIA-NN matrix file.", call. = FALSE)
    }

    quant.cols = setdiff(cn, annotation.cols)
    if (length(quant.cols) == 0) {
      stop("No run columns detected in the matrix file.", call. = FALSE)
    }
    names(quant.cols) = quant.cols

    if (isTRUE(remove.contaminants)) {
      pattern = if (is.null(contaminant.pattern)) "(^|;)(CON__|contam_|Cont_)" else contaminant.pattern
      tb = tb[!grepl(pattern, as.character(tb[[id.col]])), , drop = FALSE]
    }

    mat = .wide.to.matrix(tb, id.col, quant.cols)

    feature = if (grepl("Genes", id.col)) {
      "genes"
    } else if (grepl("Precursor", id.col)) {
      "precursors"
    } else {
      "protein groups"
    }

    return(list(counts = mat,
                log.base = 1,
                data.type = "normalized",
                normalization.method = "DIA-NN (MaxLFQ/QuantUMS, normalised matrix)",
                quantity = "matrix values",
                feature = feature))
  }




#' @title .read.spectronaut
#' @description Internal reader for long-format Spectronaut reports.
#' @param tb A data.frame.
#' @param quantity String or \code{"auto"}.
#' @param protein.id String or \code{"auto"}.
#' @param sample.id String or \code{"auto"}.
#' @param summarization One among \code{"maxlfq"}, \code{"sum"}, \code{"median"}, \code{"none"}.
#' @param q.value Numeric or NULL.
#' @param pg.q.value Numeric or NULL.
#' @param remove.contaminants Logical.
#' @param contaminant.pattern String or NULL.
#' @param install.missing One among \code{"ask"}, \code{"always"}, \code{"never"}.
#' @param verbose Logical.
#' @return A list.
#' @keywords internal

.read.spectronaut =
  function(tb, quantity, protein.id, sample.id, summarization,
           q.value, pg.q.value, remove.contaminants, contaminant.pattern,
           install.missing, verbose) {

    cn = colnames(tb)

    id.col  = .pick.column(protein.id, c("PG.ProteinGroups", "PG.ProteinAccessions", "PG.Genes"), cn, "protein.id")
    run.col = .pick.column(sample.id, c("R.FileName", "R.Raw.File.Name", "R.Label", "R.Condition"), cn, "sample.id")

    if (is.na(id.col) || is.na(run.col)) {
      stop("Columns 'PG.ProteinGroups'/'R.FileName' not found: this does not look like a long ",
           "Spectronaut report.", call. = FALSE)
    }

    tb = .filter.long.report(tb = tb, id.col = id.col,
                             precursor.q.col = .pick.column("auto", c("EG.Qvalue", "EG.QValue"), cn),
                             protein.q.col = .pick.column("auto", c("PG.Qvalue", "PG.QValue"), cn),
                             q.value = q.value, pg.q.value = pg.q.value,
                             remove.contaminants = remove.contaminants,
                             contaminant.pattern = contaminant.pattern, verbose = verbose)

    if (summarization == "maxlfq") {
      quant.col = .pick.column(quantity,
                               c("F.PeakArea", "F.NormalizedPeakArea", "FG.MS2Quantity", "FG.Quantity"),
                               cn, "quantity")
      mat = .maxlfq.rollup(tb = tb, primary.id = id.col,
                           secondary.id = c("EG.ModifiedSequence", "EG.PrecursorId", "FG.Charge",
                                            "F.FrgIon", "F.Charge", "F.FrgLossType"),
                           sample.col = run.col, quantity.col = quant.col,
                           install.missing = install.missing, verbose = verbose)
      log.base = 2
      norm.method = paste0("Spectronaut + MaxLFQ (iq) on '", quant.col, "'")
      data.type = if (grepl("Normalized", quant.col)) "normalized" else "raw"

    } else if (summarization == "none") {
      quant.col = .pick.column(quantity, c("PG.Quantity", "PG.MS2Quantity"), cn, "quantity")
      if (is.na(quant.col)) {
        stop("No protein-level quantity column found (e.g. 'PG.Quantity').", call. = FALSE)
      }
      mat = .long.to.matrix(df = unique(tb[, c(id.col, run.col, quant.col)]),
                            id.col = id.col, sample.col = run.col, quantity.col = quant.col,
                            FUN = function(x) {mean(x, na.rm = TRUE)})
      log.base = 1
      norm.method = paste0("Spectronaut '", quant.col, "'")
      data.type = "normalized"

    } else {
      quant.col = .pick.column(quantity, c("FG.Quantity", "FG.MS2Quantity", "F.PeakArea"), cn, "quantity")
      FUN = if (summarization == "sum") {
        function(x) {sum(x, na.rm = TRUE)}
      } else {
        function(x) {stats::median(x, na.rm = TRUE)}
      }
      mat = .long.to.matrix(tb, id.col, run.col, quant.col, FUN = FUN)
      log.base = 1
      norm.method = paste0("Spectronaut '", quant.col, "' (", summarization, " roll-up)")
      data.type = "raw"
    }

    return(list(counts = mat,
                log.base = log.base,
                data.type = data.type,
                normalization.method = norm.method,
                quantity = quant.col,
                feature = "protein groups"))
  }




#' @title .read.spectronaut.pivot
#' @description Internal reader for Spectronaut pivot exports, i.e. one column per run named
#'   \code{"[1] sample.raw.PG.Quantity"}.
#' @param tb A data.frame.
#' @param quantity String or \code{"auto"}.
#' @param protein.id String or \code{"auto"}.
#' @param remove.contaminants Logical.
#' @param contaminant.pattern String or NULL.
#' @return A list.
#' @keywords internal

.read.spectronaut.pivot =
  function(tb, quantity, protein.id, remove.contaminants, contaminant.pattern) {

    cn = colnames(tb)

    id.col = .pick.column(protein.id, c("PG.ProteinGroups", "PG.ProteinAccessions", "PG.Genes"), cn, "protein.id")
    if (is.na(id.col)) {
      stop("No 'PG.ProteinGroups'-like column found in this pivot report.", call. = FALSE)
    }

    quant.suffix = if (identical(tolower(quantity), "auto")) "PG.Quantity" else quantity
    suffix.regex = paste0("\\.", gsub("\\.", "\\\\.", quant.suffix), "$")

    hits = grep(suffix.regex, cn, value = TRUE)
    if (length(hits) == 0) {
      stop("No column ending in '", quant.suffix, "' was found.\nAvailable columns: ",
           paste(utils::head(cn, 20), collapse = ", "), call. = FALSE)
    }

    quant.cols = hits
    names(quant.cols) = sub(suffix.regex, "", hits)

    if (isTRUE(remove.contaminants)) {
      pattern = if (is.null(contaminant.pattern)) "(^|;)(CON__|contam_|Cont_)" else contaminant.pattern
      tb = tb[!grepl(pattern, as.character(tb[[id.col]])), , drop = FALSE]
    }

    mat = .wide.to.matrix(tb, id.col, quant.cols)

    return(list(counts = mat,
                log.base = 1,
                data.type = "normalized",
                normalization.method = paste0("Spectronaut '", quant.suffix, "' (pivot export)"),
                quantity = quant.suffix,
                feature = "protein groups"))
  }




#' @title .read.fragpipe
#' @description Internal reader for FragPipe \code{combined_protein.tsv} / \code{combined_peptide.tsv}.
#' @param tb A data.frame.
#' @param quantity String or \code{"auto"}.
#' @param protein.id String or \code{"auto"}.
#' @param remove.contaminants Logical.
#' @param contaminant.pattern String or NULL.
#' @return A list.
#' @keywords internal

.read.fragpipe =
  function(tb, quantity, protein.id, remove.contaminants, contaminant.pattern) {

    cn = colnames(tb)

    id.col = .pick.column(protein.id, c("Protein ID", "Protein", "Gene", "Peptide Sequence"), cn, "protein.id")
    if (is.na(id.col)) {
      stop("No 'Protein ID'-like column found in this FragPipe report.", call. = FALSE)
    }

    ## per-sample suffixes used by FragPipe/IonQuant
    all.suffixes = c("MaxLFQ Unique Intensity", "MaxLFQ Total Intensity", "MaxLFQ Razor Intensity",
                     "MaxLFQ Intensity", "Unique Spectral Count", "Total Spectral Count",
                     "Unique Intensity", "Total Intensity", "Razor Intensity",
                     "Spectral Count", "Intensity")

    quant.suffix = if (identical(tolower(quantity), "auto")) "MaxLFQ Intensity" else quantity
    all.suffixes = unique(c(quant.suffix, all.suffixes))

    ## a shorter suffix must not swallow the columns of a longer one ("Intensity" vs "MaxLFQ Intensity")
    longer = all.suffixes[nchar(all.suffixes) > nchar(quant.suffix)]

    hits = grep(paste0(" ", quant.suffix, "$"), cn, value = TRUE, fixed = FALSE)
    for (s in longer) {
      hits = setdiff(hits, grep(paste0(" ", s, "$"), cn, value = TRUE, fixed = FALSE))
    }

    if (length(hits) == 0) {
      available = all.suffixes[vapply(all.suffixes,
                                      function(s) {any(grepl(paste0(" ", s, "$"), cn))},
                                      FUN.VALUE = logical(1))]
      stop("No column ending in '", quant.suffix, "' was found.\nQuantity suffixes present in this file: ",
           paste(available, collapse = ", "), call. = FALSE)
    }

    quant.cols = hits
    names(quant.cols) = trimws(sub(paste0(" ", quant.suffix, "$"), "", hits))

    if (isTRUE(remove.contaminants)) {
      pattern = if (is.null(contaminant.pattern)) "^(contam_|rev_|sp\\|contam)" else contaminant.pattern
      check.col = intersect(c("Protein", "Protein ID"), cn)[1]
      if (!is.na(check.col)) {
        tb = tb[!grepl(pattern, as.character(tb[[check.col]]), ignore.case = TRUE), , drop = FALSE]
      }
    }

    mat = .wide.to.matrix(tb, id.col, quant.cols)

    return(list(counts = mat,
                log.base = 1,
                data.type = if (grepl("MaxLFQ", quant.suffix)) "normalized" else "raw",
                normalization.method = paste0("FragPipe/IonQuant '", quant.suffix, "'"),
                quantity = quant.suffix,
                feature = if (grepl("Peptide", id.col)) "peptides" else "proteins"))
  }




#' @title .read.maxquant
#' @description Internal reader for MaxQuant \code{proteinGroups.txt}.
#' @param tb A data.frame.
#' @param quantity String or \code{"auto"}.
#' @param protein.id String or \code{"auto"}.
#' @param remove.contaminants Logical.
#' @param contaminant.pattern String or NULL.
#' @param verbose Logical.
#' @return A list.
#' @keywords internal

.read.maxquant =
  function(tb, quantity, protein.id, remove.contaminants, contaminant.pattern, verbose) {

    cn = colnames(tb)

    id.col = .pick.column(protein.id, c("Majority protein IDs", "Protein IDs", "Gene names"), cn, "protein.id")
    if (is.na(id.col)) {
      stop("No 'Majority protein IDs'-like column found in this MaxQuant report.", call. = FALSE)
    }

    ## MaxQuant flags contaminants/decoys in dedicated columns
    if (isTRUE(remove.contaminants)) {
      n0 = nrow(tb)
      for (flag in c("Reverse", "Potential contaminant", "Contaminant", "Only identified by site")) {
        if (flag %in% cn) {
          tb = tb[!(as.character(tb[[flag]]) %in% c("+", "TRUE", "1")), , drop = FALSE]
        }
      }
      pattern = if (is.null(contaminant.pattern)) "(^|;)(CON__|REV__)" else contaminant.pattern
      tb = tb[!grepl(pattern, as.character(tb[[id.col]])), , drop = FALSE]

      if (isTRUE(verbose) && nrow(tb) != n0) {
        message("Removed ", n0 - nrow(tb), " contaminant/reverse/site-only entries.")
      }
    }

    quant.prefix = if (identical(tolower(quantity), "auto")) "LFQ intensity" else quantity

    hits = grep(paste0("^", quant.prefix, " "), cn, value = TRUE)
    if (length(hits) == 0) {
      stop("No column starting with '", quant.prefix, " ' was found.\n",
           "Try `quantity = \"Intensity\"` or `quantity = \"iBAQ\"`.", call. = FALSE)
    }

    quant.cols = hits
    names(quant.cols) = trimws(sub(paste0("^", quant.prefix, " "), "", hits))

    mat = .wide.to.matrix(tb, id.col, quant.cols)

    return(list(counts = mat,
                log.base = 1,
                data.type = if (grepl("LFQ", quant.prefix)) "normalized" else "raw",
                normalization.method = paste0("MaxQuant '", quant.prefix, "'"),
                quantity = quant.prefix,
                feature = "protein groups"))
  }




#' @title .read.proteome.discoverer
#' @description Internal reader for Proteome Discoverer protein-level exports
#'   (R-friendly headers recommended).
#' @param tb A data.frame.
#' @param quantity String or \code{"auto"}.
#' @param protein.id String or \code{"auto"}.
#' @param remove.contaminants Logical.
#' @param contaminant.pattern String or NULL.
#' @return A list.
#' @keywords internal

.read.proteome.discoverer =
  function(tb, quantity, protein.id, remove.contaminants, contaminant.pattern) {

    cn = colnames(tb)

    id.col = .pick.column(protein.id,
                          c("Accession", "Master.Protein.Accessions", "Protein.Accessions", "Protein Accessions"),
                          cn, "protein.id")
    if (is.na(id.col)) {
      stop("No 'Accession'-like column found in this Proteome Discoverer export.", call. = FALSE)
    }

    quant.prefix = if (identical(tolower(quantity), "auto")) "Abundances (Normalized)" else quantity

    hits = grep(paste0("^", gsub("([().|\\[\\]])", "\\\\\\1", quant.prefix)), cn, value = TRUE)
    if (length(hits) == 0) {
      quant.prefix = "Abundance"
      hits = grep("^Abundances?[:.]", cn, value = TRUE)
    }
    if (length(hits) == 0) {
      stop("No abundance column was found. Set 'quantity' explicitly.\nAvailable columns: ",
           paste(utils::head(cn, 20), collapse = ", "), call. = FALSE)
    }

    quant.cols = hits
    ## "Abundances (Normalized): F1: Sample, WT" -> "Sample, WT"
    sample.names = trimws(sub("^.*?[:.]\\s*F[0-9]+[:.]?\\s*", "", hits))
    if (any(duplicated(sample.names)) || any(sample.names == "") || identical(sample.names, hits)) {
      sample.names = make.unique(gsub("[^A-Za-z0-9._-]+", ".", hits), sep = ".")
    }
    names(quant.cols) = sample.names

    if (isTRUE(remove.contaminants)) {
      pattern = if (is.null(contaminant.pattern)) "(^|;)(CON__|contam)" else contaminant.pattern
      tb = tb[!grepl(pattern, as.character(tb[[id.col]]), ignore.case = TRUE), , drop = FALSE]
      if ("Contaminant" %in% cn) {
        tb = tb[!(as.character(tb[["Contaminant"]]) %in% c("TRUE", "True", "+", "1")), , drop = FALSE]
      }
    }

    mat = .wide.to.matrix(tb, id.col, quant.cols)

    return(list(counts = mat,
                log.base = 1,
                data.type = if (grepl("Normalized", quant.prefix)) "normalized" else "raw",
                normalization.method = paste0("Proteome Discoverer '", quant.prefix, "'"),
                quantity = quant.prefix,
                feature = "proteins"))
  }




#' @title .read.generic
#' @description Internal reader for any wide table, when the user provides the ID and quantity columns.
#' @param tb A data.frame.
#' @param id.column String, protein ID column.
#' @param sample.columns Character vector (optionally named) of the intensity columns.
#' @return A list.
#' @keywords internal

.read.generic =
  function(tb, id.column, sample.columns) {

    if (is.null(id.column) || is.null(sample.columns)) {
      stop("With `source = \"generic\"` both 'id.column' and 'sample.columns' must be provided.",
           call. = FALSE)
    }

    missing.cols = setdiff(c(id.column, unname(sample.columns)), colnames(tb))
    if (length(missing.cols) > 0) {
      stop("Column(s) not found in the file: ", paste(missing.cols, collapse = ", "), call. = FALSE)
    }

    if (is.null(names(sample.columns))) {
      names(sample.columns) = sample.columns
    }

    mat = .wide.to.matrix(tb, id.column, sample.columns)

    return(list(counts = mat,
                log.base = 1,
                data.type = "raw",
                normalization.method = NA,
                quantity = "user-defined",
                feature = "features"))
  }
