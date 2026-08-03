#' @title import.msstats
#'
#' @description Converts the protein-level output of \code{MSstats::dataProcess()} (label-free) or
#'   \code{MSstatsTMT::proteinSummarization()} (TMT) into a \code{DEprot} object. Contrary to
#'   \link{import.external}, which parses the reports written on disk by the search engines, this
#'   function takes an R object that has already been summarized, normalized and log-transformed by
#'   MSstats. The object is then built by \link{load.counts2}.
#'
#' @param object The list returned by \code{dataProcess()} / \code{proteinSummarization()}, the
#'   \code{ProteinLevelData} (or legacy \code{RunlevelData}) data.frame alone, or the path to an
#'   \code{.rds} file containing either of them.
#' @param metadata A data.frame (or path to a table) containing at least the column indicated by
#'   \code{column.id}. If \code{NULL} (default) the metadata are reconstructed from the annotation
#'   columns carried by the MSstats object itself (conditions, replicates, mixtures, channels).
#' @param type String indicating the flavour of the input: \code{"auto"} (default), \code{"lfq"}
#'   (MSstats) or \code{"tmt"} (MSstatsTMT).
#' @param protein.id String indicating the column to use as protein identifier. Default \code{"auto"}.
#' @param sample.id Character vector indicating the column(s) whose combination identifies a sample.
#'   Default \code{"auto"} uses \code{originalRUN} for label-free data and \code{Run} + \code{Channel}
#'   for TMT, which is the only combination guaranteed to be unique across mixtures and technical
#'   replicates.
#' @param sample.sep String used to paste multiple \code{sample.id} columns. Default \code{"_"}.
#' @param log.base Number indicating the base of the logarithm of the summarized abundances.
#'   Default \code{2}. Set it to \code{10} if \code{dataProcess()} was run with \code{logTrans = 10}.
#' @param annotation.columns Character vector of the columns to carry into the automatically
#'   generated metadata. Default \code{"auto"} keeps the MSstats annotation columns that are constant
#'   within each sample.
#' @param column.id String indicating the metadata column matching the sample names.
#'   Default \code{"column.id"}.
#' @param data.type String, one among \code{"auto"} (default), \code{"raw"}, \code{"normalized"},
#'   \code{"randomized"}, \code{"imputed"}. See 'Details' for how \code{"auto"} decides.
#' @param clean.sample.names Logical, whether to strip paths and raw-file extensions from the run
#'   names. Default \code{TRUE}.
#' @param subset.to.metadata Logical, whether to keep only the samples listed in the metadata.
#'   Default \code{TRUE}.
#' @param verbose Logical, whether to print progress messages. Default \code{TRUE}.
#'
#' @details
#' No package is required to read these objects: they are ordinary data.frames, so neither
#' \code{MSstats} nor \code{MSstatsTMT} needs to be installed for the import itself.
#'
#' MSstats abundances are already log-transformed, therefore the zero-to-\code{NA} conversion applied
#' to linear intensities is skipped here: a summarized abundance of 0 is a legitimate value.
#'
#' With \code{data.type = "auto"} the label-free data are loaded as \code{"imputed"} when the column
#' \code{NumImputedFeature} is present and contains non-zero values, since that column only appears
#' when \code{dataProcess()} was run with \code{MBimpute = TRUE}; otherwise they are loaded as
#' \code{"normalized"}. TMT data are always loaded as \code{"normalized"}, because
#' \code{proteinSummarization()} applies \code{global_norm} and \code{reference_norm} by default and
#' its protein-level table carries no indicator of the accelerated-failure-time imputation. Set
#' \code{data.type = "imputed"} explicitly if \code{MBimpute = TRUE} was used.
#'
#' @section A note on TMT data:
#' \code{DEprot} is designed for label-free quantitation. TMT data can be loaded and analysed, but
#' several assumptions of the package do not hold: the missingness is structured by mixture rather
#' than by sample, the fold changes are compressed by co-isolation, and the data have already been
#' normalized within and between plexes by MSstatsTMT. Re-normalizing with \link{normalize.counts}
#' or re-imputing with \link{impute.counts} is therefore discouraged. The mixture effect can instead
#' be handled with \link{harmonize.batches}, using the \code{Mixture} column of the generated metadata
#' as batch variable.
#'
#' @return A \code{DEprot} object (S4 vector).
#'
#' @import dplyr
#' @import methods
#' @importFrom data.table fread
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' \dontrun{
#' # label-free (MSstats)
#' summarized <- MSstats::dataProcess(raw)
#' dpo <- import.msstats(object = summarized)
#'
#' # TMT (MSstatsTMT)
#' summarized.tmt <- MSstatsTMT::proteinSummarization(input.pd)
#' dpo <- import.msstats(object = summarized.tmt,
#'                       type = "tmt",
#'                       data.type = "imputed")
#'
#' # correcting the mixture (plex) effect of a TMT experiment
#' dpo <- harmonize.batches(DEprot.object = dpo, batch.column = "Mixture")
#' }
#'
#' @seealso \link{import.external}, \link{load.counts2}
#'
#' @export import.msstats

import.msstats =
  function(object,
           metadata = NULL,
           type = "auto",
           protein.id = "auto",
           sample.id = "auto",
           sample.sep = "_",
           log.base = 2,
           annotation.columns = "auto",
           column.id = "column.id",
           data.type = "auto",
           clean.sample.names = TRUE,
           subset.to.metadata = TRUE,
           verbose = TRUE) {

    ### ------------------------------------------------------------------- checks
    if (missing(object)) {
      stop("Provide 'object': the output of MSstats::dataProcess() or ",
           "MSstatsTMT::proteinSummarization().", call. = FALSE)
    }

    type = match.arg(tolower(as.character(type)[1]), choices = c("auto", "lfq", "tmt"))


    ### ------------------------------------------------------------------- extract the protein table
    if (is.character(object)) {
      if (!file.exists(object)) {
        stop("File not found: ", object, call. = FALSE)
      }
      object = readRDS(object)
    }

    tb = .extract.msstats.table(object)

    if (type == "auto") {
      type = if (all(c("Channel", "Mixture") %in% colnames(tb))) "tmt" else "lfq"
      if (isTRUE(verbose)) {
        message("Detected input: MSstats", ifelse(type == "tmt", "TMT (isobaric labelling)",
                                                  " (label-free)"), ".")
      }
    }


    ### ------------------------------------------------------------------- dispatch
    imported =
      if (type == "tmt") {
        .read.msstats.tmt(tb = tb, protein.id = protein.id, sample.id = sample.id,
                          sample.sep = sample.sep, annotation.columns = annotation.columns,
                          column.id = column.id, log.base = log.base)
      } else {
        .read.msstats.lfq(tb = tb, protein.id = protein.id, sample.id = sample.id,
                          sample.sep = sample.sep, annotation.columns = annotation.columns,
                          column.id = column.id, log.base = log.base)
      }

    if (type == "tmt" && isTRUE(verbose)) {
      message("Note: these are isobaric (TMT) data, already normalized within and between mixtures ",
              "by MSstatsTMT. Re-normalization and re-imputation are discouraged; the mixture effect ",
              "can be handled with `harmonize.batches()`.")
    }


    ### ------------------------------------------------------------------- build the object
    DEprot.object =
      .finalize.import(imported = imported,
                       metadata = metadata,
                       column.id = column.id,
                       data.type = data.type,
                       log.base = log.base,
                       clean.sample.names = clean.sample.names,
                       subset.to.metadata = subset.to.metadata,
                       verbose = verbose)

    return(DEprot.object)
  } # END function




# =====================================================================================
# Internal helpers
# =====================================================================================


#' @title .extract.msstats.table
#' @description Internal helper pulling the protein-level table out of the several shapes in which
#'   MSstats/MSstatsTMT/MSstatsPTM return their results.
#' @param object A list, a data.frame or a data.table.
#' @return A data.frame.
#' @keywords internal

.extract.msstats.table =
  function(object) {

    if (is.data.frame(object)) {
      return(as.data.frame(object))
    }

    if (is.list(object)) {
      ## MSstatsPTM nests the results under $PROTEIN / $PTM
      if ("PROTEIN" %in% names(object) && is.list(object$PROTEIN)) {
        object = object$PROTEIN
      }

      ## current MSstats/MSstatsTMT name, then the legacy MSstats one
      for (slot in c("ProteinLevelData", "RunlevelData")) {
        if (slot %in% names(object)) {
          return(as.data.frame(object[[slot]]))
        }
      }

      stop("The provided list does not contain a 'ProteinLevelData' (or 'RunlevelData') element.\n",
           "Available elements: ", paste(names(object), collapse = ", "), call. = FALSE)
    }

    stop("'object' must be the list returned by dataProcess()/proteinSummarization(), the ",
         "corresponding ProteinLevelData data.frame, or the path to an .rds file containing them.",
         call. = FALSE)
  }




#' @title .build.msstats.samples
#' @description Internal helper building the sample identifier and the annotation table shared by the
#'   label-free and the TMT readers.
#' @param tb A data.frame.
#' @param sample.id Character vector or \code{"auto"}.
#' @param default.id Character vector of the columns to use when \code{sample.id = "auto"}.
#' @param sample.sep String separating the pasted components.
#' @param annotation.columns Character vector or \code{"auto"}.
#' @param default.annotation Character vector of the candidate annotation columns.
#' @param column.id String, name of the sample column in the output annotation table.
#' @return A list with the modified table (column \code{.sample.id}) and the annotation data.frame.
#' @keywords internal

.build.msstats.samples =
  function(tb,
           sample.id,
           default.id,
           sample.sep,
           annotation.columns,
           default.annotation,
           column.id) {

    cn = colnames(tb)

    ## which columns define a sample
    if (identical(tolower(as.character(sample.id)[1]), "auto")) {
      keys = intersect(default.id, cn)
      if (length(keys) == 0) {
        stop("None of the expected run/channel columns (", paste(default.id, collapse = ", "),
             ") was found in the MSstats table.", call. = FALSE)
      }
    } else {
      keys = as.character(sample.id)
      missing.keys = setdiff(keys, cn)
      if (length(missing.keys) > 0) {
        stop("The 'sample.id' column(s) ", paste(missing.keys, collapse = ", "),
             " are not present in the MSstats table.\nAvailable columns: ",
             paste(cn, collapse = ", "), call. = FALSE)
      }
    }

    ## clean the run component before pasting, so that the final names stay readable
    parts = lapply(keys, function(k) {
      x = as.character(tb[[k]])
      if (grepl("run", k, ignore.case = TRUE)) {x = .clean.run.names(x)}
      return(x)
    })

    tb$.sample.id = do.call(paste, c(parts, list(sep = sample.sep)))

    ## annotation columns that are constant within a sample
    if (identical(tolower(as.character(annotation.columns)[1]), "auto")) {
      candidates = setdiff(intersect(default.annotation, cn), keys)
    } else {
      candidates = intersect(as.character(annotation.columns), cn)
    }

    constant = candidates[vapply(X = candidates,
                                 FUN = function(k) {
                                   all(tapply(as.character(tb[[k]]), tb$.sample.id,
                                              function(v) {length(unique(v)) == 1}))
                                 },
                                 FUN.VALUE = logical(1))]

    annotation = unique(tb[, c(".sample.id", keys, constant), drop = FALSE])
    colnames(annotation)[1] = column.id
    rownames(annotation) = NULL

    if (any(duplicated(annotation[[column.id]]))) {
      ## should not happen, but never return an ambiguous metadata table
      annotation = annotation[!duplicated(annotation[[column.id]]), , drop = FALSE]
    }

    return(list(tb = tb, annotation = annotation))
  }




#' @title .read.msstats.lfq
#' @description Internal reader for the protein-level output of \code{MSstats::dataProcess()}.
#' @param tb A data.frame.
#' @param protein.id String or \code{"auto"}.
#' @param sample.id Character vector or \code{"auto"}.
#' @param sample.sep String.
#' @param annotation.columns Character vector or \code{"auto"}.
#' @param column.id String.
#' @param log.base Numeric.
#' @return A list.
#' @keywords internal

.read.msstats.lfq =
  function(tb, protein.id, sample.id, sample.sep, annotation.columns, column.id, log.base) {

    cn = colnames(tb)

    id.col = .pick.column(protein.id, c("Protein", "PROTEIN", "ProteinName"), cn, "protein.id")
    quant.col = .pick.column("auto", c("LogIntensities", "Abundance", "ABUNDANCE"), cn)

    if (is.na(id.col) || is.na(quant.col)) {
      stop("Columns 'Protein'/'LogIntensities' not found: this does not look like the ",
           "ProteinLevelData of MSstats::dataProcess().\nAvailable columns: ",
           paste(cn, collapse = ", "), call. = FALSE)
    }

    built = .build.msstats.samples(tb = tb, sample.id = sample.id,
                                   default.id = c("originalRUN", "RUN", "Run"),
                                   sample.sep = sample.sep,
                                   annotation.columns = annotation.columns,
                                   default.annotation = c("GROUP", "GROUP_ORIGINAL", "SUBJECT",
                                                          "SUBJECT_ORIGINAL", "Condition",
                                                          "BioReplicate", "originalRUN", "RUN"),
                                   column.id = column.id)

    mat = .long.to.matrix(df = built$tb, id.col = id.col, sample.col = ".sample.id",
                          quantity.col = quant.col,
                          FUN = function(x) {mean(x, na.rm = TRUE)},
                          zero.to.na = FALSE)

    ## NumImputedFeature is written only when dataProcess() was run with MBimpute = TRUE
    imputed = "NumImputedFeature" %in% cn &&
      any(suppressWarnings(as.numeric(tb[["NumImputedFeature"]])) > 0, na.rm = TRUE)

    return(list(counts = mat,
                log.base = log.base,
                data.type = if (imputed) "imputed" else "normalized",
                normalization.method = "MSstats dataProcess() summarization",
                imputation.method = if (imputed) "MSstats accelerated failure time model (MBimpute)" else NA,
                quantity = quant.col,
                feature = "proteins",
                annotation = built$annotation))
  }




#' @title .read.msstats.tmt
#' @description Internal reader for the protein-level output of
#'   \code{MSstatsTMT::proteinSummarization()}.
#' @param tb A data.frame.
#' @param protein.id String or \code{"auto"}.
#' @param sample.id Character vector or \code{"auto"}.
#' @param sample.sep String.
#' @param annotation.columns Character vector or \code{"auto"}.
#' @param column.id String.
#' @param log.base Numeric.
#' @return A list.
#' @keywords internal

.read.msstats.tmt =
  function(tb, protein.id, sample.id, sample.sep, annotation.columns, column.id, log.base) {

    cn = colnames(tb)

    id.col = .pick.column(protein.id, c("Protein", "ProteinName", "PROTEIN"), cn, "protein.id")
    quant.col = .pick.column("auto", c("Abundance", "LogIntensities"), cn)

    if (is.na(id.col) || is.na(quant.col)) {
      stop("Columns 'Protein'/'Abundance' not found: this does not look like the ProteinLevelData ",
           "of MSstatsTMT::proteinSummarization().\nAvailable columns: ",
           paste(cn, collapse = ", "), call. = FALSE)
    }

    if (!("Channel" %in% cn)) {
      stop("No 'Channel' column found: for isobaric data the sample is defined by the combination ",
           "of MS run and labelling channel.", call. = FALSE)
    }

    built = .build.msstats.samples(tb = tb, sample.id = sample.id,
                                   default.id = c("Run", "Channel"),
                                   sample.sep = sample.sep,
                                   annotation.columns = annotation.columns,
                                   default.annotation = c("Mixture", "TechRepMixture", "Condition",
                                                          "BioReplicate", "Channel", "Run"),
                                   column.id = column.id)

    mat = .long.to.matrix(df = built$tb, id.col = id.col, sample.col = ".sample.id",
                          quantity.col = quant.col,
                          FUN = function(x) {mean(x, na.rm = TRUE)},
                          zero.to.na = FALSE)

    return(list(counts = mat,
                log.base = log.base,
                data.type = "normalized",
                normalization.method = "MSstatsTMT proteinSummarization() (global + reference channel)",
                imputation.method = NA,
                quantity = quant.col,
                feature = "proteins",
                annotation = built$annotation))
  }
