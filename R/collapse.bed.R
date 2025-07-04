#' @title Merger of overlapping peaks in a provided .bed file.
#'
#' @description Merge overlapping peaks in a provided .bed file.
#'
#' @param bed Two options are possible: \cr - String with the path to a .bed file; \cr - data.frame corresponding to a bed file format (only the first 6 columns, BED6, will be kept).
#' @param maximal.distance Maximal distance between regions allowed for regions to be merged. By default \code{0}.
#' @param keep.strandness Logic value to indicate whether to force to only merge regions that are in the same strand. By default \code{FALSE}, disabled. Subordinated to not \code{NULL} value for 'only.one.strand' option.
#' @param only.one.strand Atomic string to indicate whether to force merge for one specific strand only. It must be indicated the wished strand (e.g., '+', '-', '.'). Regions in the other strand/s will be kept without any modification. By default \code{NULL}.
#' @param score.operation Applicable only if the regions contain scores. Atomic string to indicate the operation to apply to the scores of merged regions. Possible choices: 'mean', 'median', 'sum'. By default \code{"mean"}.
#' @param bed.header Logic value to define whether the .bed file contains an header or not. By default \code{FALSE}.
#' @param sep String containing the separator character for a .bed file. By default \code{"\t"}.
#' @param return.bed Logic value to define if to return the bed as a data.frame. By default \code{TRUE}. Only unique rows are kept.
#' @param export.file.name Optional: string to define the path to the file to be exported, if required. By default \code{NULL}, not exported.
#' @param export.header Logic value to define whether the header should be exported in the sorted bed file. By default \code{FALSE}.
#' @param verbose Logic value to indicate whether messages should be printed or not. By default \code{TRUE}.
#'
#' @return If required, returns a data.frame corresponding to the collapsed .bed file.
#'
#' @details The function pre-sorts the bed and keeps only unique rows and only up to 6 columns (chr, start, end, name, score, strand). \cr The names of the regions (if available) of merged regions corresponds to the concatenation of all original region's name. \cr To get more information about the bed file format see the following page: \cr \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}.
#'
#' @export collapse.bed
#'
# @import dplyr


collapse.bed = function(bed,
                        maximal.distance = 0,
                        keep.strandness = FALSE,
                        only.one.strand = NULL,
                        score.operation = "mean",
                        bed.header = FALSE,
                        sep = "\t",
                        return.bed = TRUE,
                        export.file.name = NULL,
                        export.header = FALSE,
                        verbose = TRUE) {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#


  ### -------------------------------------------------------------------------------- ###
  ###                                Libraries & Controls                              ###
  ### -------------------------------------------------------------------------------- ###
  # Load libraries
  require(dplyr)


  # Check scoreoperation
  if (!(score.operation %in% c("mean", "median", "sum"))) {
    return(warning("The 'score.operation' parameter must be one among: 'mean', 'median', 'sum'. Default value = 'mean'."))
  }



  # Read and sort bed, keeping only first 6 columns
  sorted.bed = Rseb::sort.bed(bed = bed,
                              bed.header = bed.header,
                              sep = sep,
                              verbose = F)


  # Select the first 6 columns or add dummy columns if missing column 4|5|6
  if (ncol(sorted.bed) >= 6) {
    sorted.bed = sorted.bed[,1:6]
    final.ncol = 6
  } else if (ncol(sorted.bed) == 3) {
    sorted.bed = dplyr::mutate(sorted.bed,
                               name = ".",
                               score = 0,
                               strand = "+")
    final.ncol = 3
  } else if (ncol(sorted.bed) == 4) {
    sorted.bed = dplyr::mutate(sorted.bed,
                               score = 0,
                               strand = "+")
    final.ncol = 4
  } else if (ncol(sorted.bed) == 5) {
    sorted.bed = dplyr::mutate(sorted.bed, strand = "+")
    final.ncol = 5
  }

  names(sorted.bed) = c("chr", "start", "end", "name", "score", "strand")


  # in case of only one strand, check if the strand exists
  if (!is.null(only.one.strand)) {
    if (!(only.one.strand %in% sorted.bed$strand) & (verbose == T)) {
      only.one.strand = NULL
      message(paste0("The 'only.one.strand' indicate is '",
                     only.one.strand,
                     "', but no regions are displaying this strand.\n",
                     "For this reason, the 'only.one.strand' option has been coerced to 'NULL' value.\n",
                     "The collapsing will be performed ignoring the strandness."))
    }
  }


  # Check that END.position > START.position
  if (length((unique((sorted.bed$end - sorted.bed$start) < 0))) > 1) {
    return(warning("The input bed file contains regions in which the END boundary position is lower than the START one. Collapsing interrupted with no output."))
  }
  ### -------------------------------------------------------------------------------- ###



  ### -------------------------------------------------------------------------------- ###
  ###                                    Function                                      ###
  ### -------------------------------------------------------------------------------- ###
  # Add the maximal.distance to each side of each peak
  sorted.bed = dplyr::mutate(.data = sorted.bed,
                             start = start - maximal.distance,
                             end = end + maximal.distance)


  #F# define the function for the overlap's merge
  merge.peaks = function(sorted.bed) {
    merged.bed = data.frame()

    for (c in unique(sorted.bed$chr)) {
      current.chr.bed = dplyr::filter(.data = sorted.bed, chr == c)

      # Add check columns
      current.chr.bed = dplyr::mutate(.data = current.chr.bed,
                                      analyse = T,
                                      ID = 1:nrow(current.chr.bed))

      # Identify the overlaps
      for (i in 1:nrow(current.chr.bed)) {
        if (current.chr.bed$analyse[i] == T) {
          current_overlap = dplyr::filter(.data = current.chr.bed,
                                          (start == current.chr.bed$start[i]) | (start <= current.chr.bed$end[i]),
                                          analyse == T)

          if (nrow(current_overlap) == 1) {
            merged.bed = rbind(merged.bed, current_overlap[,1:6])
            current.chr.bed$analyse[i] = F

          } else {
            merged.region = data.frame(chr = c,
                                       start = min(current_overlap$start),
                                       end = max(current_overlap$end),
                                       name = paste(unique(current_overlap$name), collapse = "_"),
                                       score = ifelse(test = score.operation == "mean",
                                                      yes = mean(current_overlap$score),
                                                      no = ifelse(test = score.operation == "median",
                                                                  yes = median(current_overlap$score),
                                                                  no = sum(current_overlap$score))),
                                       strand = ifelse(test = length(unique(current_overlap$strand)) == 1,
                                                       yes = unique(current_overlap$strand),
                                                       no = "."))

            merged.bed = rbind(merged.bed, merged.region)

            current.chr.bed$analyse[current.chr.bed$ID %in% current_overlap$ID] = FALSE
          }
        }


      } #> END overlap for loop

    } #> END chr for loop

    return(merged.bed[,1:6])
  } #F# END merge.peaks




  collapsed.bed = data.frame()

  # Compute overlap merge for only one strand (if only.one.strand != NULL)
  if (!is.null(only.one.strand)) {
    strand.collapsed.bed = data.frame()
    continue = T
    to.collapse.bed = dplyr::filter(sorted.bed, strand == only.one.strand)

    # Repeat the overlapping until no anymore overlaps are found
    while (continue == T) {
      strand.collapsed.bed = merge.peaks(sorted.bed = to.collapse.bed)

      continue = ifelse(test = nrow(strand.collapsed.bed) == nrow(to.collapse.bed),
                        yes = F,
                        no = T)

      to.collapse.bed = strand.collapsed.bed
    }

    collapsed.bed = Rseb::sort.bed(rbind(strand.collapsed.bed[,1:6],
                                   dplyr::filter(sorted.bed, strand != only.one.strand),
                                   verbose = F))
  }


  # Compute overlap keeping strandness (if only.one.strand == NULL & keep.strandness == TRUE)
  if (is.null(only.one.strand) & keep.strandness == TRUE) {
    bed_list = list()

    for (s in 1:length(unique(sorted.bed$strand))) {
      strand.collapsed.bed = data.frame()
      continue = T
      to.collapse.bed = dplyr::filter(sorted.bed, strand == unique(sorted.bed$strand)[s])

      # Repeat the overlapping until no anymore overlaps are found
      while (continue == T) {
        strand.collapsed.bed = merge.peaks(sorted.bed = to.collapse.bed)

        continue = ifelse(test = nrow(strand.collapsed.bed) == nrow(to.collapse.bed),
                          yes = F,
                          no = T)

        to.collapse.bed = strand.collapsed.bed
      }

      bed_list[[s]] = strand.collapsed.bed
    }

    collapsed.bed = Rseb::sort.bed(purrr::reduce(.x = bed_list, .f = rbind), verbose = F)
  }



  # Compute overlap without strand resctriction
  if (is.null(only.one.strand) & keep.strandness == FALSE) {
    to.collapse.bed = sorted.bed
    continue = T

    while (continue == T) {
      collapsed.bed = merge.peaks(sorted.bed = to.collapse.bed)

      continue = ifelse(test = nrow(collapsed.bed) == nrow(to.collapse.bed),
                        yes = F,
                        no = T)

      to.collapse.bed = Rseb::sort.bed(collapsed.bed, verbose = F)
    }
  }


  # Remove the maximal.distance at the extremities added at the beginning
  collapsed.bed$start = collapsed.bed$start + maximal.distance
  collapsed.bed$end = collapsed.bed$end - maximal.distance


  ### -------------------------------------------------------------------------------- ###
  ###                                      Export                                      ###
  ### -------------------------------------------------------------------------------- ###
  # Export bed if required
  if (!is.null(export.file.name) & (verbose == T)) {
    write.table(x = unique(collapsed.bed[,1:final.ncol]),
                file = export.file.name,
                quote = F,
                sep = sep,
                col.names = export.header)

    message(paste("The collapsed bed file has been exported as",
                  export.file.name))
  }

  # Return the bed
  if (return.bed == T) {return(format(unique(collapsed.bed[,1:final.ncol]), scientific = FALSE))}
}
