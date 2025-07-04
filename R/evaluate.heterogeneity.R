#' @title Evaluate genomic heterogeneity among samples.
#'
#' @description This tools evaluates what is the fraction of peaks covered by each sample provided in a peaks dataset obtained by merging all the peaks together or provided by the user. The peaks in the reference dataset are ranked by number of samples in which are present and average score all over the samples. This function uses the deeptools function multiBigwigSummary.
#'
#' @param bigWig.list A string vector with bigwig paths (same order than paths).
#' @param peak.list A list of GRanges objects (not GRglist) or data.frames or a string vector with the path to bed files (same order than bigwigs).
#' @param labels The labels to use for the samples (same order than bigwigs/peaks). Default: the basename of the bigWig.list.
#' @param reference.peaks Default: \code{NULL}, the peaks of all samples provided are merged and collapsed together.
#' @param distribution.line.color Color to use for the distribution line. Default: \code{"#1c30a3"} (dark blue).
#' @param distribution.line.size Line size of the distribution plot. Default: \code{1}.
#' @param distribution.line.type Line type of the distribution plot. Default: \code{1}.
#' @param distribution.n.vertical.divisions Number of sectors in which divide the distribution plot (vertical lines will be plotted). Default: \code{NULL} (no divisions).
#' @param distribution.as.percentage Logical value to define whether the distribution plot should show percentage of sample coverage rather than number of samples. Default: \code{FALSE}.
#' @param heatmap.color Color to use for the heatmaps; a gradient from this color to white will be used. Default: \code{"#1c30a3"} (dark blue).
#' @param heatmap.zMax Maximum of the heatmap scale. Default: \code{NA}.
#' @param heatmap.log1p.scale Logic value to define whether the heatmap scale should display log1p values. Default: \code{TRUE}.
#' @param bar.color Color to use for the barplot showing the fraction of reference peaks present in each sample. Default is to use the 'distribution.line.color'.
#' @param widths.proportion Two-elements numeric vector to be passed to \code{plot_grid} rel_width.
#' @param heights.proportion Two-elements numeric vector to be passed to \code{plot_grid} rel_height.
#' @param min.percentage.reference Numeric value within 0-100 to define which percentage of 'reference' dataset must overlap with a 'sample'. If the value is lower than 0 or greater than 100, will be coerced to 0 or 100 respectively. Default value: \code{0}.
#' @param min.percentage.test Numeric value within 0-100 to define which percentage of 'sample' datasets must overlap with a region in the 'reference' one. If the value is lower than 0 or greater than 100, will be coerced to 0 or 100 respectively. Default value: \code{0}.
#' @param min.bases.overlap Integer, greater than 0, value to indicate the minimal number of bases to consider as minimum overlap between two regions. Non integer values will be rounded at integer, while number lower that 1 will be coerced to 1. Default value: \code{1}.
#' @param multiBigWigSummary.path Path/command to run deeptools multiBigWigSummary tool. Default: \code{"multiBigWigSummary"}.
#'
#' @return The function returns a list containing:
#' \itemize{
#'   \item \code{metadata}: a list with the following elements \itemize{
#'                          \item \code{sample.id}: vector with the labels of the samples
#'                          \item \code{bigWig.file}: string vector with the path to each bigwig file
#'                          \item \code{bed.file}: list of peaks associated to each sample
#'                          \item \code{n.peaks}: vector with number of peaks in each sample
#'                          }
#'   \item \code{count.matrix}: data.frame with presence (1) or absence (0) of each peak per each sample;
#'   \item \code{score.matrix}: data.frame with average score at each peak per each sample;
#'   \item \code{plot.list}: list of single separated plots: counts.distribution, fraction.peaks.per.sample, scores.heatmap;
#'   \item \code{multiplot}: the multiplot generated from the plot.list.
#'  }
#'
#'
#' @details To know more about the deepTools's function \code{multiBigwigSummary} see the package manual at the following link: \cr \url{https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html}.
#'
#' @export evaluate.heterogeneity



evaluate.heterogeneity = function(bigWig.list,
                                  peak.list,
                                  labels = sub(pattern = "(.*)\\..*$",
                                               replacement = "\\1",
                                               basename(bigWig.list)),
                                  reference.peaks = NULL,
                                  distribution.line.color = "#1c30a3",
                                  distribution.line.size = 1,
                                  distribution.line.type = 1,
                                  distribution.n.vertical.divisions = NULL,
                                  distribution.as.percentage = F,
                                  heatmap.color = "#1c30a3",
                                  heatmap.zMax = NA,
                                  heatmap.log1p.scale = TRUE,
                                  bar.color = distribution.line.color,
                                  widths.proportion = c(0.25,1),
                                  heights.proportion = c(1,1),
                                  min.percentage.reference = 0 ,
                                  min.percentage.test = 0,
                                  min.bases.overlap = 1,
                                  multiBigWigSummary.path = "multiBigWigSummary") {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#
  
  
  require(IRanges)
  require(GenomicRanges)
  require(dplyr)
  require(ggplot2)
  #require(scales)


  # Check the length of the variables
  if (length(unique(c(length(unique(bigWig.list)), length(unique(peak.list)), length(unique(labels))))) != 1) {
    return(warning("The length of elements in 'bigWig.list', 'peak.list' and 'labels' does not correpond."))
  }


  # Reading the bed files
  bed.list = list()

  if ("character" %in% class(peak.list)) {
    for (i in 1:length(peak.list)) {
      bed.list[[i]] = data.frame(data.table::fread(peak.list[i]), stringsAsFactors = F)[,1:3]
      names(bed.list[[i]])[1:3] = c("seqnames", "start", "end")
      bed.list[[i]] = IRanges::reduce(GenomicRanges::makeGRangesFromDataFrame(bed.list[[i]]))
    }
  } else if ("list" %in% class(peak.list)) {
    for (i in 1:length(peak.list)) {
      if ("data.frame" %in% class(peak.list[[i]])) {
        names(peak.list[[i]])[1:3] = c("seqnames", "start", "end")
        bed.list[[i]] = IRanges::reduce(GenomicRanges::makeGRangesFromDataFrame(as.data.frame(peak.list[[i]])))
      } else if ("GRanges" %in% class(peak.list[[i]])){
        bed.list[[i]] = IRanges::reduce(peak.list[[i]])
      }
    }
  }

  names(bed.list) = labels


  # Combine all peaks in a unique file
  if (is.null(reference.peaks)) {
    all.peaks.merge = IRanges::reduce(GenomicRanges::makeGRangesFromDataFrame(do.call(rbind, lapply(bed.list, data.frame))))
  } else {
    if ("data.frame" %in% class(reference.peaks)) {
      names(reference.peaks)[1:3] = c("seqnames", "start", "end")
      all.peaks.merge = IRanges::reduce(GenomicRanges::makeGRangesFromDataFrame(as.data.frame(reference.peaks)))
    } else if ("character" %in% class(reference.peaks)) {
      all.peaks.merge = IRanges::reduce(rtracklayer::import.bed(reference.peaks))
    } else if ("GRanges" %in% class(reference.peaks)) {
      all.peaks.merge = IRanges::reduce(reference.peaks)
    } else {
      return(warning("The 'reference.peaks' cannot be read. Accepted formats: 'GRanges', 'character' (path), 'data.frame'/'data.table'"))
    }
  }




  # Find overlaps in each sample
  overlaps.counts.table = dplyr::mutate(as.data.frame(all.peaks.merge)[,1:3],
                                        seqnames = gsub("chr","",seqnames))

  for (i in 1:length(bed.list)) {
    regions.with.overlap = (Rseb::intersect.regions(reference.regions = all.peaks.merge,
                                                    test.regions = bed.list[[i]],
                                                    min.percentage.reference = min.percentage.reference,
                                                    min.percentage.test = min.percentage.test,
                                                    min.bases.overlap = min.bases.overlap))$overlaps.reference

    overlaps.counts.table = dplyr::left_join(x = overlaps.counts.table,
                                             y = dplyr::mutate(as.data.frame(regions.with.overlap)[,1:3],
                                                               seqnames = gsub("chr","",seqnames),
                                                               "{labels[i]}" := 1),
                                             by = c("seqnames", "start", "end"))
  }

  overlaps.counts.table[is.na(overlaps.counts.table)] = 0

  # Computing fraction of samples that show a given peak
  overlaps.counts.table =
    overlaps.counts.table %>%
    dplyr::mutate(fraction.samples = rowSums(overlaps.counts.table[,-c(1:3)], na.rm = T) / (ncol(overlaps.counts.table)-3)) %>%
    dplyr::arrange(desc(fraction.samples))



  # Compute the average score at each peak
  overlaps.scores.table = dplyr::mutate(as.data.frame(all.peaks.merge)[,1:3],
                                        seqnames = gsub("chr","",seqnames))

  # write the temporary merged bed
  write.table(x = as.data.frame(all.peaks.merge)[,1:3],
              file = "~/temp_merge_all_peaks.bed",
              sep = "\t", quote = F, col.names = F, row.names = F)

  # Compute Scores at bed region
  system(paste(multiBigWigSummary.path, "BED-file",
               "-b", paste(bigWig.list, collapse = " "),
               "-o ~/temp_matrix.npz",
               "--outRawCounts ~/temp_matrix.txt",
               "--BED ~/temp_merge_all_peaks.bed",
               "--labels", paste(labels, collapse = " "),
               "-p max"))

  # Import scores multiBigWig matrix
  score.matrix = data.table::fread("~/temp_matrix.txt")
  colnames(score.matrix) = gsub("#|'", "", colnames(score.matrix))
  colnames(score.matrix)[1] = "seqnames"

  score.matrix =
    dplyr::left_join(x = score.matrix %>% dplyr::mutate(average.score = rowMeans(as.matrix(score.matrix[,-c(1:3)]), na.rm = T)),
                     y = overlaps.counts.table[,c(1:3,ncol(overlaps.counts.table))],
                     by = c("seqnames", "start", "end")) %>%
    dplyr::arrange(desc(fraction.samples), desc(average.score)) %>%
    dplyr::filter(average.score > 0)

  score.matrix = dplyr::mutate(score.matrix, rank = 1:nrow(score.matrix))


  # adding the rank to the counts
  overlaps.counts.table =
    dplyr::left_join(x = overlaps.counts.table,
                     y = data.frame(score.matrix)[,c(1:3,ncol(score.matrix))],
                     by = c("seqnames", "start", "end")) %>%
    dplyr::arrange(rank)



  ######################
  # Creating the plots #
  ######################

  ### reshaping the tables
  # Counts per sample
  number.peaks.per.sample = colSums(overlaps.counts.table[,-c(1:3,(ncol(overlaps.counts.table)-1),ncol(overlaps.counts.table))])
  total.counts.per.sample =
    data.frame(sample = factor(names(number.peaks.per.sample), levels = rev(labels)),
               n.peaks = number.peaks.per.sample,
               fraction.total.peaks = number.peaks.per.sample / nrow(data.frame(overlaps.counts.table)))

  # Scores per sample
  scores.heatmap =
    reshape2::melt(data = score.matrix,
                   id = c("seqnames", "start", "end", "rank", "fraction.samples", "average.score"),
                   value.name = "signal") %>%
    dplyr::mutate(variable = factor(variable, levels = rev(labels))) %>%
    dplyr::rename(sample = variable)



  ####### Building the plots
  # Distribution plot
  if (isFALSE(distribution.as.percentage)) {
    peak.presence.distribution =
      ggplot(data = overlaps.counts.table,
             aes(x = rank,
                 y = fraction.samples * length(labels)))
  } else {
    peak.presence.distribution =
      ggplot(data = overlaps.counts.table,
             aes(x = rank,
                 y = fraction.samples * 100))
  }

  if (!is.null(distribution.n.vertical.divisions)) {
    division.step = (1/distribution.n.vertical.divisions) * nrow(score.matrix)

    peak.presence.distribution =
      peak.presence.distribution +
      geom_vline(xintercept = seq(division.step,(nrow(score.matrix)-division.step), division.step),
                 linetype = "dotted",
                 color = "gray70")
  }

  peak.presence.distribution =
    peak.presence.distribution +
    geom_line(color = distribution.line.color,
              linewidth = distribution.line.size,
              linetype = distribution.line.type,
              show.legend = F) +
    ylab(ifelse(distribution.as.percentage, yes = "% Samples", no = "Sample number")) +
    xlab(NULL) +
    scale_y_continuous(breaks = scales::pretty_breaks(),
                       limits = c(0,
                                  ifelse(distribution.as.percentage, yes = 100, no = length(labels)))) +
    coord_cartesian(expand = F) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(color = "#000000"),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(color = "#000000"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = "#000000"),
          panel.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), "lines"))



  peak.signal.heatmap =
    ggplot(data = scores.heatmap,
           aes(x = rank,
               y = sample,
               fill = signal)) +
    geom_tile(na.rm = T)
  
  if (heatmap.log1p.scale == TRUE) {
    peak.signal.heatmap =
      peak.signal.heatmap +
      scale_fill_gradient(name = "log1p(signal)",
                          low = "white",
                          high = heatmap.color,
                          limits = c(0, heatmap.zMax),
                          na.value = heatmap.color,
                          trans = "log1p")
  } else {
    peak.signal.heatmap =
      peak.signal.heatmap +
      scale_fill_gradient(name = "Signal",
                          low = "white",
                          high = heatmap.color,
                          limits = c(0, heatmap.zMax),
                          na.value = heatmap.color)
  }
  
  peak.signal.heatmap =
    peak.signal.heatmap +
    ylab(NULL) +
    scale_y_discrete(position = "right") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    coord_cartesian(expand = F) +
    xlab("Rank") +
    theme_classic() +
    theme(axis.text = element_text(color = "#000000"),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(color = "#000000"),
          panel.border = element_rect(colour = "#000000", fill=NA, linewidth=1),
          panel.background = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), "lines"))

  heatmap.rank.breaks = c(1, ggplot_build(peak.signal.heatmap)$layout$panel_params[[1]]$x.sec$breaks, nrow(scores.heatmap))
  heatmap.rank.breaks = sort(heatmap.rank.breaks[!is.na(heatmap.rank.breaks)])

  peak.signal.heatmap =
    peak.signal.heatmap +
    scale_x_continuous(breaks = heatmap.rank.breaks,
                       labels = scales::comma)



  fraction.peaks.per.sample =
    ggplot(total.counts.per.sample,
           aes(x = fraction.total.peaks,
               y = sample)) +
    geom_bar(data = data.frame(sample = factor(names(number.peaks.per.sample), levels = rev(labels)),
                               total = 1),
             aes(x = total,
                 y = sample),
             position = "stack",
             stat = "identity",
             inherit.aes = F,
             show.legend = F,
             color = NA,
             fill = colorspace::lighten(bar.color, amount = 0.925)) +
    geom_bar(position="stack",
             stat="identity",
             show.legend = F,
             fill = bar.color) +
    geom_bar(data = data.frame(sample = factor(names(number.peaks.per.sample), levels = rev(labels)),
                               total = 1),
             aes(x = total,
                 y = sample),
             position = "stack",
             stat = "identity",
             inherit.aes = F,
             show.legend = F,
             color = "#000000",
             fill = NA) +
    ylab(NULL) +
    xlab("Fraction of all peaks") +
    scale_x_reverse(expand = c(0, 0), limits = c(1,0)) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(color = "#000000"),
          panel.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), "lines"))


  # Assembly of the multiplot
  multiplot = cowplot::plot_grid(plotlist = list(NULL,
                                                 peak.presence.distribution,
                                                 fraction.peaks.per.sample,
                                                 peak.signal.heatmap),
                                 rel_widths = widths.proportion,
                                 rel_heights = heights.proportion,
                                 align = "hv",
                                 axis = "lrtb")


  # Remove temporary files
  invisible(file.remove("~/temp_merge_all_peaks.bed"))
  invisible(file.remove("~/temp_matrix.txt"))
  invisible(file.remove("~/temp_matrix.npz"))


  # Build the output
  return(list(metadata = list(sample.id = labels,
                              bigWig.file = bigWig.list,
                              bed.file = peak.list,
                              n.peaks = as.vector(sapply(bed.list, function(x){return(nrow(data.frame(x)))}))),
              count.matrix = dplyr::arrange(overlaps.counts.table, rank),
              score.matrix = dplyr::arrange(score.matrix, rank),
              plot.list = list(counts.distribution = peak.presence.distribution,
                               fraction.peaks.per.sample = fraction.peaks.per.sample,
                               scores.heatmap = peak.signal.heatmap),
              multiplot = multiplot))

} # END function
