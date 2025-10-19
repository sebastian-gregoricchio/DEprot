#' @title contrast.LFQ
#'
#' @description Plots a scatter of the average LFQ values of each condition in a contrast. Dots will be colored by FoldChange and the size will correspond to the -log10(Padj).
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast Integer indicating the position of the contrast to use. Default: \code{1}.
#' @param show.only.differential Logical value indicating whether only up or down regulated proteins should be plotted. Default: \code{FALSE}.
#' @param show.only.significant Logical value indicating whether only proteins with a significant p-value adjusted should be plotted. Notice that this is subordinated to `show.only.differential` (differential proteins pass the P-adj threshold by default). Default: \code{FALSE}.
#' @param log2FC.positive.color String indicating any R-supported color for the positive foldchage values. Default: \code{"firebrick"}.
#' @param log2FC.negative.color String indicating any R-supported color for the negative foldchage values. Default: \code{"steelblue4"}.
#' @param log2FC.scale.min Numeric value indicating the minimum of the foldchange scale (lower values will be colored with the minimum of the scale). Default: \code{NULL}.
#' @param log2FC.scale.max Numeric value indicating the maximum of the foldchange scale (lower values will be colored with the maximum of the scale). Default: \code{NULL}.
#' @param identical.axes Logical value indicating whether the axes should be forces to be identical between x and y. Default: \code{TRUE}.
#' @param dot.labels String vector indicating labels to show on the plot that should correspond to \code{prot.id} column values. Default: \code{NULL} (no labels shown).
#' @param protein.names.pattern Character indicating a regular expression to remove from the protein IDs. Default: \code{""}, no alterations in the protein IDs.
#' @param labels.in.boxes Logical value indicating whether the labels should be visualized as boxes. Default: \code{FALSE}.
#' @param label.font.size Numeric value indicating the size to use for the dot labels. Default: \code{3}.
#' @param label.max.overlaps Numeric value indicating the maximum number of overlaps allowed between labels. Default: \code{100}.
#' @param min.segment.length.labels Numeric value indicating the minimal length of the segments that connect the labels to the points. Default: \code{0} (segment always shown).
#'
#'
#'
#' @return A scatter plot of class ggplot2.
#'
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#' @import ggrepel
#' @importFrom ggpubr theme_pubr
#'
#' @examples
#' contrast.LFQ(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma,
#'              contrast = 1)
#'
#'
#' @export contrast.LFQ



contrast.LFQ =
  function(DEprot.analyses.object,
           contrast = 1,
           show.only.differential = FALSE,
           show.only.significant = FALSE,
           log2FC.positive.color = "firebrick",
           log2FC.negative.color = "steelblue4",
           log2FC.scale.min = NULL,
           log2FC.scale.max = NULL,
           identical.axes = TRUE,
           dot.labels = NULL,
           protein.names.pattern = "",
           labels.in.boxes = FALSE,
           label.font.size = 3,
           label.max.overlaps = 100,
           min.segment.length.labels = 0) {

    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      stop("The input must be an object of class 'DEprot.analyses'.")
      #return(invisible())
    }

    ### check and collect contrast
    if (!(contrast %in% 1:length(DEprot.analyses.object@analyses.result.list))) {
      stop(paste0("The `contrast` indicated (",contrast,") is not available in the DEprot.analyses object."))
    } else {
      results = DEprot::get.results(DEprot.analyses.object = DEprot.analyses.object, contrast = contrast)
      contrast.info = DEprot.analyses.object@contrasts[[contrast]]

      if (DEprot.analyses.object@differential.analyses.params$padj.method == "effective FDR") {
        results  = dplyr::mutate(.data = results, padj = FDR)
        size.label = "-log<sub>10</sub>(FDR)"
      } else {
        size.label = "-log<sub>10</sub>(P<sub>adj</sub>)"
      }
    }


    ### Make table for plotting
    # build table
    tb = data.frame(prot.id = results$prot.id,
                    lfq.A = results[,3],
                    lfq.B = results[,4],
                    log2FC = results[,5],
                    pval = results$padj,
                    status = results$diff.status)

    # filter not significant, if needed
    if (show.only.differential == TRUE) {
      tb = dplyr::filter(.data = tb, !(status %in% c("null", "unresponsive")))
    } else if (show.only.significant == TRUE) {
      tb = dplyr::filter(.data = tb, pval < DEprot.analyses.object@differential.analyses.params$padj.th)
    }

    # adjust the over scale values
    log2FC.scale = range(tb$log2FC)

    if (!is.null(log2FC.scale.min)) {
      tb$log2FC[tb$log2FC < log2FC.scale.min] = log2FC.scale.min
      log2FC.scale[1] = log2FC.scale.min
    }

    if (!is.null(log2FC.scale.max)) {
      tb$log2FC[tb$log2FC > log2FC.scale.max] = log2FC.scale.max
      log2FC.scale[2] = log2FC.scale.max
    }



    ### Make the plot
    lfq_plot =
      ggplot(data = tb,
             aes(x = lfq.B,
                 y = lfq.A,
                 color = log2FC,
                 size = -log10(pval))) +
      geom_point(alpha = 0.75) +
      scale_color_gradient2(name = "log<sub>2</sub>(Fold Change)",
                            low = log2FC.negative.color,
                            mid = "white",
                            high = log2FC.positive.color,
                            midpoint = 0,
                            limits = log2FC.scale) +
      labs(size = size.label) +
      xlab(paste0("log<sub>2</sub>(LFQ ", contrast.info$var.2, ")")) +
      ylab(paste0("log<sub>2</sub>(LFQ ", contrast.info$var.1, ")")) +
      ggtitle(label = paste0("**",contrast.info$metadata.column,"**"),
              subtitle = paste0(contrast.info$var.1, " *vs* ", contrast.info$var.2)) +
      geom_abline(slope = 1, intercept = 0, colour = "gray", linetype = "dashed") +
      ggpubr::theme_pubr(legend = "right") +
      theme(aspect.ratio = 1,
            axis.title.x = ggtext::element_markdown(),
            axis.title.y = ggtext::element_markdown(),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5),
            legend.title = ggtext::element_markdown())


    # Add symmetrric axises, if needed
    if (identical.axes == TRUE) {
      build = ggplot_build(lfq_plot)
      max = max(c(build$layout$panel_params[[1]]$x.range, build$layout$panel_params[[1]]$y.range))
      min = min(c(build$layout$panel_params[[1]]$x.range, build$layout$panel_params[[1]]$y.range))
      lfq_plot = lfq_plot + scale_x_continuous(limits = c(min,max)) + scale_y_continuous(limits = c(min,max))
    }


    # add labels if required
    if (!is.null(dot.labels)) {
      filt.tb = dplyr::filter(.data = tb, prot.id %in% dot.labels)

      if (nrow(filt.tb) > 0) {
        if (labels.in.boxes == TRUE) {
          lfq_plot =
            lfq_plot +
            ggrepel::geom_label_repel(data = filt.tb,
                                      aes(label = gsub(protein.names.pattern,"",prot.id)),
                                      color = "black",
                                      show.legend = F,
                                      min.segment.length = min.segment.length.labels,
                                      max.overlaps = label.max.overlaps,
                                      size = label.font.size)
        } else {
          lfq_plot =
            lfq_plot +
            ggrepel::geom_text_repel(data = filt.tb,
                                     aes(label = gsub(protein.names.pattern,"",prot.id)),
                                     color = "black",
                                     show.legend = FALSE,
                                     min.segment.length = min.segment.length.labels,
                                     max.overlaps = label.max.overlaps,
                                     size = label.font.size)
        }
      }
    }


    ### Return the object
    return(lfq_plot)

  } # END function













