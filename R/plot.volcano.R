#' @title plot.volcano
#'
#' @description Plots a volcano plot log2(FoldChange) x -log10(p.adjusted) of differential expression results
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast Number indicating the position of the contrast to use for the plotting.
#' @param up.color String indicating the color to use for up-regulated proteins in the plots. Default: \code{"indianred"}.
#' @param down.color String indicating the color to use for up-regulated proteins in the plots. Default: \code{"steelblue"}.
#' @param unresponsive.color String indicating the color to use for unresponsive proteins in the plots. Default: \code{"purple"}.
#' @param null.color String indicating the color to use for null proteins in the plots. Default: \code{"gray"}.
#' @param point.size Numeric value indicating the size of the dots. Default: \code{2}.
#' @param point.alpha Numeric value between 0 and 1 to indicate the transparency (alpha) of the dots. Default: \code{0.5}.
#' @param title String indicating the title to use. Default: \code{NULL} (automatic title).
#' @param use.uncorrected.pvalue Logical value indicating whether it should be used the normal p-value instead of the adjusted one (differential proteins numbers are recomputed). Default: \code{FALSE}, padj is used.
#' @param symmetric.x Logical values indicating whether the x-axis scale should be symmetric or not. Default: \code{TRUE}.
#' @param dot.labels String vector indicating labels to show on the plot that should correspond to \code{prot.id} column values. Default: \code{NULL} (no labels shown).
#' @param protein.names.pattern Character indicating a regular expression to remove from the protein IDs. Default: \code{""}, no alterations in the protein IDs.
#' @param labels.in.boxes Logical value indicating whether the labels should be visualized as boxes. Default: \code{FALSE}.
#' @param label.font.size Numeric value indicating the size to use for the dot labels. Default: \code{2}.
#' @param label.max.overlaps Numeric value indicating the maximum number of overlaps allowed between labels. Default: \code{100}.
#' @param min.segment.length.labels Numeric value indicating the minimal length of the segments that connect the labels to the points. Default: \code{0} (segment always shown).
#'
#' @return A ggplot object.
#'
#' @name plot.volcano
#'
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import ggtext
#' @importFrom purrr pmap
#'
#' @author Sebastian Gregoricchio
#'
#' @export plot.volcano

plot.volcano =
  function(DEprot.analyses.object,
           contrast = 1,
           up.color = "indianred",
           down.color = "steelblue",
           unresponsive.color = "purple",
           null.color = "gray",
           point.size = 2,
           point.alpha = 0.5,
           title = NULL,
           use.uncorrected.pvalue = FALSE,
           symmetric.x = TRUE,
           dot.labels = NULL,
           protein.names.pattern = "",
           labels.in.boxes = FALSE,
           label.font.size = 2,
           label.max.overlaps = 100,
           min.segment.length.labels = 0) {

    # ### Libraries
    # require(dplyr)
    # require(ggplot2)


    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      stop("The input must be an object of class 'DEprot.analyses'.")
      #return()
    } else {
      # identify whether the analyses have been performed using prolfqua
      is.prolfqua = "strategy" %in% names(DEprot.analyses.object@differential.analyses.params)

      if (isTRUE(is.prolfqua)) {
        padj.label = "-log<sub>10</sub>(*FDR*)"
      } else {
        padj.label = "-log<sub>10</sub>(*P<sub>adj</sub>*)"
      }
    }

    ### check and collect contrast
    if (is.numeric(contrast)) {
      if (contrast <= length(DEprot.analyses.object@analyses.result.list)) {
        data = DEprot.analyses.object@analyses.result.list[[contrast]]$results
        n.diff = DEprot.analyses.object@analyses.result.list[[contrast]]$n.diff
        contrasts.info = DEprot.analyses.object@contrasts[[contrast]]

        # Change the name for 'FDR' column into 'padj' for prolfqua analyses
        if (isTRUE(is.prolfqua)) {data = data %>% dplyr::rename(padj = FDR)}

      } else {
        stop("The 'contrast' indicated is not available.")
        #return()
      }
    } else {
      stop("The 'contrast' must be a numeric value.")
      #return()
    }


    ## Define params and counts
    colors.plots = c(up.color, down.color, unresponsive.color, null.color)
    names(colors.plots) = c(contrasts.info$var.1, contrasts.info$var.2, "unresponsive", "null")


    ### Make a table for plotting
    diff.tb = data.frame(prot.id = data$prot.id,
                         log2.Fold_group1.vs.group2 = data[,5],
                         padj = data$padj,
                         diff.status = factor(data$diff.status, levels = names(colors.plots)))



    ### Define pvalue
    if (use.uncorrected.pvalue == TRUE) {
      diff.tb$padj = data$p.value

      de.status =
        function(FC, p) {
          ifelse(p < DEprot.analyses.object@differential.analyses.params$padj.th,
                 yes = ifelse(abs(FC) >= log2(DEprot.analyses.object@differential.analyses.params$linear.FC.th),
                              yes = ifelse(sign(FC) == 1,
                                           yes = contrasts.info$var.1,
                                           no = contrasts.info$var.2),
                              no = "null"),
                 no = ifelse(FC >= log2(DEprot.analyses.object@differential.analyses.params$linear.FC.unresp.range[1]) & FC <= log2(DEprot.analyses.object@differential.analyses.params$linear.FC.unresp.range[2]),
                             yes = "unresponsive",
                             no = "null"))
        }

      diff.tb$diff.status =
        factor(unlist(purrr::pmap(.l = list(FC = diff.tb$log2.Fold_group1.vs.group2,
                                            p = diff.tb$padj),
                                  .f = function(FC,p){de.status(FC,p)})),
               levels = names(colors.plots))

      n.diff =
        data.frame(diff.tb %>%
                     dplyr::group_by(diff.status, .drop = FALSE) %>%
                     dplyr::summarise(n = n(),
                                      median.FoldChange = median(log2.Fold_group1.vs.group2)))
    }


    ### Make volcano
    volcano =
      ggplot(diff.tb,
             aes(x = log2.Fold_group1.vs.group2,
                 y = -log10(padj),
                 color = diff.status)) +
      geom_point(stroke = NA,
                 alpha = point.alpha,
                 size = point.size) +
      scale_color_manual(values = colors.plots,
                         name = "Differential\nstatus",
                         drop = FALSE) +
      geom_hline(yintercept = -log10(DEprot.analyses.object@differential.analyses.params$padj.th), linetype = 2, color = "gray40") +
      geom_vline(xintercept = c(-1,1)*log2(DEprot.analyses.object@differential.analyses.params$linear.FC.th), linetype = 2, color = "gray40") +
      ylab(ifelse(use.uncorrected.pvalue == FALSE, yes = padj.label, no = "-log~10~(*P*)")) +
      xlab(paste0("log~2~(Fold Change<sub>",contrasts.info$var.1,"</sup>&frasl;<sub>",contrasts.info$var.2,"</sub></sub>)")) +
      ggtitle(ifelse(is.null(title), yes = paste0("**",contrasts.info$var.1, "** *vs* **", contrasts.info$var.2, "**"),no = title)) +
      guides(color = guide_legend(override.aes = list(size = max(c(point.size, 3))))) +
      theme_classic() +
      theme(axis.text = ggtext::element_markdown(color = "black"),
            axis.title = ggtext::element_markdown(color = "black"),
            axis.ticks = element_line(color = "black"),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            aspect.ratio = 1)



    # add labels if required
    if (!is.null(dot.labels)) {
      filt.tb = dplyr::filter(.data = diff.tb, prot.id %in% dot.labels)

      if (nrow(filt.tb) > 0) {
        if (labels.in.boxes == TRUE) {
          volcano =
            volcano +
            ggrepel::geom_label_repel(data = filt.tb,
                                      aes(label = gsub(protein.names.pattern,"",prot.id)),
                                      color = "black",
                                      show.legend = F,
                                      min.segment.length = min.segment.length.labels,
                                      max.overlaps = label.max.overlaps,
                                      size = label.font.size)
        } else {
          volcano =
            volcano +
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



    # make x-axis symmetric (if required)
    if (symmetric.x == TRUE) {
      x.max = max(abs(ggplot_build(volcano)$layout$panel_params[[1]]$x.range))

      volcano = volcano + xlim(c(-1,1)*x.max)
    }

    volcano =
      volcano +
      annotate(geom = "text",
               x = -Inf, y = +Inf,
               color = colors.plots[2],
               hjust = -0.2, vjust = 1.5,
               label = paste0("n = ",n.diff[n.diff$diff.status == contrasts.info$var.2,"n"])) +
      annotate(geom = "text",
               x = +Inf, y = +Inf,
               hjust = 1.2, vjust = 1.5,
               color = colors.plots[1],
               label = paste0("n = ",n.diff[n.diff$diff.status == contrasts.info$var.1,"n"]))


    ### Export plot
    return(volcano)
  } # END of function

