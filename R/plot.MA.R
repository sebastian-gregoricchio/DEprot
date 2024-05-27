#' @title plot.MA
#'
#' @description Plots a MA plot log2(basemean) x log2(FoldChange) of differential expression results
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast Number indicating the position of the contrast to use for the plotting.
#' @param up.color String indicating the color to use for up-regulated proteins in the plots. Default: \code{"indianred"}.
#' @param down.color String indicating the color to use for up-regulated proteins in the plots. Default: \code{"steelblue"}.
#' @param density.colors List of colors, passed to \code{scale_fill_gradientn}, to use for the density gradient. Default: \code{"colorRampPalette(colors = RColorBrewer::brewer.pal(9, "Blues"))(101)"}.
#' @param point.size Numeric value indicating the size of the dots. Default: \code{2}.
#' @param point.alpha Numeric value between 0 and 1 to indicate the transparency (alpha) of the dots. Default: \code{0.5}.
#' @param title String indicating the title to use. Default: \code{NULL} (automatic title).
#' @param use.uncorrected.pvalue Logical value indicating whether it should be used the normal p-value instead of the adjusted one (differential proteins numbers are recomputed). Default: \code{FALSE}, padj is used.
#' @param symmetric.x Logical values indicating whether the x-axis scale should be symmetric or not. Default: \code{TRUE}.
#'
#' @return A ggplot object.
#'
#' @export plot.MA

plot.MA =
  function(DEprot.analyses.object,
           contrast = 1,
           up.color = "indianred",
           down.color = "steelblue",
           density.colors = colorRampPalette(colors = RColorBrewer::brewer.pal(9, "Blues"))(101),
           point.size = 2,
           point.alpha = 0.5,
           title = NULL,
           use.uncorrected.pvalue = FALSE,
           symmetric.y = TRUE) {

    ### Libraries
    require(dplyr)
    require(ggplot2)


    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      return(warning("The input must be an object of class 'DEprot.analyses'."))
    }

    ### check and collect contrast
    if (is.numeric(contrast)) {
      if (contrast <= length(DEprot.analyses.object@analyses.result.list)) {
        data = DEprot.analyses.object@analyses.result.list[[contrast]]$results
        n.diff = DEprot.analyses.object@analyses.result.list[[contrast]]$n.diff
        contrasts.info = DEprot.analyses.object@contrasts[[contrast]]
      } else {
        return(warning("The 'contrast' indicated is not available."))
      }
    } else {
      return(warning("The 'contrast' must be a numeric value."))
    }


    ## Define params and counts
    colors.plots = c(up.color, down.color)
    names(colors.plots) = c(contrasts.info$var.1, contrasts.info$var.2)


    ### Make a table for plotting
    diff.tb = data.frame(log2.Fold_group1.vs.group2 = data[,5],
                         padj = data$padj,
                         basemean.log2 = data$basemean.log2,
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
               levels = c(names(colors.plots), "unresponsive", "null"))

      n.diff =
        data.frame(diff.tb %>%
                     dplyr::group_by(diff.status, .drop = F) %>%
                     dplyr::summarise(n = n(),
                                      median.FoldChange = median(log2.Fold_group1.vs.group2)))
    }



    ### Make MA-plot
    ma.plot =
      ggplot() +
      stat_density_2d(data = diff.tb,
                      aes(x = basemean.log2,
                          y = log2.Fold_group1.vs.group2,
                          fill = after_stat(count)),
                      geom = "raster",
                      contour = FALSE,
                      show.legend = T,
                      n = 200,
                      adjust = 5) +
      scale_fill_gradientn(colours = density.colors, name = "Count")


    if (sum((n.diff %>% dplyr::filter(!(diff.status %in% c("unresponsive", "null"))))$n, na.rm = T) > 0) {
      ma.plot =
        ma.plot +
        geom_point(data = dplyr::filter(diff.tb, diff.status %in% c(contrasts.info$var.1,contrasts.info$var.2)),
                   mapping = aes(x = basemean.log2,
                                 y = log2.Fold_group1.vs.group2,
                                 color = diff.status),
                   alpha = point.alpha,
                   size = point.size,
                   stroke = NA,
                   show.legend = T,
                   inherit.aes = F) +
        scale_color_manual(values = colors.plots, name = "Differential\nstatus", drop = FALSE)
    }

    ma.plot =
      ma.plot +
      geom_hline(yintercept = c(-1,1)*log2(DEprot.analyses.object@differential.analyses.params$linear.FC.th), linetype = 2, color = "gray40") +
      geom_hline(yintercept = 0, linetype = 1, color = "steelblue") +
      theme_classic() +
      xlab("log~2~(Base Mean)") +
      ylab(paste0("log~2~(Fold Change<sub>",contrasts.info$var.1,"</sup>&frasl;<sub>",contrasts.info$var.2,"</sub>)")) +
      ggtitle(ifelse(is.null(title), yes = paste0("**",contrasts.info$var.1, "** *vs* **", contrasts.info$var.2, "**"), no = title)) +
      scale_x_continuous(expand = c(0,0)) +
      annotate(geom = "text",
               x = -Inf, y = -Inf,
               color = colors.plots[2],
               hjust = -0.2, vjust = -0.5,
               label = paste0("n = ",n.diff[n.diff$diff.status == contrasts.info$var.2,"n"])) +
      annotate(geom = "text",
               x = -Inf, y = +Inf,
               hjust = -0.2, vjust = 1.5,
               color = colors.plots[1],
               label = paste0("n = ",n.diff[n.diff$diff.status == contrasts.info$var.1,"n"])) +
      theme(axis.ticks = element_line(color = "black"),
            axis.text = element_text(color = "black"),
            axis.title = ggtext::element_markdown(color = "black"),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            aspect.ratio = 0.6)


    if (symmetric.y == TRUE) {
      y.max = max(abs(ggplot_build(ma.plot)$layout$panel_params[[1]]$y.range))
      ma.plot = ma.plot + scale_y_continuous(expand = c(0,0), limits = c(-1,1)*y.max)
    }

    ### Export plot
    return(ma.plot)
  } # END of function

