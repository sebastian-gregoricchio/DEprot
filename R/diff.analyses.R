#' @title diff.analyses
#'
#' @description Allows for the computation of differential analyses using multiple t/Wilcox tests. Includes means, Fold Changes, and p-values.
#'
#' @param DEprot.object An object of class \code{DEprot}.
#' @param contrast.list List of 3-elements vectors indicating (in order): metadata_column, variable_1, variable_2.
#' @param replicate.column String indicating the name of a column from the metadata table in which are stored the replicate IDs. This column is used only if \code{paired.test = TRUE}. Default: \code{NULL}.
#' @param linear.FC.th Number indicating the (absolute) fold change threshold (linear scale) to use to define differential proteins. Default: \code{2}.
#' @param linear.FC.unresp.range A numeric 2-elements vector indicating the range (linear scale) used to define the unresponsive fold changes. Default: \code{c(1/1.1, 1.1)}.
#' @param padj.th Numeric value indicating the p.adjusted threshold to apply to the differential analyses. Default: \code{0.05}.
#' @param padj.method String indicating the method to use to correct the p-values. One among: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: \code{BH}.
#' @param stat.test String indicating the type of statistic test to use. One among: "t-test" and "wilcoxon". Default: \code{"t.test"}.
#' @param paired.test Logical value indicating whether paired statistical test should be performed. For each contrast it will be verified that replicate IDs are not repeated within a group and, that replicate number and the identifiers coincide between the two groups. Default: \code{FALSE}.
#' @param up.color String indicating the color to use for up-regulated proteins in the plots. Default: \code{"indianred"}.
#' @param down.color String indicating the color to use for up-regulated proteins in the plots. Default: \code{"steelblue"}.
#' @param unresponsive.color String indicating the color to use for unresponsive proteins in the plots. Default: \code{"purple"}.
#' @param null.color String indicating the color to use for null proteins in the plots. Default: \code{"gray"}.
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'imp'. Default: \code{"imputed"}.
#' @param overwrite.analyses Logical value to indicate whether overwrite analyses already generated. Default: \code{FALSE}.
#'
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @importFrom RColorBrewer brewer.pal
#' @importFrom purrr pmap
#' @import ggtext
#' @importFrom grDevices colorRampPalette
#' @importFrom stats t.test wilcox.test p.adjust t.test wilcox.test
#'
#' @author Sebastian Gregoricchio
#'
#' @export diff.analyses

diff.analyses =
  function(DEprot.object,
           contrast.list,
           replicate.column = NULL,
           linear.FC.th = 2,
           linear.FC.unresp.range = c(1/1.1, 1.1),
           padj.th = 0.05,
           padj.method = "BH", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
           stat.test = "t.test",
           paired.test = FALSE,
           up.color = "indianred",
           down.color = "steelblue",
           unresponsive.color = "purple",
           null.color = "gray",
           which.data = "imputed",
           overwrite.analyses = FALSE) {

    # ### Packages
    # require(dplyr)
    # require(ggplot2)
    # require(patchwork)


    ### check object and extract metadata table
    if (!("DEprot" %in% class(DEprot.object))) {
      if (!("DEprot.analyses" %in% class(DEprot.object))) {
        stop("The input must be an object of class 'DEprot'.")
        #return(DEprot.object)
      }
    }

    meta.tb = DEprot.object@metadata


    ### Check 'replicate.column' presence and eventually the paired.test variable
    if (!is.null(replicate.column[[1]])) {
      if (!(replicate.column %in% colnames(meta.tb))) {
        stop("The 'replicate.column' is not present in the metadata of the object.")
        #return(DEprot.object)
      } else {
        meta.tb[,replicate.column[[1]]] = as.character(meta.tb[,replicate.column[[1]]])
      }
    } else if (paired.test == TRUE) {
      warning("The option 'paired-test' is TRUE, however no 'replicate.column' is provided: 'paired.test' will ignored.")
      paired.test = FALSE
    }



    ### Check contrasts
    if (all(is.na(DEprot.object@contrasts))) {
      if ("character" %in% class(contrast.list)) {
        contrasts = list(contrast.list)
      } else if ("list" %in% class(contrast.list)) {
        contrasts = contrast.list
      } else {
        stop("The 'contrast.list' must be a list of 3-elements vectors indicating: metadata_column, variable_1, variable_2.")
        #return(DEprot.object)
      }

    } else if (overwrite.analyses == FALSE) {
      DEprot.object@contrasts
      stop("The DEprot object contains already a contrast list.\n",
           "       To overwrite the contrast list set the parameter `overwrite.analyses = TRUE`")
      #return(DEprot.object)
    } else {
      if ("character" %in% class(contrast.list)) {
        contrasts = list(contrast.list)
      } else if ("list" %in% class(contrast.list)) {
        contrasts = contrast.list
      } else {
        stop("The 'contrast.list' must be a list of 3-elements vector indicating: metadata_column, variable_1, variable_2.")
        #return(DEprot.object)
      }
    }


    # check the presence of columns and variables in the contrast
    contrasts.info = list()

    for (i in 1:length(contrasts)) {
      if (!(contrasts[[i]][1] %in% colnames(meta.tb))) {
        stop(paste0("The column indicated in the contrast #", i, " ('",contrasts[[i]][1],"'), it is not available in the metadata table."))
        #return(DEprot.object)
      } else if (!(contrasts[[i]][2] %in% meta.tb[,contrasts[[i]][1]])) {
        stop(paste0("In the contrast #", i, " ('",contrasts[[i]][1],"'), the first variable ('",contrasts[[i]][2],"') is not available."))
        #return(DEprot.object)
      } else if (!(contrasts[[i]][3] %in% meta.tb[,contrasts[[i]][1]])) {
        stop(paste0("In the contrast #", i, " ('",contrasts[[i]][1],"'), the first variable ('",contrasts[[i]][3],"') is not available."))
        #return(DEprot.object)
      } else {
        # Determine whether it can be run in paired mode
        if (paired.test == TRUE) {
          # re order table by condition and replicate
          meta.tb =
            meta.tb %>%
            dplyr::arrange(.data[[contrasts[[i]][1]]],
                           .data[[replicate.column]])

          reps.group.1 = meta.tb[meta.tb[,contrasts[[i]][1]] == contrasts[[i]][2], replicate.column]
          reps.group.2 = meta.tb[meta.tb[,contrasts[[i]][1]] == contrasts[[i]][3], replicate.column]

          # check duplicated rep IDs
          if (length(reps.group.1) != length(unique(reps.group.1))) {
            stop("At least on replicate ID in the 'replicate.column' is duplicated.")
            #return(DEprot.object)
          }

          if (length(reps.group.2) != length(unique(reps.group.2))) {
            stop("At least on replicate ID in the 'replicate.column' is duplicated.")
            #return(DEprot.object)
          }

          # check whether the reps are the same among the two groups
          if (length(reps.group.1) == length(reps.group.2)) {
            if (all(reps.group.1 == reps.group.2)) {
              paired.test.contrast = TRUE
            } else {
              message(paste0("For contrast #", i, " replicate IDs are not equivalent among groups.\nPaired-test is coerced to `FALSE` for this contrast."))
              paired.test.contrast = FALSE
            }
          } else {
            message(paste0("For contrast #", i, " the number of replicates is different among groups.\nPaired-test is coerced to `FALSE` for this contrast."))
            paired.test.contrast = FALSE
          }
        } else {
          paired.test.contrast = FALSE
        }


        # Collect info
        contrasts.info[[i]] = list(metadata.column = contrasts[[i]][1],
                                   var.1 = contrasts[[i]][2],
                                   var.2 = contrasts[[i]][3],
                                   group.1 = meta.tb[meta.tb[,contrasts[[i]][1]] == contrasts[[i]][2],"column.id"],
                                   group.2 = meta.tb[meta.tb[,contrasts[[i]][1]] == contrasts[[i]][3],"column.id"],
                                   paired.test = paired.test.contrast)
        names(contrasts.info)[i] = paste0(contrasts[[i]][1], "_", contrasts[[i]][2], ".vs.", contrasts[[i]][3])
      }
    }



    ### Check and extract table
    if (tolower(which.data) == "raw") {
      if (!is.null(DEprot.object@raw.counts)) {
        mat = DEprot.object@raw.counts
        data.used = "raw"
      } else {
        stop(paste0("Use of RAW counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        #return(DEprot.object)
      }
    } else if (tolower(which.data) %in% c("norm", "normalized", "normal")) {
      if (!is.null(DEprot.object@norm.counts)) {
        mat = DEprot.object@norm.counts
        data.used = "normalized"
      } else {
        stop(paste0("Use of NORMALIZED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        #return(DEprot.object)
      }
    } else if (tolower(which.data) %in% c("imputed", "imp", "impute")) {
      if (!is.null(DEprot.object@imputed.counts)) {
        mat = DEprot.object@imputed.counts
        data.used = "imputed"
      } else {
        stop(paste0("Use of IMPUTED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        #return(DEprot.object)
      }
    } else {
      stop(paste0("The 'which.data' value is not recognized.\n",
                  "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
      #return(DEprot.object)
    }


    ## Convert table to log2
    if (!is.numeric(DEprot.object@log.base)) {
      message("The log.base is not numeric, linear counts are assumed. Counts matrix will be converted to log2+1 values to analyze the data.")
      mat.log2 = log2(mat + 1)
    } else if (as.numeric(DEprot.object@log.base) != 2) {
      message("The log.base is not 2, counts will be converted to log2 values to analyze the data.")
      mat.log2 = log2(DEprot.object@log.base^mat)
    } else {
      mat.log2 = mat
    }



    ############
    ### Compute differential analyses
    linear.FC.unresp.range = sort(linear.FC.unresp.range)
    diff.analyses.list = list()

    for (i in 1:length(contrasts.info)) {
      # means and FoldChange
      diff.tb =
        data.frame(prot.id = rownames(mat.log2),
                   basemean.log2 = rowMeans(mat.log2[,c(contrasts.info[[i]]$group.1, contrasts.info[[i]]$group.2)], na.rm = TRUE),
                   log2.mean.group1 = rowMeans(mat.log2[,contrasts.info[[i]]$group.1], na.rm = TRUE),
                   log2.mean.group2 = rowMeans(mat.log2[,contrasts.info[[i]]$group.2], na.rm = TRUE)) %>%
        dplyr::mutate(log2.Fold_group1.vs.group2 = log2.mean.group1 - log2.mean.group2)

      ## Wilcoxon/t.test pval
      # split the matrix in vectors
      pval.list = c()
      if (tolower(stat.test) %in% c("t.test", "ttest", "t-test", "t", "student")) {
        for (k in 1:nrow(mat.log2)){
          pval.list[k] = suppressWarnings(t.test(x = as.vector(mat.log2[k,contrasts.info[[i]]$group.1]),
                                                 y = as.vector(mat.log2[k,contrasts.info[[i]]$group.2]),
                                                 paired = contrasts.info[[i]]$paired.test,
                                                 exact = TRUE))$p.value
        }
      } else {
        for (k in 1:nrow(mat.log2)){
          pval.list[k] = suppressWarnings(wilcox.test(x = as.vector(mat.log2[k,contrasts.info[[i]]$group.1]),
                                                      y = as.vector(mat.log2[k,contrasts.info[[i]]$group.2]),
                                                      paired = contrasts.info[[i]]$paired.test,
                                                      exact = TRUE))$p.value
        }
      }


      diff.tb$p.value = pval.list
      diff.tb$padj = p.adjust(p = diff.tb$p.value, method = padj.method)


      ## Define signif status
      de.status =
        function(FC, p) {
          ifelse(p < padj.th,
                 yes = ifelse(abs(FC) >= log2(linear.FC.th),
                              yes = ifelse(sign(FC) == 1,
                                           yes = contrasts.info[[i]]$var.1,
                                           no = contrasts.info[[i]]$var.2),
                              no = "null"),
                 no = ifelse(FC >= log2(linear.FC.unresp.range[1]) & FC <= log2(linear.FC.unresp.range[2]),
                             yes = "unresponsive",
                             no = "null"))
        }

      diff.tb$diff.status =
        unlist(purrr::pmap(.l = list(FC = diff.tb$log2.Fold_group1.vs.group2,
                                     p = diff.tb$padj),
                           .f = function(FC,p){de.status(FC,p)}))

      diff.tb =
        diff.tb %>%
        dplyr::mutate(diff.status = factor(diff.status,
                                           levels = c(contrasts.info[[i]]$var.2,
                                                      contrasts.info[[i]]$var.1,
                                                      "unresponsive",
                                                      "null")))


      ## Run PCA
      PCA.data = DEprot::perform.PCA(DEprot.object = DEprot.object,
                                     sample.subset = c(contrasts.info[[i]]$group.1, contrasts.info[[i]]$group.2),
                                     which.data = which.data)

      scatter.PC1.2 = DEprot::plot.PC.scatter(DEprot.PCA.object = PCA.data, PC.x = 1, PC.y = 2, color.column = contrasts.info[[i]]$metadata.column, shape.column = switch(isTRUE(paired.test)+1, NULL, replicate.column), plot.zero.lines = FALSE) + geom_hline(yintercept = 0, colour = "gray", linetype = 2) + theme(legend.position = "none")
      if (length(unique(sign(PCA.data@PCs$PC1))) > 1){scatter.PC1.2 + geom_vline(xintercept = 0, colour = "gray", linetype = 2)}
      scatter.PC2.3 = DEprot::plot.PC.scatter(DEprot.PCA.object = PCA.data, PC.x = 3, PC.y = 2, color.column = contrasts.info[[i]]$metadata.column, shape.column = switch(isTRUE(paired.test)+1, NULL, replicate.column))



      ## Run Correlations
      corr.data.pearson = DEprot::plot.correlation.heatmap(DEprot.object = DEprot.object,
                                                           correlation.method = "pearson",
                                                           sample.subset = c(contrasts.info[[i]]$group.1, contrasts.info[[i]]$group.2),
                                                           which.data = which.data,
                                                           correlation.scale.limits = c(NA,1))

      corr.data.spearman = DEprot::plot.correlation.heatmap(DEprot.object = DEprot.object,
                                                            correlation.method = "spearman",
                                                            sample.subset = c(contrasts.info[[i]]$group.1, contrasts.info[[i]]$group.2),
                                                            which.data = which.data,
                                                            correlation.scale.limits = c(NA,1))



      ## Define params and counts
      colors.plots = c(up.color, down.color, null.color, unresponsive.color)
      names(colors.plots) = c(contrasts.info[[i]]$var.1, contrasts.info[[i]]$var.2, "null", "unresponsive")


      ### Summarize the number of diff proteins
      n.diff =
        data.frame(diff.tb %>%
                     dplyr::group_by(diff.status, .drop = F) %>%
                     dplyr::summarise(n = n(),
                                      median.FoldChange = median(log2.Fold_group1.vs.group2)))


      ## Make volcano plot
      volcano =
        ggplot(diff.tb,
               aes(x = log2.Fold_group1.vs.group2,
                   y = -log10(padj),
                   color = diff.status)) +
        geom_point(stroke = NA,
                   alpha = 0.5,
                   size = 2) +
        scale_color_manual(values = colors.plots,
                           name = "Differential\nstatus",
                           drop = FALSE) +
        geom_hline(yintercept = -log10(padj.th), linetype = 2, color = "gray40") +
        geom_vline(xintercept = c(-1,1)*log2(linear.FC.th), linetype = 2, color = "gray40") +
        ylab("-log~10~(*P~adj~*)") +
        #xlab(paste0("log~2~(Fold Change<sub>",contrasts.info[[i]]$var.1, "/", contrasts.info[[i]]$var.2,"</sub>)")) +
        xlab(paste0("log~2~(Fold Change<sub>",contrasts.info[[i]]$var.1,"</sup>&frasl;<sub>",contrasts.info[[i]]$var.2,"</sub></sub>)")) +
        ggtitle(paste0("**",contrasts.info[[i]]$var.1, "** *vs* **", contrasts.info[[i]]$var.2, "**")) +
        theme_classic() +
        theme(axis.text = ggtext::element_markdown(color = "black"),
              axis.title = ggtext::element_markdown(color = "black"),
              axis.ticks = element_line(color = "black"),
              plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
              aspect.ratio = 1)

      x.max = max(abs(ggplot_build(volcano)$layout$panel_params[[1]]$x.range))

      volcano =
        volcano +
        xlim(c(-1,1)*x.max) +
        annotate(geom = "text",
                 x = -Inf, y = +Inf,
                 color = colors.plots[2],
                 hjust = -0.2, vjust = 1.5,
                 label = paste0("n = ",n.diff[n.diff$diff.status == contrasts.info[[i]]$var.2,"n"])) +
        annotate(geom = "text",
                 x = +Inf, y = +Inf,
                 hjust = 1.2, vjust = 1.5,
                 color = colors.plots[1],
                 label = paste0("n = ",n.diff[n.diff$diff.status == contrasts.info[[i]]$var.1,"n"]))



      ## Make MA-plot
      ma.plot =
        ggplot() +
        stat_density_2d(data = diff.tb,
                        aes(x = basemean.log2,
                            y = log2.Fold_group1.vs.group2,
                            fill = after_stat(count)),
                        geom = "raster",
                        contour = FALSE,
                        show.legend = TRUE,
                        n = 200,
                        adjust = 5) +
        scale_fill_gradientn(colours = colorRampPalette(colors = RColorBrewer::brewer.pal(9, "Blues"))(101),
                             name = "Count")


      if (sum((n.diff %>% dplyr::filter(!(diff.status %in% c("unresponsive", "null"))))$n, na.rm = T) > 0) {
        ma.plot =
          ma.plot +
          geom_point(data = dplyr::filter(diff.tb, diff.status %in% c(contrasts.info[[i]]$var.1,contrasts.info[[i]]$var.2)),
                     mapping = aes(x = basemean.log2,
                                   y = log2.Fold_group1.vs.group2,
                                   color = diff.status),
                     alpha = 0.5,
                     size = 2,
                     stroke = NA,
                     show.legend = TRUE,
                     inherit.aes = FALSE) +
          scale_color_manual(values = colors.plots, name = "Differential\nstatus", drop = FALSE)
      }

      ma.plot =
        ma.plot +
        geom_hline(yintercept = c(-1,1)*log2(linear.FC.th), linetype = 2, color = "gray40") +
        geom_hline(yintercept = 0, linetype = 1, color = "steelblue") +
        theme_classic() +
        xlab("log~2~(Base Mean)") +
        ylab(paste0("log~2~(Fold Change<sub>",contrasts.info[[i]]$var.1,"</sup>&frasl;<sub>",contrasts.info[[i]]$var.2,"</sub></sub>)")) +
        ggtitle(paste0("**",contrasts.info[[i]]$var.1, "** *vs* **", contrasts.info[[i]]$var.2, "**")) +
        scale_x_continuous(expand = c(0,0)) +
        annotate(geom = "text",
                 x = -Inf, y = -Inf,
                 color = colors.plots[2],
                 hjust = -0.2, vjust = -0.5,
                 label = paste0("n = ",n.diff[n.diff$diff.status == contrasts.info[[i]]$var.2,"n"])) +
        annotate(geom = "text",
                 x = -Inf, y = +Inf,
                 hjust = -0.2, vjust = 1.5,
                 color = colors.plots[1],
                 label = paste0("n = ",n.diff[n.diff$diff.status == contrasts.info[[i]]$var.1,"n"])) +
        theme(axis.ticks = element_line(color = "black"),
              axis.text = element_text(color = "black"),
              axis.title = ggtext::element_markdown(color = "black"),
              plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
              aspect.ratio = 0.6)

      y.max = max(abs(ggplot_build(ma.plot)$layout$panel_params[[1]]$y.range))
      ma.plot = ma.plot + scale_y_continuous(expand = c(0,0), limits = c(-1,1)*y.max)


      ### rename columns
      colnames(diff.tb) = gsub("group1", contrasts.info[[i]]$var.1, colnames(diff.tb))
      colnames(diff.tb) = gsub("group2", contrasts.info[[i]]$var.2, colnames(diff.tb))
      colnames(n.diff) = gsub("group1", contrasts.info[[i]]$var.1, colnames(n.diff))
      colnames(n.diff) = gsub("group2", contrasts.info[[i]]$var.2, colnames(n.diff))


      ### Build list
      diff.analyses.list[[i]] = list(results = diff.tb,
                                     n.diff =  n.diff,
                                     PCA.data = PCA.data,
                                     PCA.plots = (scatter.PC1.2 | scatter.PC2.3) / PCA.data@cumulative.PC.plot,
                                     correlations = corr.data.pearson@heatmap | corr.data.spearman@heatmap,
                                     volcano = volcano,
                                     MA.plot = ma.plot)
    }

    names(diff.analyses.list) = names(contrasts.info)


    ### Make DEprot.analyses object
    DEprot.object.analyses =
      new(Class = "DEprot.analyses",
          metadata = DEprot.object@metadata,
          raw.counts = DEprot.object@raw.counts,
          norm.counts =  DEprot.object@norm.counts,
          imputed.counts = DEprot.object@imputed.counts,
          log.base = DEprot.object@log.base,
          log.transformed = DEprot.object@log.transformed,
          imputed = DEprot.object@imputed,
          imputation = DEprot.object@imputation,
          normalized = DEprot.object@normalized,
          normalization.method = DEprot.object@normalization.method,
          boxplot.raw = DEprot.object@boxplot.raw,
          boxplot.norm = DEprot.object@boxplot.norm,
          boxplot.imputed = DEprot.object@boxplot.imputed,
          analyses.result.list = diff.analyses.list,
          contrasts = contrasts.info,
          differential.analyses.params = list(linear.FC.th = linear.FC.th,
                                              linear.FC.unresp.range = linear.FC.unresp.range,
                                              padj.th = padj.th,
                                              padj.method = padj.method,
                                              counts.used = data.used))

    return(DEprot.object.analyses)
  } # END function
