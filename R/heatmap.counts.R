#' @title heatmap.counts
#'
#' @description Plots an heatmap of the counts (raw, normalized, or imputed). It is possible to perform a scaling (z-score) by row or by column.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'imp'. Default: \code{"imputed"}.
#' @param contrast Numeric vector indicating the position of the contrast to use for the plotting. Only differential proteins in this contrast will be shown. Option available only for an object of class \code{DEprot.analyses}. Default: \code{NULL} (non differential protein selection).
#' @param top.n Numeric value indicated the top differentially expressed proteins to consider for the contrast selected. The rank is based on the product of log2Fc and -log10Padj. Option available only for an object of class \code{DEprot.analyses}. Default: \code{NULL} (all differential proteins of that contrast).
#' @param sample.subset Character vector indicating a subset of samples to display. The identifiers must correspond to a IDs in the \code{column.id} column of the object's metadata. Default: \code{NULL} (all samples are shown).
#' @param protein.subset Character vector indicating a subset of proteins to display. The identifiers must correspond to the full row.names of the counts table (equivalent to the \code{prot.id} column of the fold change table of \code{DEprot.analyses} object). This options is can be used in combination with \code{contrast} and \code{top.n}. Default: \code{NULL} (all proteins are shown).
#' @param group.by.metadata.column String indicating a column from the metadata table. This column will be used to define sample groups, and for each group it will be computed a mean of the counts. Default: \code{NULL} (no groups).
#' @param scale String indicating whether Z-scores should be computed. Possible choices: "row" or "column". Default: \code{NULL} (no scaling).
#' @param clust.rows Logical value indicating whether heatmap rows (proteins) should be clustered. Default: \code{TRUE}.
#' @param clust.columns Logical value indicating whether heatmap columns (samples or groups) should be clustered. Default: \code{TRUE}.
#' @param distance.method The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given. Default: \code{"euclidean"}.
#' @param clustering.method The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (UPGMA), "mcquitty" (WPGMA), "median" (WPGMC) or "centroid" (UPGMC). Default: \code{"complete"}.
#' @param palette List of colors to use for the color gradient of the heatmap. This parameters is used when raw-counts are plotted. For scaled data see the high, low and mid parameters. Default: \code{RColorBrewer\:\:brewer.pal(n = 9, name = "Blues")}.
#' @param high.color String indicating the color to use for positive Z-score for scaled data. Default: \code{"indianred"}.
#' @param low.color String indicating the color to use for negative Z-score for scaled data. Default: \code{"\#2166AC"} (blue).
#' @param mid.color String indicating the color to use for the 0 Z-score value for scaled data. Default: \code{"white"}.
#' @param na.color String indicating the color to use for the NA values in the heatmap. Default: \code{"gray"}.
#' @param color.limits Numeric vector of length 2 indicating lower and upper limit of the color scale values. Default: \code{c(NA,NA)}, automatic limits applied.
#' @param cell.border.color String indicating the color to use for the individual cells border. Default: \code{NA} (no border).
#' @param cell.border.width Numeric value indicating the width of the cell border line. Ignored when \code{cell.border.color = NA}. Default: \code{0.5}.
#' @param show.protein.names Logical value to indicate whether the protein names should be displayed. Default: \code{FALSE}.
#' @param protein.names.pattern Character indicating a regular expression to remove from the protein IDs. Default: \code{NULL}, no alterations in the protein IDs.
#' @param title String indicating the title to use, markdown formatting supported. Default: \code{NULL} (automatic title).
#' @param use.uncorrected.pvalue Logical value indicating whether it should be used the normal p-value instead of the adjusted one (differential proteins numbers are recomputed). Default: \code{FALSE}, padj is used.
#'
#' @return A \code{DEprot.contrast.heatmap} object, which contains the heatmap and the hclust object used to order the rows.
#'
#' @import dplyr
#' @import ggplot2
#' @import legendry
#' @import ggdendro
#' @importFrom RColorBrewer brewer.pal
#' @importFrom purrr pmap
#' @importFrom reshape2 melt
#' @import ggtext
#' @importFrom stats dist hclust
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' require(legendry)
#'
#' # Counts per each sample
#' heatmap.counts(DEprot.object = DEprot::test.toolbox$diff.exp.limma,
#'                top.n = 5,
#'                contrast = 1,
#'                scale = "column")
#'
#' # Counts averaged by group
#' heatmap.counts(DEprot.object = DEprot::test.toolbox$diff.exp.limma,
#'                top.n = 5,
#'                contrast = 1,
#'                scale = "column",
#'                group.by.metadata.column = "combined.id")
#'
#'
#' @export heatmap.counts

heatmap.counts =
  function(DEprot.object,
           which.data = "imputed",
           contrast = NULL,
           top.n = NULL,
           sample.subset = NULL,
           protein.subset = NULL,
           group.by.metadata.column = NULL,
           scale = NULL,
           clust.rows = TRUE,
           clust.columns = TRUE,
           distance.method = "euclidean",
           clustering.method = "complete",
           palette = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
           high.color = "firebrick",
           low.color = "#2166AC",
           mid.color = "white",
           na.color = "gray",
           color.limits = c(NA,NA),
           cell.border.color = NA,
           cell.border.width = 0.5,
           show.protein.names = FALSE,
           protein.names.pattern = NULL,
           title = NULL,
           use.uncorrected.pvalue = FALSE) {

    # ### Libraries
    # require(dplyr)
    # require(ggplot2)
    # require(legendry)

    ### Internal functions
    de.status =
      function(FC, p, contrasts.info, params) {
        ifelse(p < params$padj.th,
               yes = ifelse(abs(FC) >= log2(params$linear.FC.th),
                            yes = ifelse(sign(FC) == 1,
                                         yes = contrasts.info$var.1,
                                         no = contrasts.info$var.2),
                            no = "null"),
               no = ifelse(FC >= log2(params$linear.FC.unresp.range[1]) & FC <= log2(params$linear.FC.unresp.range[2]),
                           yes = "unresponsive",
                           no = "null"))
      }


    check.matrix =
      function(m){
        warn = "Upon subsetting, no values to show are left."
        if (!is.logical(m)) {
          if (nrow(m) == 0 | ncol(m) == 0) {
            stop(warn)
          }
        } else {
          stop(warn)
        }
      }



    is.nan_df = function(data.frame) {do.call(cbind, lapply(data.frame, is.nan))}



    scale_matrix =
      function(mat, scale){
        if(!(scale %in% c("none", "row", "column"))){
          stop("scale argument shoud take values: 'none', 'row' or 'column'")
        }

        scale_rows = function(x){
          m = apply(x, 1, mean, na.rm = T)
          s = apply(x, 1, sd, na.rm = T)
          return((x - m) / s)
        }

        mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
        return(mat)
      }

    ######################################################################################

    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.object))) {
      if (!("DEprot" %in% class(DEprot.object))) {
        stop("The input must be an object of class 'DEprot' or 'DEprot.analyses'.")
      }
    }


    ### check grouping column
    if (!is.null(group.by.metadata.column)) {
      if (!(group.by.metadata.column %in% colnames(DEprot.object@metadata))) {
        stop(paste0("The 'group.by.metadata.column' is not present in the metadata of the object provided.\n",
                    "Available column IDs: ", paste0(colnames(DEprot.object@metadata), collapse = ", ")))
      } else {
        meta = DEprot.object@metadata
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
      return(DEprot.object)
    }



    ### Filter table of counts
    if (!is.null(sample.subset)) {
      mat.filtered = mat[,which(colnames(mat) %in% sample.subset), drop=FALSE]
    } else {
      mat.filtered = mat
    }
    check.matrix(mat.filtered)


    if (!is.null(protein.subset)) {
      mat.filtered = mat.filtered[which(rownames(mat.filtered) %in% protein.subset),,drop=FALSE]
    } else {
      mat.filtered = mat.filtered
    }
    check.matrix(mat.filtered)





    ### check and collect contrast (if possible: requires DEprot.analyses-class object)
    if (!is.null(contrast)) {
      if ("DEprot.analyses" %in% class(DEprot.object)) {
        if (is.numeric(contrast)) {
          if (contrast <= length(DEprot.object@analyses.result.list)) {
            fc.data = DEprot.object@analyses.result.list[[contrast]]$results
            colnames(fc.data)[5] = "log2.Fold_group1.vs.group2"

            # Change the column name for 'FDR' column into 'padj' for prolfqua analyses
            if ("strategy" %in% names(DEprot.object@differential.analyses.params)) {
              fc.data = fc.data %>% dplyr::rename(padj = FDR)
            }

            ## recompute up-downs if required
            if (use.uncorrected.pvalue == TRUE) {
              fc.data$padj = fc.data$p.value

              fc.data$diff.status =
                unlist(purrr::pmap(.l = list(FC = fc.data$log2.Fold_group1.vs.group2,
                                             p = fc.data$padj),
                                   .f = function(FC,p){
                                     de.status(FC, p,
                                               contrasts.info = DEprot.object@contrasts[[contrast]],
                                               params = DEprot.object@differential.analyses.params)}))
            }


           ## select top.n proteins
            if (!is.null(top.n)) {
              ## Use log2FC * -log10(padj) as ranking score
              fc.data =
                fc.data %>%
                dplyr::mutate(contrast = names(DEprot.object@contrasts)[contrast],
                              ranking.score = log2.Fold_group1.vs.group2 * -log10(fc.data$padj)) %>%
                dplyr::arrange(desc(abs(ranking.score))) %>%
                dplyr::mutate(rank = 1:nrow(fc.data))

              diff.prot = (fc.data %>% dplyr::filter(!(diff.status %in% c("null", "unresponsive"))))$prot.id
              diff.prot = unique(diff.prot[1:min(c(length(diff.prot), top.n))])

            } else {
              diff.prot = unique((fc.data %>% dplyr::filter(!(diff.status %in% c("null", "unresponsive"))))$prot.id)
            }


            ## filter the matrix for these proteins
            mat.filtered = mat.filtered[which(rownames(mat.filtered) %in% diff.prot),,drop=FALSE]
            check.matrix(mat.filtered)


          } else {
            stop("The 'contrast' indicated is not available.")
          }
        } else {
          stop("The 'contrast' must be a numeric value.")
        }
      } else {
        stop("The 'contrast' cannot be determined from a 'DEprot' object, provide a 'DEprot.analyses' object instead.")
      }
    }



    ## Average/group data if required
    if (!is.null(group.by.metadata.column)) {
       if (!is.null(sample.subset)) {
         meta = meta %>% dplyr::filter(column.id %in% sample.subset)
       }
      groups = unique(unique(meta[,group.by.metadata.column]))

      groups.tb = meta[,c("column.id", group.by.metadata.column)]
      colnames(groups.tb)[2] = "group"

      new.data = list()
      for (i in 1:length(groups)) {
        samples.in.group = (dplyr::filter(.data = groups.tb, group == groups[i]))$column.id
        new.data[[i]] = rowMeans(mat.filtered[,samples.in.group, drop = FALSE], na.rm = TRUE)
      }

      new.data = data.frame(new.data)
      colnames(new.data) = make.names(groups)
      rownames(new.data) = rownames(mat.filtered)
      mat.filtered = as.matrix(new.data)
    }




    ## Compute z-score if required
    if (!is.null(scale)) {
      if (tolower(scale) %in% c("col", "columns", "col", "column", "c", "cols")) {
        if (nrow(mat.filtered) > 1) {
          final.mat = scale_matrix(mat = (DEprot.object@log.base^mat.filtered)-1, scale = "column")
          final.mat[is.nan_df(final.mat)] = 0
          scale.legend.title = "**Z-score**<br>(by column)"
        } else {
          warning("Scaling by column is not possible because the number of rows is only 1.")
          final.mat = mat.filtered # no scaling applied
          scale.legend.title = paste(which.data, "counts")
        }
      } else if (tolower(scale) %in% c("row", "rows", "r")) {
        if (ncol(mat.filtered) > 1) {
          final.mat = scale_matrix(mat = (DEprot.object@log.base^mat.filtered)-1, scale = "row")
          final.mat[is.nan_df(final.mat)] = 0
          scale.legend.title = paste0("**Z-score**<br>", which.data, " counts<br>(by row)")
        } else {
          warning("Scaling by row is not possible because the number of columns is only 1.")
          final.mat = mat.filtered # no scaling applied
          scale.legend.title = paste(which.data, "counts")
        }
      } else {
        final.mat = mat.filtered # no scaling applied
        scale.legend.title = paste(which.data, "counts")
      }
    } else {
      final.mat = mat.filtered # no scaling applied
      scale.legend.title = paste(which.data, "counts")
    }



    ###############
    ## reshape matrix for plotting
    plotting.matrix =
      data.frame(final.mat) %>%
      dplyr::mutate(prot.id = rownames(final.mat)) %>%
      reshape2::melt(value.name = "score", id.vars = "prot.id", variable.name = "sample")


    if (!is.null(protein.subset)) {
      plotting.matrix = dplyr::mutate(.data = plotting.matrix, prot.id = factor(prot.id, levels = rev(unique(protein.subset))))
    }


    ### remove protein pattern if required
    if (!is.null(protein.names.pattern)) {
      if ("character" %in% class(protein.names.pattern)) {
        plotting.matrix = plotting.matrix %>% dplyr::mutate(prot.id = gsub(protein.names.pattern, "", prot.id))
      } else {
        stop("The 'protein.names.pattern' must be a character indicating a regular expression to remove from the protein IDs.")
      }
    }



    ## generate basic plot
    heatmap =
      ggplot(data = plotting.matrix,
             aes(x = sample,
                 y = prot.id,
                 fill = score)) +
      geom_tile(color = cell.border.color,
                linewidth = cell.border.width,
                show.legend = TRUE) +
      ggtitle(title) +
      theme_classic() +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
            legend.title = ggtext::element_markdown(color = "black"),
            axis.title = element_blank())


    ## Perform clustering
    if (clust.rows == TRUE & nrow(final.mat) > 1) {
      row.clust = hclust(d = stats::dist(x = final.mat,
                                         method = distance.method),
                         method = clustering.method)
      row.clust$call = "hclust(d = dist(x = counts.matrix, method = distance.method), method = clustering.method)"
      if (!is.null(protein.names.pattern)) {row.clust$labels = gsub(protein.names.pattern, "", row.clust$labels)}
      heatmap = heatmap + legendry::scale_y_dendro(clust = row.clust, name = "Protein ID", expand = c(0,0))
    } else {
      row.clust = NULL
      heatmap = heatmap + scale_y_discrete(name = "Protein ID", expand = c(0,0))
    }


    if (clust.columns == TRUE & ncol(final.mat) > 1) {
      columns.clust = hclust(d = stats::dist(x = t(final.mat),
                                             method = distance.method),
                         method = clustering.method)
      columns.clust$call = "hclust(d = dist(x = t(counts.matrix), method = distance.method), method = clustering.method)"
      heatmap =
        heatmap +
        legendry::scale_x_dendro(clust = columns.clust, name = "Sample ID", expand = c(0,0), position = "top") +
        theme(axis.text.x = ggtext::element_markdown(color = "black", angle = 90, hjust = 0))
    } else {
      columns.clust = NULL
      heatmap =
        heatmap +
        scale_x_discrete(name = "Sample ID", expand = c(0,0), position = "top") +
        theme(axis.text.x = ggtext::element_markdown(color = "black", angle = 30, hjust = 0))
    }


    ## Chose whether display proteins
    if (show.protein.names == TRUE) {
      heatmap = heatmap + theme(axis.text.y = ggtext::element_markdown(color = "black"))
    } else {
      heatmap = heatmap + theme(axis.text.y = element_blank())
    }



    ## re-set the color scale
    if (!is.null(scale)) {
      if (tolower(scale) %in% c("col", "columns", "col", "column", "c", "cols", "row", "rows", "r")) {
        heatmap =
          heatmap +
          scale_fill_gradient2(name = scale.legend.title,
                               low = low.color,
                               mid = mid.color,
                               high = high.color,
                               midpoint = 0,
                               na.value = na.color,
                               limits = color.limits)
      }
    } else {
      heatmap = heatmap + scale_fill_gradientn(name = scale.legend.title, colors = palette, na.value = na.color, limits = color.limits)
    }


    ##############################################
    ### Export heatmap object
    DEprot.counts.heatmap.object =
      new(Class = "DEprot.counts.heatmap",
          heatmap = heatmap,
          row.cluster = row.clust,
          column.cluster = columns.clust)

    return(DEprot.counts.heatmap.object)
  } # END of function

