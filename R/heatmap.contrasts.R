#' @title heatmap.contrasts
#'
#' @description Plots an heatmap of the log2(FoldChange) for differentially expressed proteins
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrasts Numeric vector indicating the position of the contrast to use for the plotting. Default: \code{NULL} (all contrasts).
#' @param top.n Numeric value indicated the top differentially expressed proteins to consider for each contrast selected. The rank is based on the product of log2Fc and -log10Padj. Default: \code{NULL} (all differential proteins).
#' @param distance.method The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given. Default: \code{"euclidean"}.
#' @param clustering.method The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (UPGMA), "mcquitty" (WPGMA), "median" (WPGMC) or "centroid" (UPGMC). Default: \code{"complete"}.
#' @param high.color String indicating the color to use for up-regulated protein FoldChanges in the plot. Default: \code{"firebrick"}.
#' @param low.color String indicating the color to use for down-regulated protein FoldChanges in the plot. Default: \code{"\#2166AC"} (blue).
#' @param mid.color String indicating the color to use for the 0 in the plots. Default: \code{"white"}.
#' @param na.color String indicating the color to use for the NA values in the heatmap. Default: \code{"gray"}.
#' @param color.limits Numeric vector of length 2 indicating lower and upper limit of the color scale values. Default: \code{c(NA,NA)}, automatic limits applied.
#' @param cell.border.color String indicating the color to use for the individual cells border. Default: \code{NA} (no border).
#' @param cell.border.width Numeric value indicating the width of the cell border line. Ignored when \code{cell.border.color = NA}. Default: \code{0.5}.
#' @param title String indicating the title to use, markdown formatting supported. Default: \code{NULL} (automatic title).
#' @param show.protein.names Logical value to indicate whether the protein names should be displayed. Default: \code{FALSE}.
#' @param protein.names.pattern Character indicating a regular expression to remove from the protein IDs. Default: \code{NULL}, no alterations in the protein IDs.
#' @param use.uncorrected.pvalue Logical value indicating whether it should be used the normal p-value instead of the adjusted one (differential proteins numbers are recomputed). Default: \code{FALSE}, padj is used.
#'
#' @return A \code{DEprot.contrast.heatmap} object, which contains the heatmap (ggplot) and the hclust object used to order the rows.
#'
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#' @import legendry
#' @import ggdendro
#' @importFrom purrr pmap
#' @importFrom reshape2 dcast
#' @importFrom stats dist hclust
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' heatmap.contrasts(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, top.n = 5)
#'
#'
#' @export heatmap.contrasts

heatmap.contrasts =
  function(DEprot.analyses.object,
           contrasts = NULL,
           top.n = NULL,
           distance.method = "euclidean",
           clustering.method = "complete",
           high.color = "firebrick",
           low.color = "#2166AC",
           mid.color = "white",
           na.color = "gray",
           color.limits = c(NA,NA),
           cell.border.color = NA,
           cell.border.width = 0.5,
           title = NULL,
           show.protein.names = FALSE,
           protein.names.pattern = NULL,
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

    ######################################################################################

    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      stop("The input must be an object of class 'DEprot.analyses'.")
    }

    ### check and collect contrast
    if (!is.null(contrasts)) {
      if (!all(contrasts %in% 1:length(DEprot.analyses.object@analyses.result.list))) {
        stop("At least one of the 'contrasts' indicated is not available in the results list.")
      } else {
        contr.list = contrasts
      }
    } else {
      contr.list = 1:length(DEprot.analyses.object@analyses.result.list)
    }




    ### Collect the data required for plotting
    fc.tb.list = list()
    differential_proteins = c()

    for (i in 1:length(contr.list)) {
      fc.tb.list[[i]] = DEprot.analyses.object@analyses.result.list[[contr.list[i]]]$results
      colnames(fc.tb.list[[i]])[5] = "log2.Fold_group1.vs.group2"

      # Change the column name for 'FDR' column into 'padj' for prolfqua analyses
      if ("strategy" %in% names(DEprot.analyses.object@differential.analyses.params)) {
        fc.tb.list[[i]] = fc.tb.list[[i]] %>% dplyr::rename(padj = FDR)
      }


      ## redefine the diff genes
      if (use.uncorrected.pvalue == TRUE) {
        fc.tb.list[[i]]$padj = fc.tb.list[[i]]$p.value

        fc.tb.list[[i]]$diff.status =
          unlist(purrr::pmap(.l = list(FC = fc.tb.list[[i]]$log2.Fold_group1.vs.group2,
                                       p = fc.tb.list[[i]]$padj),
                             .f = function(FC,p){
                               de.status(FC, p,
                                         contrasts.info = DEprot.analyses.object@contrasts[[contr.list[i]]],
                                         params = DEprot.analyses.object@differential.analyses.params)}))
      }

      ## Use log2FC * -log10(padj) as ranking score
      fc.tb.list[[i]] =
        fc.tb.list[[i]] %>%
        dplyr::mutate(contrast = names(DEprot.analyses.object@contrasts)[contr.list[i]],
                      ranking.score = log2.Fold_group1.vs.group2 * -log10(fc.tb.list[[i]]$padj)) %>%
        dplyr::arrange(desc(abs(ranking.score))) %>%
        dplyr::mutate(rank = 1:nrow(fc.tb.list[[i]]))


      ## Define differential proteins for this contrast
      diff.prot = (fc.tb.list[[i]] %>% dplyr::filter(!(diff.status %in% c("null", "unresponsive"))))$prot.id

      if (is.null(top.n)) {
        differential_proteins = unique(c(differential_proteins, diff.prot))
      } else {
        differential_proteins = unique(c(differential_proteins, diff.prot[1:min(c(length(diff.prot), top.n))]))
      }
    } # end FC tables



    ### Combine foldchange tables
    combined.fc =
      do.call(rbind,
              lapply(fc.tb.list,
                     function(x)(
                       x %>%
                         dplyr::filter(prot.id %in% differential_proteins) %>%
                         dplyr::select(prot.id, basemean.log2,
                                       log2.Fold_group1.vs.group2, p.value, padj, diff.status,
                                       contrast, ranking.score,rank)
                     ))) %>%
      dplyr::mutate(contrast = factor(contrast, levels = names(DEprot.analyses.object@contrasts)[contr.list]))

    if (nrow(combined.fc) == 0) {stop("No data to be shown.")}


    ### remove protein pattern if required
    if (!is.null(protein.names.pattern)) {
      if ("character" %in% class(protein.names.pattern)) {
        combined.fc = combined.fc %>% dplyr::mutate(prot.id = gsub(protein.names.pattern, "", prot.id))
      } else {
        stop("The 'protein.names.pattern' must be a character indicating a regular expression to remove from the protein IDs.")
      }
    }



    ### Make heatmap
    heatmap =
      ggplot(data = combined.fc,
             aes(x = contrast,
                 y = prot.id,
                 fill = log2.Fold_group1.vs.group2)) +
      geom_tile(color = cell.border.color,
                linewidth = cell.border.width,
                show.legend = T) +
      scale_x_discrete(name = "Contrast ID", expand = c(0,0), position = "top") +
      scale_fill_gradient2(name = "log~2~(Fold Change)",
                           low = low.color,
                           mid = mid.color,
                           high = high.color,
                           midpoint = 0,
                           na.value = na.color,
                           limits = color.limits) +
      ggtitle(title) +
      theme_classic() +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = ggtext::element_markdown(color = "black", angle = 30, hjust = 0),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
            legend.title = ggtext::element_markdown(color = "black"),
            axis.title = element_blank())

    if (length(unique(combined.fc$prot.id)) > 1) {
      ### Define row cluster
      combined.matrix =
        reshape2::dcast(data = combined.fc[,c("prot.id","contrast","log2.Fold_group1.vs.group2")],
                        formula = prot.id ~ contrast,
                        value.var = "log2.Fold_group1.vs.group2")

      rownames(combined.matrix) = combined.matrix$prot.id

      row.clust = hclust(d = stats::dist(x = as.matrix(dplyr::select(.data = combined.matrix, -prot.id)),
                                         method = distance.method),
                         method = clustering.method)
      if (!is.null(protein.names.pattern)) {row.clust$labels = gsub(protein.names.pattern, "", row.clust$labels)}

      row.clust$call = "hclust(d = dist(x = foldchange.matrix, method = distance.method), method = clustering.method)"
      heatmap = heatmap + legendry::scale_y_dendro(clust = row.clust, name = "Protein ID", expand = c(0,0))
    } else {
      heatmap = heatmap + scale_y_discrete(name = "Protein ID", expand = c(0,0))
      row.clust = NULL
    }



    if (show.protein.names == TRUE) {
      heatmap = heatmap + theme(axis.text.y = ggtext::element_markdown(color = "black"))
    } else {
      heatmap = heatmap + theme(axis.text.y = element_blank())
    }



    ### Export heatmap.contrast object
    DEprot.contrast.heatmap.object =
      new(Class = "DEprot.contrast.heatmap",
          heatmap = heatmap,
          cluster = row.clust)

    return(DEprot.contrast.heatmap.object)
  } # END of function

