#' @title compare.ranking
#'
#' @description Compare the two possible methods for the gene ranking used to perform Gene Set Enrichment Analyses (GSEA).
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast Number indicating the position of the contrast to use for the plotting.
#' @param color.up String indicating any R-supported color that will be used for dots with a positive fold change. Default: \code{"indianred"}.
#' @param color.down String indicating any R-supported color that will be used for dots with a negative fold change. Default: \code{"steelblue"}.
#' @param regression.line.color String indicating any R-supported color to use for the regression line and error. Default: \code{"purple"}.
#'
#' @return Two ggplots combined in an object of class \code{patchwork}.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom ggpubr stat_cor
#' @importFrom patchwork wrap_plots
#'
#' @export compare.ranking


compare.ranking =
  function(DEprot.analyses.object,
           contrast,
           color.up = "indianred",
           color.down = "steelblue",
           regression.line.color = "purple") {

    # ### Libraries
    # require(dplyr)
    # require(ggplot2)


    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      stop("The input must be an object of class 'DEprot.analyses'.")
      #return(invisible())
    }

    ### check and collect contrast
    if (is.numeric(contrast)) {
      if (contrast <= length(DEprot.analyses.object@analyses.result.list)) {
        data = DEprot.analyses.object@analyses.result.list[[contrast]]$results
        contrasts.info = DEprot.analyses.object@contrasts[[contrast]]
      } else {
        stop("The 'contrast' indicated is not available.")
        #return(invisible())
      }
    } else {
      stop("The 'contrast' must be a numeric value.")
      #return(invisible())
    }



    ### rank genes: foldchange mode
    gene_list_fc = data[,5]
    names(gene_list_fc) = data$prot.id
    gene_list_fc = sort(gene_list_fc, decreasing = TRUE)



    ### rank genes: correlation mode
    # extract tables by group
    counts_var1 = DEprot.analyses.object@imputed.counts[,contrasts.info$group.1]
    counts_var2 = DEprot.analyses.object@imputed.counts[,contrasts.info$group.2]
    group_idx = c(rep(1, ncol(counts_var2)), rep(2, ncol(counts_var1)))


    # run correlation
    corr_scores = sapply(1:nrow(counts_var1), function(x){suppressWarnings(cor.test(x = group_idx, y = c(counts_var2[x,],counts_var1[x,]), method = "spearman"))$estimate}, USE.NAMES = F)
    names(corr_scores) = rownames(counts_var1)
    gene_list_cor = sort(corr_scores, decreasing = TRUE)



    ### combine tables
    merged_tb =
      merge(data.frame(gene_list_fc),
            data.frame(gene_list_cor),
            by = "row.names",
            all = TRUE)

    merged_tb =
      merged_tb %>%
      dplyr::arrange(desc(gene_list_fc)) %>%
      dplyr::mutate(rank.fc = 1:nrow(merged_tb)) %>%
      dplyr::arrange(desc(gene_list_cor)) %>%
      dplyr::mutate(rank.cor = 1:nrow(merged_tb),
                    sign = ifelse(gene_list_fc > 0, yes = "up", no = "down"))



    ### define 0-score for ranks
    zero_fc = max(merged_tb[merged_tb$gene_list_fc > 0,]$rank.fc)
    zero_cor = max(merged_tb[merged_tb$gene_list_cor > 0,]$rank.cor)



    ### plot rank values
    rank_corr =
      ggplot(data = merged_tb,
           aes(x = rank.fc,
               y = rank.cor)) +
      geom_point(mapping = aes(color = sign),
                 stroke = NA,
                 size = 3,
                 alpha = 0.25,
                 show.legend = FALSE) +
      geom_vline(xintercept = zero_fc, color = "gray", linetype = 2) +
      geom_hline(yintercept = zero_cor, color = "gray", linetype = 2) +
      scale_color_manual(values = c("up" = color.up, "down" = color.down)) +
      ylab("Rank on Spearman's correlation") +
      xlab("Rank on log~2~(Fold Change)") +
      ggtitle(label = paste0(contrasts.info$metadata.column,": **",contrasts.info$var.1,"** *vs* **", contrasts.info$var.2,"**"),
              subtitle = "Ranks correlation") +
      geom_smooth(method = "glm", formula = y~x, color = regression.line.color, fill = regression.line.color) +
      scale_x_reverse() +
      scale_y_reverse() +
      ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho") +
      ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", vjust = 3) +
      theme_classic() +
      theme(aspect.ratio = 1,
            axis.text = element_text(color = "black"),
            axis.title = ggtext::element_markdown(color = "black"),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(color = "black", hjust = 0.5),
            axis.ticks = element_line(color = "black"),
            axis.line = element_blank(),
            panel.background = element_rect(fill = NA, color = "black", linewidth = 1))



    score_corr =
      ggplot(data = merged_tb,
             aes(x = gene_list_fc,
                 y = gene_list_cor)) +
      geom_point(mapping = aes(color = sign),
                 stroke = NA,
                 size = 3,
                 alpha = 0.25,
                 show.legend = F) +
      geom_vline(xintercept = 0, color = "gray", linetype = 2) +
      geom_hline(yintercept = 0, color = "gray", linetype = 2) +
      scale_color_manual(values = c("up" = color.up, "down" = color.down)) +
      ylab("Spearman's correlation coefficient (\u03c1)") +
      xlab("log~2~(Fold Change)") +
      ggtitle(label = paste0(contrasts.info$metadata.column,": **",contrasts.info$var.1,"** *vs* **", contrasts.info$var.2,"**"),
              subtitle = "Scores correlation") +
      geom_smooth(method = "glm", formula = y~x, color = regression.line.color, fill = regression.line.color) +
      ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho") +
      ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", vjust = 3) +
      theme_classic() +
      scale_y_continuous(limits = c(-1,1)) +
      theme(aspect.ratio = 1,
            axis.text = element_text(color = "black"),
            axis.title = ggtext::element_markdown(color = "black"),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(color = "black", hjust = 0.5),
            axis.ticks = element_line(color = "black"),
            axis.line = element_blank(),
            panel.background = element_rect(fill = NA, color = "black", linewidth = 1))


    ## return combined plot
    return(patchwork::wrap_plots(rank_corr, score_corr, nrow = 1))
  }
