#' @title geneset.enrichment
#'
#' @description Perform Gene Set Enrichment Analyses (GSEA) or OverRepresentation Analyses (ORA) on a proteomics differential analyses.
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast Number indicating the position of the contrast to use for the plotting.
#' @param enrichment.type String indicating the type of analyses to perform. One among: GSEA, ORA. No default.
#' @param diff.status.category String indicating a diff.status among the ones present in the results table of the specific contrast. Used only one 'ORA' is performed. Default: \code{NULL}.
#' @param gsub.pattern.prot.id String indicating a pattern to be passed to gsub and to remove from the prot.id. Default: \code{NULL} (non changes in the IDs).
#' @param pvalueCutoff Numeric value indicating the adjusted pvalue cutoff on enrichment tests to report. Default: \code{0.05}.
#' @param qvalueCutoff Numeric value indicating the qvalue cutoff on enrichment tests to report as significant (only for ORA). Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported. Default: \code{0.05}.
#' @param pAdjustMethod String indicating the method to use for the p-value adjustment. One mong "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: \code{"BH"}.
#' @param dotplot.n Numeric value indicating the maximum number of categories to plot in the dotplot. Default: \code{10}.
#'
#' @return An object of class \code{DEprot.enrichResult}.

#' @export geneset.enrichment

geneset.enrichment =
  function(DEprot.analyses.object,
           contrast,
           TERM2GENE,
           enrichment.type,
           diff.status.category = NULL,
           gsub.pattern.prot.id = NULL,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           dotplot.n = 10) {

    ### Libraries
    require(dplyr)
    require(ggplot2)
    #require(clusterProfiler)
    #require(aPEAR)

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      warning("The 'clusterProfiler' (Bioconductor) package must be installed to use this function.")

      ### Ask for installing
      install = readline("Do you want to install `clusterProfiler`? [yes/no] ")
      if (tolower(install) %in% c("yes","y","yeah","yep","yo")) {
        BiocManager::install("clusterProfiler")
        library(clusterProfiler)
      } else {
        return(DEprot.object)
      }
    } else {
      require(clusterProfiler)
    }



    if (!requireNamespace("aPEAR", quietly = TRUE)) {
      warning("The 'aPEAR' (CRAN) package must be installed to use this function.")

      ### Ask for installing
      install = readline("Do you want to install `aPEAR`? [yes/no] ")
      if (tolower(install) %in% c("yes","y","yeah","yep","yo")) {
        install.packages("aPEAR")
        library(aPEAR)
      } else {
        return(DEprot.object)
      }
    } else {
      require(aPEAR)
    }

    ######################################################################################
    # Functions

    plot_NES = function(gsea.object,
                        pos.NES.label = "+NES",
                        neg.NES.label = "-NES",
                        pos.NES.color = "steelblue",
                        neg.NES.color = "orange",
                        string.pattern.to.remove = "HALLMARK|GOBP",
                        alpha.range = c(0.3,1),
                        add.counts = TRUE,
                        perc.bleeding.x = 8,
                        axes.text.size = 10,
                        title = "NES enrichments") {

      # libraries
      require(dplyr)
      require(ggplot2)


      # extract results and clean
      result =
        gsea.object@result %>%
        dplyr::arrange(desc(NES)) %>%
        dplyr::mutate(alias = gsub("_", " ", gsub(string.pattern.to.remove, "", ID)),
                      dataset = ifelse(NES >= 0, yes = pos.NES.label, no = neg.NES.label)) %>%
        dplyr::mutate(alias = factor(alias, levels = rev(alias)))

      geneSets_size = data.frame(sapply(gsea.object@geneSets, length), stringsAsFactors = F)
      geneSets_size$ID = rownames(geneSets_size)
      colnames(geneSets_size)[1] = "geneSet_size"

      result = dplyr::left_join(result, geneSets_size, by = "ID")


      # define colors
      NES.colors = c(pos.NES.color, neg.NES.color)
      names(NES.colors) = c(pos.NES.label, neg.NES.label)

      max_abs.NES = ceiling(max(abs(result$NES), na.rm = T))

      # Generating the plot
      NES.plot =
        ggplot(data = result,
               aes(x = NES,
                   y = alias,
                   fill = factor(dataset, levels = c(pos.NES.label, neg.NES.label)))) +
        ggplot2::geom_bar(aes(alpha = -log10(p.adjust)),
                          stat = "identity",
                          show.legend = T,
                          width = 0.8) +
        scale_alpha_continuous(range = alpha.range) +
        scale_fill_manual(values = NES.colors, name = "dataset", drop = F) +
        ylab(NULL) +
        ggtitle(title) +
        theme_classic() +
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text = element_text(color = "black",
                                       size = axes.text.size),
              panel.background = element_blank(),
              plot.title = ggtext::element_markdown(hjust = 0.5))


      # Add counts
      if (add.counts == TRUE) {
        if (nrow(result %>% dplyr::filter(NES >= 0)) > 0) {
          NES.plot =
            NES.plot +
            geom_text(data = result %>% dplyr::filter(NES >= 0),
                      aes(x = NES,
                          y = alias,
                          label = paste0(setSize,"/",geneSet_size)),
                      vjust = 0.5,
                      hjust = -0.2)
        }

        if (nrow(result %>% dplyr::filter(NES < 0)) > 0) {
          NES.plot =
            NES.plot +
            geom_text(data = result %>% dplyr::filter(NES < 0),
                      aes(x = NES,
                          y = alias,
                          label = paste0(setSize,"/",geneSet_size)),
                      vjust = 0.5,
                      hjust = 1.2)
        }
      } else {
        perc.bleeding.x = 0
      }

      NES.plot =
        NES.plot +
        #xlim(c(-1,1)*(max_abs.NES + (perc.bleeding.x/100)*max_abs.NES)) +
        scale_x_continuous(limits = c(-1,1)*(max_abs.NES + (perc.bleeding.x/100)*max_abs.NES),
                           expand = c(0,0))


      # export data
      return(NES.plot)

    } #END plot_NES function



    ######################################################################################


    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      return(warning("The input must be an object of class 'DEprot.analyses'."))
    }

    ### check and collect contrast
    if (is.numeric(contrast)) {
      if (contrast <= length(DEprot.analyses.object@analyses.result.list)) {
        data = DEprot.analyses.object@analyses.result.list[[contrast]]$results
        contrasts.info = DEprot.analyses.object@contrasts[[contrast]]
      } else {
        return(warning("The 'contrast' indicated is not available."))
      }
    } else {
      return(warning("The 'contrast' must be a numeric value."))
    }


    ### Check enrichment.type
    if (!(toupper(enrichment.type) %in% c("ORA", "GSEA", "OVER"))) {
      return(warning("The 'enrichment.type' must be one among: 'ORA', 'GSEA'."))
    }


    ### Check TERM2GENE
    if (!("data.frame" %in% class(TERM2GENE))) {
      return(warning("The 'enrichment.type' must be one among: 'ORA', 'GSEA'."))
    } else {
      if (!("gs_name" %in% colnames(TERM2GENE))) {
        return(warning("The 'TERM2GENE' must be a 2-column data.frame. One of this columns must be named 'gs_name' (gs: gene set)"))
      } else if (ncol(TERM2GENE) > 2) {
        return(warning("The 'TERM2GENE' must be a 2-column data.frame. One of this columns must be named 'gs_name' (gs: gene set)"))
      }
    }


    ### Check diff.status.category
    if (toupper(enrichment.type) != "GSEA") {
      if (!is.null(diff.status.category)) {
        if (!(diff.status.category %in% unique(data$diff.status))) {
          return(warning(paste0("The 'diff.status.category' must be a value of the 'diff.status' column in the results.\n\n",
                                "For contrast #", contrast, ": ", paste(unique(data$diff.status), collapse = ", "),".")))
        }
      } else {
        return(warning(paste0("The 'diff.status.category' must be a value of the 'diff.status' column in the results.\n\n",
                              "For contrast #", contrast, ": ", paste(unique(data$diff.status), collapse = ", "),".")))
      }
    }


    ### gsub the protein.id if necessary
    if (!is.null(gsub.pattern.prot.id)) {
      data$prot.id = gsub(gsub.pattern.prot.id, "", data$prot.id)
    }




    #####################################################################################

    #### perform enrichments
    if (toupper(enrichment.type) == "GSEA") {
      gene_list = data[,5]
      names(gene_list) = data$prot.id
      gene_list = sort(gene_list, decreasing = T)


      enrichment.discovery = tryCatch(clusterProfiler::GSEA(geneList = gene_list,
                                                            pvalueCutoff = pvalueCutoff,
                                                            pAdjustMethod = pAdjustMethod,
                                                            TERM2GENE = TERM2GENE),
                                      error = function(x)(return(NULL)))


      ## plot networks
      protein.network = tryCatch(clusterProfiler::cnetplot(enrichment.discovery), error = function(x)(return(NULL)))

      pathway.network.clusters = tryCatch(aPEAR::findPathClusters(enrichment.discovery@result),error = function(x)(return(NULL)))
      pathway.network.plot = tryCatch(aPEAR::enrichmentNetwork(enrichment.discovery@result, drawEllipses = TRUE, colorBy = "NES"),
                                      error = function(x)(return(NULL)))

      NES.plot = tryCatch(plot_NES(gsea.object = enrichment.discovery,
                                   pos.NES.label = contrasts.info$var.1,
                                   neg.NES.label = contrasts.info$var.2,
                                   title = paste0(contrasts.info$metadata.column,": **", contrasts.info$var.1, "** *vs* **", contrasts.info$var.2,"**")),
                          error = function(x)(return(NULL)))

      dotplot_fold.enrichment = NULL

    } else {
      enrichment.discovery = clusterProfiler::enricher(gene = dplyr::filter(.data = data,
                                                                            diff.status == diff.status.category)$prot.id,
                                                       pvalueCutoff = pvalueCutoff,
                                                       qvalueCutoff = qvalueCutoff,
                                                       pAdjustMethod = pAdjustMethod,
                                                       TERM2GENE = TERM2GENE)

      ## plot networks
      protein.network = tryCatch(clusterProfiler::cnetplot(enrichment.discovery), error = function(x)(return(NULL)))

      pathway.network.clusters = tryCatch(aPEAR::findPathClusters(enrichment.discovery@result),error = function(x)(return(NULL)))
      pathway.network.plot = tryCatch(aPEAR::enrichmentNetwork(enrichment.discovery@result, colorBy = 'pvalue', colorType = 'pval', pCutoff = -5),
                                      error = function(x)(return(NULL)))

      NES.plot = NULL


      dotplot_fold.enrichment =
        tryCatch(clusterProfiler::dotplot(enrichment.discovery,
                                          x = "FoldEnrichment",
                                          showCategory = dotplot.n) +
                   ggtitle(paste0(contrasts.info$metadata.column,": **", contrasts.info$var.1, "** *vs* **", contrasts.info$var.2,"**")) +
                   viridis::scale_fill_viridis(option = "mako", direction = -1, begin = 0.3) +
                   theme(plot.title = ggtext::element_markdown(hjust = 0.5),
                         axis.ticks.y = element_blank()),
                 error = function(x)(return(NULL)))
    }


    ## plot dotplot
    dotplot_gene.ratio =
      tryCatch(clusterProfiler::dotplot(enrichment.discovery,
                                        x = "GeneRatio",
                                        showCategory = dotplot.n) +
                 ggtitle(paste0(contrasts.info$metadata.column,": **", contrasts.info$var.1, "** *vs* **", contrasts.info$var.2,"**")) +
                 viridis::scale_fill_viridis(option = "rocket", direction = -1, begin = 0.3) +
                 theme(plot.title = ggtext::element_markdown(hjust = 0.5),
                       axis.ticks.y = element_blank()))

    ## Pathway network list
    pathway.network = list(clusters = pathway.network.clusters,
                           plot = pathway.network.plot)



    ##################################################
    ### build object
    DEprot.enrichResult.object =
      new(Class = "DEprot.enrichResult",
          enrichment.discovery = enrichment.discovery,
          protein.network = protein.network,
          pathway.network = pathway.network,
          NES.plot = NES.plot,
          dotplot_gene.ratio = dotplot_gene.ratio,
          dotplot_fold.enrichment = dotplot_fold.enrichment,
          parameters = list(enrichment.type = enrichment.type,
                            diff.status.category = diff.status.category,
                            gsub.pattern.prot.id = gsub.pattern.prot.id,
                            pvalueCutoff = pvalueCutoff,
                            qvalueCutoff = qvalueCutoff,
                            pAdjustMethod = pAdjustMethod,
                            dotplot.n = dotplot.n))

    return(DEprot.enrichResult.object)
  } #END function
