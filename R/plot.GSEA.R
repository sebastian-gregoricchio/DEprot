#' @title plot.GSEA
#'
#' @description Function to plot GSEA results (see \href{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}{clusterprofiler}).
#'
#' @param gsea.results Object of class \code{DEprot.enrichResult}, generated using \link{geneset.enrichment}, or \code{gseaResult}, generate by \href{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}{clusterprofiler}.
#' @param geneset.id Numeric value or a string identifying the Nth geneSet (numeric) or a specific id (string) if geneSet in the result table. Default value: \code{NULL}, which returns the ordered list of available geneSets.
#' @param enrichment.geom String indicating the type of graph to use to plot the enrichment scores. Possible options: 'line', 'lines', 'dot', 'dots', 'point', 'points' (case insensitive). Default: \code{'line'}.
#' @param enrichment.color String indicating any R-supported color to be used for the enrichment score plot. Default: \code{'green'}.
#' @param enrichment.geom.size Numeric value indicating the size of the line, or dots, used in the enrichment score plot. Default: \code{1}.
#' @param enrichment.plot.zero.line Logical value to indicated whether to plot an horizontal line at 0 in the enrichment score plot. Default: \code{TRUE}.
#' @param enrichment.zero.line.color String indicating any R-supported color to be used for the 0-line in the enrichment score plot (active when \code{enrichment.plot.zero.line = TRUE}). Default: \code{'gray'}.
#' @param enrichment.zero.line.width Numeric value indicating the line width of the 0-line in the enrichment score plot (active when \code{enrichment.plot.zero.line = TRUE}). Default: \code{0.5}.
#' @param enrichment.annotations.vjust.offset Numeric value to add to the vjust (vertical positioning) of the enrichment plot annotations (P, Padj, q, NES, set size). Positive values will shift-down the annotations. Default: \code{0} (base line).
#' @param geneset.segments.width Numeric value indicating the line width of the geneSet vertical segments. Default: \code{0.3}.
#' @param geneset.segments.color String indicating any R-supported color to be used for the geneSet segments. Default: \code{'black'}.
#' @param add.condition.labels Logic value indicating whether the labels of the conditions compared should be plotted. Ignored if the object is not of class \code{DEprot.enrichResult}. Default: \code{TRUE}.
#' @param condition.label.positive String indicating the label to use for the condition showing a positive ranking score (left side of the plot). Ignored if \code{add.condition.labels = FALSE}. Default: \code{NULL}, none or automatic.
#' @param condition.label.negative String indicating the label to use for the condition showing a negative ranking score (right side of the plot). Ignored if \code{add.condition.labels = FALSE}. Default: \code{NULL}, none or automatic.
#' @param rank.max.color String indicating any R-supported color to be used for the max rank dotted lines and annotation. Default: \code{'indianred'}.
#' @param ranking.color String indicating any R-supported color to be used for the ranked list plot (histogram). Default: \code{'gray'}.
#' @param gradient.colors Two-values string vector indicating the shadows of palettes to use for the genset gradient. Possible values: 'Blues', 'Greens', 'Greys', 'Oranges', 'Purples', 'Reds'. Default: \code{c('Reds', 'Blues')}.
#' @param title.position String indicating the position of the title: 'left', 'center', 'right'. Default: \code{'center'}.
#' @param title String indicating the title to use. Default: \code{NA}, this will automatically use the geneset name chosen. Use \code{NULL} to do not plot the title.
#' @param image.file.name String indicating the full path for the export of a pdf file of the combined plot. Default: \code{NULL}, no plot will be exported.
#' @param image.width Numeric value to indicate the width (in inches) to use for the exported pdf file. Active only when \code{image.file.name} is not \code{NULL}. Default: \code{7}.
#' @param image.height Numeric value to indicate the height (in inches) t use for the exported pdf file. Active only when \code{image.file.name} is not \code{NULL}. Default: \code{5}.
#' @param return.all.objects Logical value to indicate whether the function should return only the combined plot (ggplot object), or all the different panels and the combined plot in a list. Default: \code{FALSE} (only combined plot).
#'
#' @return Either a ggplot-object with the final combined plot, or a list with the three panels separated and the combined plot: \code{list(enrichment.panel, geneset.panel, rank.panel, combined.plot)}.
#'
#' @seealso The internal function \code{gsInfo} is derived from the package \href{https://github.com/YuLab-SMU/enrichplot/blob/bf86dcdc74316c7ee8d9e7f5fcac7e88392204e6/R/gseaplot.R#L86}{enrichplot::gsInfo}.
#'
#' @aliases plot.GSEA
#'
#' @name plot.GSEA
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices dev.off pdf
#' @importFrom utils getFromNamespace
#'
#' @author Sebastian Gregoricchio
#'
#' @export plot.GSEA


plot.GSEA =
  function(gsea.results,
           geneset.id = NULL,
           enrichment.geom = "line",
           enrichment.color = "green",
           enrichment.geom.size = 1,
           enrichment.plot.zero.line = TRUE,
           enrichment.zero.line.color = "gray",
           enrichment.zero.line.width = 0.5,
           enrichment.annotations.vjust.offset = 0,
           geneset.segments.width = 0.3,
           geneset.segments.color = "black",
           add.condition.labels = TRUE,
           condition.label.positive = NULL,
           condition.label.negative = NULL,
           rank.max.color = "indianred",
           ranking.color = "gray",
           gradient.colors = c("Reds", "Blues"),
           title.position = "center",
           title = NA,
           image.file.name = NULL,
           image.width = 7,
           image.height = 5,
           return.all.objects = FALSE) {

    # # Libraries
    # require(ggplot2)



    #####################################
    # gsInfo from enrichplot
    gsInfo <- function(object, geneSetID) {
      geneList <- object@geneList

      if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

      geneSet <- object@geneSets[[geneSetID]]
      exponent <- object@params[["exponent"]]
      gseaScores <- getFromNamespace("gseaScores", "DOSE")
      df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
      df$ymin <- 0
      df$ymax <- 0
      pos <- df$position == 1
      h <- diff(range(df$runningScore))/20
      df$ymin[pos] <- -h
      df$ymax[pos] <- h
      df$geneList <- geneList
      if (length(object@gene2Symbol) == 0) {
        df$gene <- names(geneList)
      } else {
        df$gene <- object@gene2Symbol[names(geneList)]
      }

      df$Description <- object@result[geneSetID, "Description"]
      return(df)
    }


    #####################################


    # collect results
    if ("DEprot.enrichResult" %in% class(gsea.results)) {
      # Check if it is a GSEA analysis
      if (toupper(gsea.results@parameters$enrichment.type) != "GSEA") {
        stop("This function can plot only analyses performed using the 'GSEA' mode.")
        #return()
      }
      results = data.frame(gsea.results@enrichment.discovery@result)
      contrast = gsea.results@parameters$contrast
      gsea.results = gsea.results@enrichment.discovery
      mode = "DEprot"
    } else if ("gseaResult" %in% class(gsea.results)) {
      results = data.frame(gsea.results@result)
      mode = "clusterProfiler"
    } else {
      stop("The input of this functions must be an object of class 'DEprot.enrichResult' - as generated using the function `geneset.enruchment` - or ''.")
      #return()
    }


    # check geneset.id
    if (is.numeric(geneset.id)) {
      if (geneset.id <= nrow(results) & geneset.id > 0) {gene.set = results$ID[geneset.id]
      } else {
        print(data.frame(gene.set.ID = results$ID))
        stop("The genset.id must be a string or a numeric value among the ones above.")
        #return(invisible())
      }
    } else if (is.character(geneset.id)) {
      if (geneset.id %in% results$ID) {gene.set = geneset.id
      } else {
        print(data.frame(gene.set.ID = results$ID))
        stop("The genset.id must be a string or a numeric value among the ones above.")
        #return(invisible())
      }
    } else {
      print(data.frame(gene.set.ID = results$ID))
      stop("The genset.id must be a string or a numeric value among the ones above.")
      return(invisible())
    }

    # Check graphic parameters
    if (!(tolower(enrichment.geom) %in% c("line", "lines", "dot", "dots", "point", "points"))) {
      stop("The 'enrichment.geom' must be one among: 'line', 'lines', 'dot', 'dots', 'point', 'points'.")
    }

    if (!is.null(title)) {
      if (is.na(title)) {
        title = gene.set
      }
    }

    if (tolower(title.position) == "right") {
      title.position = 1
    } else if (tolower(title.position) == "left") {
      title.position = 0
    } else {
      title.position = 0.5
    }

    # Get plotting info
    info = gsInfo(object = gsea.results, geneSetID = gene.set)

    # Enrichment panel
    if (enrichment.geom %in% c("line", "lines")) {
      enrichment.panel =
        ggplot(data = info,
               aes(x = x,
                   y = runningScore)) +
        geom_line(color = enrichment.color,
                  linewidth = enrichment.geom.size)
    } else {
      enrichment.panel =
        ggplot(data = info,
               aes(x = x,
                   y = runningScore)) +
        geom_point(color = enrichment.color,
                   size = enrichment.geom.size)
    }

    enrichment.panel =
      enrichment.panel +
      ggtitle(title) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      geom_vline(xintercept = results$rank[results$ID == gene.set],
                 linetype = "dotted",
                 color = rank.max.color) +
      xlab(NULL) +
      ylab("Running Enrichment Score (ES)") +
      theme_classic() +
      annotate(geom = "text",
               x = +Inf,
               y = +Inf,
               label = deparse(bquote(italic("P =")~.(signif(results$pvalue[results$ID == gene.set], digits = 2)))),
               #label = paste0("P = ", signif(results$pvalue[results$ID == gene.set], digits = 2)),
               hjust = 1,
               vjust = 1+enrichment.annotations.vjust.offset,
               parse=TRUE) +
      annotate(geom = "text",
               x = +Inf,
               y = +Inf,
               label = deparse(bquote(italic(P[adj])~italic("=")~.(signif(results$p.adjust[results$ID == gene.set], digits = 3)))),
               #label = paste0("P adj = ", signif(results$p.adjust[results$ID == gene.set], digits = 3)),
               hjust = 1,
               vjust = 2.25+enrichment.annotations.vjust.offset,
               parse=TRUE) +
      annotate(geom = "text",
               x = +Inf,
               y = +Inf,
               label = deparse(bquote(italic("q =")~.(signif(results[results$ID == gene.set,grep("q.*val", colnames(results))], digits = 3)))),
               #label = paste0("q = ", signif(results$qvalues[results$ID == gene.set], digits = 3)),
               hjust = 1,
               vjust = 4+enrichment.annotations.vjust.offset,
               parse=TRUE) +
      annotate(geom = "text",
               x = +Inf,
               y = +Inf,
               label = paste0("NES = ", signif(results$NES[results$ID == gene.set], digits = 3)),
               hjust = 1,
               vjust = 7.5+enrichment.annotations.vjust.offset) +
      annotate(geom = "text",
               x = +Inf,
               y = +Inf,
               label = paste0("set size: ", results$setSize[results$ID == gene.set]),
               hjust = 1,
               vjust = 9.5+enrichment.annotations.vjust.offset) +
      theme(axis.text.y = element_text(color = "black"),
            axis.ticks.y = element_line(color = "black"),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = ggplot2::margin(c(0,0,0,0)),
            plot.title = element_text(hjust = title.position))

    if (enrichment.plot.zero.line == TRUE) {
      enrichment.panel =
        enrichment.panel +
        geom_hline(yintercept = 0,
                   linetype = 1,
                   linewidth = enrichment.zero.line.width,
                   color = enrichment.zero.line.color)
    }

    #enrichment.panel.build = ggplot_build(enrichment.panel)
    #enrichment.panel.build$layout$panel_params[[1]]$x.range

    # Gene set panel
    geneset.panel =
      ggplot(data = info) +
      geom_segment(aes(x = x,
                       xend = x,
                       y = ymin,
                       yend = ymax),
                   color = geneset.segments.color,
                   linewidth = geneset.segments.width) +
      scale_x_continuous(expand = c(0,0),
                         limits = range(info$x)) +
      scale_y_continuous(expand = c(0,0)) +
      geom_vline(xintercept = results$rank[results$ID == gene.set],
                 linetype = "dotted",
                 color = rank.max.color) +
      xlab(NULL) +
      ylab(NULL) +
      theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line.x = element_blank(),
            plot.margin = ggplot2::margin(c(0,0,0,0)))

    # Computing percentiles gradients adn add it ot the
    v = seq(1, sum(info$position), length.out=9)
    inv = findInterval(rev(cumsum(info$position)), v)
    if (min(inv) == 0) inv <- inv + 1
    col = c(rev(RColorBrewer::brewer.pal(5, gradient.colors[2])),
            RColorBrewer::brewer.pal(5, gradient.colors[1]))

    xmin = which(!duplicated(inv))
    gradient.df = data.frame(ymin = min(info$ymin),
                             ymax = min(info$ymin) * 0.5, #calculated on bottom half, so use 2*wanted fraction
                             xmin = xmin,
                             xmax = xmin + as.numeric(table(inv)[as.character(unique(inv))]),
                             col = col[unique(inv)])

    gradient.df$xmax[nrow(gradient.df)] = max(info$x)

    geneset.panel =
      geneset.panel +
      geom_rect(data = gradient.df,
                aes(xmin = xmin,
                    xmax = xmax,
                    ymin = ymin,
                    ymax = ymax,
                    fill = I(col)),
                alpha = 0.85,
                inherit.aes = FALSE) +
      theme(plot.margin = ggplot2::margin(c(0,0,0,0)))

    # rank panel
    ## compute 0-cross
    info$abs = abs(info$geneList)
    zero.cross = info$x[info$abs == min(info$abs)]
    if (length(zero.cross) > 1) {zero.cross = zero.cross[ceiling(length(zero.cross) /2)]}
    max.rank = results$rank[results$ID == gene.set]
    hjust.rank = ifelse(test = max.rank < (nrow(info)/5),
                        yes = -0.1,
                        no = 1.1)

    rank.panel =
      ggplot(data = info,
             aes(x = x,
                 y = geneList)) +
      #geom_area(fill = ranking.color) +
      geom_bar(stat = "identity", fill = ranking.color) +
      geom_vline(xintercept = results$rank[results$ID == gene.set],
                 linetype = "dotted",
                 color = rank.max.color) +
      geom_vline(xintercept = zero.cross-0.5,
                 linetype = "dashed",
                 color = ranking.color) +
      scale_x_continuous(expand = c(0,0),
                         limits = range(info$x)) +
      scale_y_continuous(expand = c(0,0))+
      xlab("Rank in ordered Dataset") +
      ylab("Ranked list metric") +
      theme_classic() +
      annotate(geom = "text",
               x = max.rank,
               y = -Inf,
               label = paste0("Rank max: ", max.rank),
               hjust = hjust.rank,
               vjust = -0.5,
               color = rank.max.color) +
      annotate(geom = "text",
               x = zero.cross,
               y = 0,
               label = paste0("Zero score at ", zero.cross),
               hjust = 0.5,
               vjust = -1) +
      theme(axis.text = element_text(color = "black"),
            axis.ticks = element_line(color = "black"),
            plot.margin = ggplot2::margin(c(0,0,0,0)))


    # add condition labels to rank panel
    if (add.condition.labels == TRUE) {
      if (mode == "DEprot") {
        rank.panel =
          rank.panel +
          annotate(geom = "text",
                   x = -Inf,
                   y = +Inf,
                   label = ifelse(test = is.null(condition.label.positive),
                                  yes = contrast$var.1,
                                  no = as.character(condition.label.positive)),
                   hjust = -0.05,
                   vjust = 2,
                   color = gradient.df$col[1]) +
          annotate(geom = "text",
                   x = +Inf,
                   y = +Inf,
                   label = ifelse(test = is.null(condition.label.negative),
                                  yes = contrast$var.2,
                                  no = as.character(condition.label.negative)),
                   hjust = 1.05,
                   vjust = 2,
                   color = gradient.df$col[nrow(gradient.df)])
      } else {
        if (!is.null(condition.label.positive)) {
          rank.panel =
            rank.panel +
            annotate(geom = "text",
                     x = -Inf,
                     y = +Inf,
                     label = as.character(condition.label.positive),
                     hjust = -0.1,
                     vjust = 2,
                     color = gradient.df$col[1])
        }

        if (!is.null(condition.label.negative)) {
          rank.panel =
            rank.panel +
            annotate(geom = "text",
                     x = +Inf,
                     y = +Inf,
                     label = as.character(condition.label.negative),
                     hjust = 1.1,
                     vjust = 2,
                     color = gradient.df$col[nrow(gradient.df)])
        }
      }
    }



    # COMBINING the plots
    combined.plot =
      patchwork::wrap_plots(enrichment.panel, geneset.panel, rank.panel,
                            ncol = 1,
                            heights = c(1,0.3,0.6))

    # Export image if required
    if (!is.null(image.file.name)) {
      pdf(file = image.file.name, width = image.width, height = image.height)
      print(combined.plot)
      invisible(dev.off())
    }

    if (return.all.objects == TRUE) {
      return(list(enrichment.panel = enrichment.panel,
                  geneset.panel = geneset.panel,
                  rank.panel = rank.panel,
                  combined.plot = combined.plot))
    } else {
      return(combined.plot)
    }

  } # END FUNCTION
