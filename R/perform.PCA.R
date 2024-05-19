#' @title perform.PCA
#'
#' @description This function performs principal component analyses (PCA).
#'
#' @param DEprot.object An object of class \code{DEprot}.
#' @param column.subset String vector indicating the column names to keep in the counts table (or 'column.id' in the metdata table). Default: \code{NULL} (no subsetting).
#' @param use.normalized.data Logical value to indicate whether to use the normalized data or not. Default: \code{TRUE}.
#'
#' @return A \code{DEprot.PCA}, containing the PC values (\code{PCs}) and the importance summary (\code{importance}).
#'
#' @export perform.PCA

perform.PCA =
  function(DEprot.object,
           column.subset = NULL,
           use.normalized.data = TRUE) {

    ### Libraries
    require(dplyr)
    require(ggplot2)


    ### check object
    if (!("DEprot" %in% class(DEprot.object))) {
      return(warning("The input must be an object of class 'DEprot'."))
    }

    ### Check and extract table
    if (use.normalized.data == TRUE) {
      if (!is.null(DEprot.object$norm.counts)) {
        mat = DEprot.object$norm.counts
      } else {
        return(warning("Use of normalized counts have been required, but not normalized data are available."))
      }
    } else {
      if (!is.null(DEprot.object$raw.counts)) {
        mat = DEprot.object$raw.counts
      } else {
        return(warning("Use of raw counts have been required, but not raw data are available."))
      }
    }

    ### subset table
    if (!is.null(column.subset)) {
      mat = mat[,which(colnames(mat) %in% column.subset)]
      PCA.meta = dplyr::filter(DEprot.object$metadata, column.id %in% column.subset)
    } else {
      PCA.meta = DEprot.object$metadata
    }


    ### Compute PCA
    pc = prcomp(mat, center = TRUE, scale. = TRUE)
    pc.summary = summary(pc)


    ### Combine PCA with metadata
    combo.tb = dplyr::left_join(dplyr::mutate(data.frame(pc.summary$rotation),
                                              column.id = rownames(data.frame(pc.summary$rotation))),
                                PCA.meta,
                                by = "column.id")


    ### Plot cumulative plot
    importance.tb =
      data.frame(t(pc.summary$importance)) %>%
      dplyr::mutate(PC = factor(gsub("PC", "", rownames(t(pc.summary$importance))),
                                levels = gsub("PC", "", rownames(t(pc.summary$importance)))),
                    Percentage.of.Variance = Proportion.of.Variance*100)

    cumulative_plot =
      ggplot(data = importance.tb,
             aes(x = PC)) +
      geom_bar(mapping = aes(y = Proportion.of.Variance),
               stat = "identity",
               fill = "steelblue") +
      geom_line(mapping = aes(y = Cumulative.Proportion,
                              group = 1),
               color = "navyblue",
               linetype = 1) +
      geom_point(mapping = aes(y = Cumulative.Proportion,
                              group = 1),
                stroke = NA,
                size = 2,
                color = "navyblue") +
      ylab("Proportion of variance") +
      xlab("Principal Component (PC)") +
      theme_classic() +
      scale_y_continuous(expand = c(0,0)) +
      theme(axis.ticks.x = element_blank(),
            axis.text = element_text(color = "black"),
            axis.line.x = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.grid.major.y = element_line(color = "gray", linewidth = 0.1),
            panel.grid.minor.y = element_line(color = "gray", linewidth = 0.1),
            aspect.ratio = 0.5)


    DEprot.PCA.object =
      structure(list(PCs = combo.tb,
                     importance = importance.tb),
                # classes
                class = c("DEprot", "DEprot.PCA"),
                # packgs
                package = "DEprot",
                # attributes
                cumulative.PC.plot = cumulative_plot,
                PCA.metadata = PCA.meta)


    return(DEprot.PCA.object)
  } # END function
