#' @title perform.PCA
#'
#' @description This function performs principal component analyses (PCA).
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param sample.subset String vector indicating the column names (samples) to keep in the counts table (the 'column.id' in the metadata table). Default: \code{NULL} (no subsetting).
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'imp'. Default: \code{"imputed"}.
#' @param n.PCs Integer number indicating the number of PCs to be computed. This is used only when NAs are present in the the data set. Default: \code{10}.
#'
#' @return A \code{DEprot.PCA}, containing the PC values (\code{PCs}) and the importance summary (\code{importance}).
#'
#' @export perform.PCA

perform.PCA =
  function(DEprot.object,
           sample.subset = NULL,
           which.data = "imputed",
           n.PCs = 10) {

    ### Libraries
    require(dplyr)
    require(ggplot2)
    # require(pcaMethods)


    ### check object
    if (!("DEprot" %in% class(DEprot.object)) & !("DEprot.analyses" %in% class(DEprot.object))) {
      return(warning("The input must be an object of class 'DEprot'."))
    }

    ### Check and extract table
    if (tolower(which.data) == "raw") {
      if (!is.null(DEprot.object@raw.counts)) {
        mat = DEprot.object@raw.counts
        data.used = "raw"
      } else {
        warning(paste0("Use of RAW counts was required, but not available.\n",
                       "Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        return(DEprot.object)
      }
    } else if (tolower(which.data) %in% c("norm", "normalized", "normal")) {
      if (!is.null(DEprot.object@norm.counts)) {
        mat = DEprot.object@norm.counts
        data.used = "normalized"
      } else {
        warning(paste0("Use of NORMALIZED counts was required, but not available.\n",
                       "Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        return(DEprot.object)
      }
    } else if (tolower(which.data) %in% c("imputed", "imp", "impute")) {
      if (!is.null(DEprot.object@imputed.counts)) {
        mat = DEprot.object@imputed.counts
        data.used = "imputed"
      } else {
        warning(paste0("Use of IMPUTED counts was required, but not available.\n",
                       "Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        return(DEprot.object)
      }
    } else {
      warning(paste0("The 'which.data' value is not recognized.\n",
                     "Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
      return(DEprot.object)
    }



    ## Convert table to log2
    if (!is.numeric(DEprot.object@log.base)) {
      message("The log.base is not numeric, linear counts are assumed. Counts matrix will be converted to log2(score+1) values to analyze the data.")
      mat.log2 = log2(mat + 1)
    } else if (as.numeric(DEprot.object@log.base) != 2) {
      message("The log.base is not 2, counts will be converted to log2 values to analyze the data.")
      mat.log2 = log2(DEprot.object@log.base^mat)
    } else {
      mat.log2 = mat
    }

    ### subset table
    if (!is.null(sample.subset)) {
      mat = mat[,which(colnames(mat) %in% sample.subset)]
      PCA.meta = dplyr::filter(DEprot.object@metadata, column.id %in% sample.subset)
    } else {
      PCA.meta = DEprot.object@metadata
    }

    ### Convert all NaN/NA
    mat[is.nan(mat)] = NA


    ### Compute PCA
    # -----------------------------
    # Data without NAs (e.g., imputed data)
    if (!(TRUE %in% is.na(mat))) {
      pc = prcomp(mat, center = TRUE, scale. = TRUE)
      pc.summary = summary(pc)


      ### Combine PCA with metadata
      combo.tb = dplyr::left_join(dplyr::mutate(data.frame(pc.summary$rotation),
                                                column.id = rownames(data.frame(pc.summary$rotation))),
                                  PCA.meta,
                                  by = "column.id")


      ### Create importance tb
      importance.tb =
        data.frame(t(pc.summary$importance)) %>%
        dplyr::mutate(PC = factor(gsub("PC", "", rownames(t(pc.summary$importance))),
                                  levels = gsub("PC", "", rownames(t(pc.summary$importance)))),
                      Percentage.of.Variance = Proportion.of.Variance*100)

      # ---------------------------
      # For data with NA/NaN, such as raw/normalized data
    } else {
      pc = pcaMethods::pca(object = t(mat), method = "nipals", nPcs = round(n.PCs,0), center = T, scale. = T)


      ### Combine PCA with metadata
      combo.tb = dplyr::left_join(dplyr::mutate(data.frame(pc@scores),
                                                column.id = rownames(data.frame(pc@scores))),
                                  PCA.meta,
                                  by = "column.id")

      ### Create importance tb
      importance.tb =
        data.frame(Standard.deviation = pc@sDev,
                   Proportion.of.Variance = pc@R2,
                   Cumulative.Proportion = pc@R2cum,
                   PC = factor(gsub("PC", "", colnames(data.frame(pc@scores))),
                               levels = gsub("PC", "", colnames(data.frame(pc@scores)))),
                   Percentage.of.Variance = pc@R2 * 100)
    }



    ##################### PLOTS ##########################
    ### Plot cumulative plot
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
      new(Class = "DEprot.PCA",
          PCA.metadata = PCA.meta,
          sample.subset = paste(PCA.meta$column.id, collapse = ", "),
          data.used = which.data,
          prcomp = pc,
          PCs = combo.tb,
          importance = importance.tb,
          cumulative.PC.plot = cumulative_plot)

    return(DEprot.PCA.object)
  } # END function
