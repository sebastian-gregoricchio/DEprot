#' @title plot.counts
#'
#' @description Plots a volcano plot log2(FoldChange) x -log10(p.adjusted) of differential expression results
#'
#' @param DEprot.object An object of class \code{DEprot.analyses}.
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'imp'. Default: \code{"imputed"}.
#' @param violin.color String indicating the color to use for the violin plots. Default: \code{"darkorange"}.
#' @param max.line.color String indicating the color to use for up-regulated proteins in the plots. Default: \code{"indianred"}.
#' @param min.line.color String indicating the color to use for up-regulated proteins in the plots. Default: \code{"steelblue"}.
#' @param title String indicating the title to use. Default: \code{NULL} (automatic title).
#' @param convert.log2 Logical value to define whether counts should be log2 transformed. Default: \code{TRUE}.
#'
#' @return A ggplot object.
#'
#' @export plot.counts

plot.counts =
  function(DEprot.object,
           which.data = "imputed",
           violin.color = "darkorange",
           max.line.color = "indianred",
           min.line.color = "steelblue",
           title = NULL,
           convert.log2 = TRUE) {

    ### Libraries
    require(dplyr)
    require(ggplot2)
    require(patchwork)


    ### check object and extract metadata table
    if (!("DEprot" %in% class(DEprot.object))) {
      if (!("DEprot.analyses" %in% class(DEprot.object))) {
        warning("The input must be an object of class 'DEprot' or 'DEprot.analyses'.")
        return(DEprot.object)
      }
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



    ### Convert to log if required
    if (convert.log2 == TRUE & DEprot.object@log.base != 2) {
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
      log.base = DEprot.object@log.base
    } else {
      mat.log2 = mat
      log.base = DEprot.object@log.base
    }



    ### Generate boxplot of counts
    # melt counts table
    melt.cnt =
      suppressMessages(reshape2::melt(as.data.frame(mat.log2))) %>%
      dplyr::mutate(variable = factor(variable, levels = colnames(mat.log2)))

    # compute stats
    cnt.stats =
      melt.cnt %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(min = min(value, na.rm = T),
                       max = max(value, na.rm = T))

    boxplot =
      ggplot() +
      geom_violin(data = melt.cnt,
                  mapping = aes(x = variable,
                                y = value,
                                group = variable),
                  width = 0.75,
                  alpha = 0.75,
                  fill = violin.color,
                  color = NA) +
      geom_boxplot(data = melt.cnt,
                   mapping = aes(x = variable,
                                 y = value,
                                 group = variable),
                   fill = "white",
                   color = colorspace::darken(violin.color, amount = 0.3),
                   width = 0.15,
                   outlier.color = "black",
                   outlier.stroke = NA,
                   outlier.size = 2,
                   outlier.alpha = 0.25) +
      geom_line(data = data.frame(cnt.stats),
                mapping = aes(x = variable,
                              y = max,
                              group = 1),
                color = max.line.color,
                linetype = 2,
                inherit.aes = F) +
      geom_line(data = data.frame(cnt.stats),
                mapping = aes(x = variable,
                              y = min,
                              group = 1),
                color = min.line.color,
                linetype = 2,
                inherit.aes = F) +
      ylab(ifelse(is.na(log.base),
                  yes = "Intensity",
                  no = paste0(ifelse(log.base == exp(1),
                                     yes = "ln", no = paste0("log<sub>",log.base,"</sub>")),
                              "(Intensity)"))) +
      ggtitle(title) +
      xlab("Sample") +
      theme_classic() +
      theme(axis.text.y = element_text(color = "black"),
            axis.text.x = element_text(color = "black", hjust = 1, angle = 30),
            axis.title = ggtext::element_markdown(color = "black"),
            axis.ticks.y = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            aspect.ratio = 10/ncol(mat))


    ### Export plot
    return(boxplot)
  } # END of function

