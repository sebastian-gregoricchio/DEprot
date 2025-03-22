#' @title expression.boxplot
#'
#' @description Plots a boxplot of the expression of a specific protein. Samples can be groups depending on a metadata column.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param protein.id String indicating a protein for which plot the expression. The identifier must correspond to the full row.name of the counts table (equivalent to the \code{prot.id} column of the fold change table of \code{DEprot.analyses} object).
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'imp'. Default: \code{"imputed"}.
#' @param sample.subset Character vector indicating a subset of samples to display. The identifiers must correspond to a IDs in the \code{column.id} column of the object's metadata. Default: \code{NULL} (all samples are shown).
#' @param shape.column String indicating a column from the metadata table. This column will be used as factor for the shape of the points on the boxplot. Default: \code{NULL}: no different shapes.
#' @param group.by.metadata.column String indicating a column from the metadata table. This column will be used to define sample groups, and for each group it will be computed a mean of the counts. Default: \code{"column.id"} (no groups).
#' @param group.levels Ordered string vector indicating the order to use for the groups. Default: \code{NULL}, counts table order will be applied
#' @param scale.expression Logic value indicating whether Z-scores should be computed. Default: \code{FALSE} (no scaling).
#' @param x.label.angle Numeric value indicating the rotation angle to use for the x-axis labels. Default: \code{30}.
#'
#' @return A boxplot of class ggplot2.
#'
#' @export expression.boxplot




expression.boxplot =
  function(DEprot.object,
           protein.id,
           which.data = "imputed",
           sample.subset = NULL,
           shape.column = NULL,
           group.by.metadata.column = "column.id",
           group.levels = NULL,
           scale.expression = FALSE,
           x.label.angle = 30) {


    ### Libraries
    require(dplyr)
    require(ggplot2)


    ### Internal functions
    check.matrix =
      function(m){
        warn = "Upon subsetting, no values to show are left."
        if (!is.logical(m)) {
          if (nrow(m) == 0 | ncol(m) == 0) {
            return(return(warning(warn)))
          }
        } else {
          return(return(warning(warn)))
        }
      }



    is.nan_df = function(data.frame) {do.call(cbind, lapply(data.frame, is.nan))}

    ######################################################################################

    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.object))) {
      if (!("DEprot" %in% class(DEprot.object))) {
        return(warning("The input must be an object of class 'DEprot' or 'DEprot.analyses'."))
      }
    }


    ### check grouping column
    if (!is.null(group.by.metadata.column)) {
      if (!(group.by.metadata.column %in% colnames(DEprot.object@metadata))) {
        return(warning(paste0("The 'group.by.metadata.column' is not present in the metadata of the object provided.\n",
                              "Available column IDs: ", paste0(colnames(DEprot.object@metadata), collapse = ", "))))
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



    ### Filter table of counts (samples and protein)
    if (!is.null(sample.subset)) {
      mat.filtered = mat[,which(colnames(mat) %in% sample.subset), drop=F]
    } else {
      mat.filtered = mat
    }
    check.matrix(mat.filtered)


    if (protein.id %in% rownames(mat.filtered)) {
      mat.filtered = mat.filtered[rownames(mat.filtered) == protein.id,,drop=F]
    } else {
      warning(paste0("The protein '", protein.id,"' is not present in the dataset."))
      return(DEprot.object)
    }



    ### reshape table
    exp.tb = as.data.frame(t(mat.filtered))
    colnames(exp.tb)[1] = "expression"
    exp.tb$column.id = rownames(exp.tb)

    # scale/center (z.score)
    if (scale.expression == T) {
      exp.tb = dplyr::mutate(exp.tb, expression = (expression - mean(expression, na.rm = T)) / sd(expression, na.rm = T))
    }



    ### add group column
    if (is.null(group.by.metadata.column)) {group.by.metadata.column = "column.id"}

    if (group.by.metadata.column != "column.id") {
      exp.tb =
        dplyr::left_join(x = exp.tb,
                         y = meta[,c("column.id", group.by.metadata.column)],
                         by = "column.id")

      colnames(exp.tb)[ncol(exp.tb)] = "group"
    } else {
      exp.tb$group = exp.tb$column.id
    }



    ### add replicate column
    if (!is.null(shape.column)) {
      if (shape.column %in% colnames(meta)) {
        exp.tb =
          dplyr::left_join(x = exp.tb,
                           y = meta[,c("column.id", shape.column)],
                           by = "column.id")
        colnames(exp.tb)[ncol(exp.tb)] = "shape"
      } else {
        warning("The 'shape.column' provided is not present in the metadata table.")
      }
    }



    ## add levels to group.column
    if (!is.null(group.levels)) {
      if (all(unique(exp.tb$group) %in% unique(group.levels))) {
        exp.tb = dplyr::mutate(.data = exp.tb, group = factor(group, levels = group.levels))
      } else {
        warning("The 'group.levels' do not include all the groups in the 'group.by.metadata.column'.")
      }
    }




    ### Generate boxplot
    boxplot =
      ggplot(data = exp.tb,
             aes(x = group,
                 y = expression,
                 fill = group,
                 color = group))

    if (scale.expression == T) {
      boxplot = boxplot + geom_hline(yintercept = 0)
    }

    if (group.by.metadata.column != "column.id") {
      boxplot =
        boxplot +
        geom_boxplot(alpha = 0.25,
                     outliers = F,
                     show.legend = F)
    }


    if (!is.null(shape.column)) {
      boxplot =
        boxplot +
        geom_point(aes(shape = factor(shape)),
                   #stroke = NA,
                   size = 3,
                   alpha = 0.5,
                   position = position_jitter(width = 0.15),
                   show.legend = T) +
        guides(shape = guide_legend(title = shape.column))}
    else {
      boxplot =
        boxplot +
        geom_point(stroke = NA,
                   size = 3,
                   alpha = 0.5,
                   position = position_jitter(width = 0.15),
                   show.legend = F)
    }


    boxplot =
      boxplot +
      ggtitle(paste0("**",protein.id,"**")) +
      xlab(NULL) +
      ylab(ifelse(test = scale.expression == T,
                  yes = paste0("centered log~",DEprot.object@log.base,"~(expression)"),
                  no = paste0("log~",DEprot.object@log.base,"~(expression)"))) +
      ggpubr::stat_compare_means(method = "kruskal", show.legend = F) +
      guides(color = "none", fill = "none") +
      theme_classic() +
      theme(axis.title = ggtext::element_markdown(color = "black"),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            axis.text.x = element_text(color = "black", angle = x.label.angle, hjust = ifelse(x.label.angle %in% c(0), yes = 0.5, no = 1)),
            axis.text.y = element_text(color = "black"),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = "black"))



    ### return plot
    return(boxplot)
  } # END function
