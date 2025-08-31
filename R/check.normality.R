#' @title check.normality
#'
#' @description This function performs principal component analyses (PCA).
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param p.threshold Numeric value indicating the p-value threshold to use to define whether a distribution is normal (based on Anderson-Darling normality test). Default: \code{0.05}.
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'imp'. Default: \code{"imputed"}.
#' @param verbose Logical value indicating whether messages should be printed. Default: \code{TRUE}.
#'
#' @return A \code{DEprot.normality} object containing: the \code{norm.statement} indicating whether all the samples are normally distributed or not, the \code{norm.AD.tests} list of Anderson-Darling normality test outputs, a list of \code{qqplots} and \code{densities} profiles for each sample.
#'
#' @import ggplot2
#' @importFrom nortest ad.test
#' @import ggtext
#' @import ggpubr
#'
#' @examples
#' norm <- check.normality(DEprot.object = DEprot::test.toolbox$dpo.imp)
#' patchwork::wrap_plots(norm@qqplots$Sample_A, norm@densities$Sample_A)
#'
#' @export check.normality


check.normality =
  function(DEprot.object,
           p.threshold = 0.05,
           which.data = "imputed",
           verbose = TRUE) {

    # ### Libraries
    # require(ggplot2)


    ### check object
    if (!("DEprot" %in% class(DEprot.object)) & !("DEprot.analyses" %in% class(DEprot.object))) {
      stop("The input must be an object of class 'DEprot'.")
      #return(invisible())
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




    ### Collect norm-test and plots for each sample in the counts
    AD.test = list()
    qqplot = list()
    density = list()

    for (i in 1:ncol(mat)) {
      AD.test[[i]] = nortest::ad.test(mat[,i])
      pval.label = paste0(gsub("e-","\u00D710<sup>-", as.character(formatC(AD.test[[i]]$p.value, format = "e", digits = 2))),"</sup>")

      ### Q-Q plots
      qqplot[[i]] =
        ggpubr::ggqqplot(mat[,i],
                         title = paste0("**",colnames(mat)[i], "**"),
                         stroke = NA,
                         color = ifelse(AD.test[[i]]$p.value < p.threshold, yes = "black", no = "firebrick")) +
        ggtext::geom_richtext(data = data.frame(x = -Inf,
                                                y = +Inf,
                                                label = paste0("Anderson-Darling, *P* = ", pval.label)),
                              mapping = aes(x = x,
                                            y = y,
                                            label = label),
                              fill = NA,
                              label.color = ifelse(AD.test[[i]]$p.value < p.threshold, yes = NA, no = "firebrick"),
                              hjust = -0.1,
                              vjust = 1.5,
                              inherit.aes = FALSE) +
        ylab("Empirical") +
        theme(plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
              aspect.ratio = 1,
              axis.line = element_blank(),
              panel.background = element_rect(fill = NA, linewidth = 1, color = "black"))


      ## Change line color in qqplot
      qqplot[[i]]$layers[[2]]$aes_params$colour = ifelse(AD.test[[i]]$p.value < p.threshold, yes = "steelblue", no = "firebrick")
      qqplot[[i]]$layers[[3]]$aes_params$fill = ifelse(AD.test[[i]]$p.value < p.threshold, yes = "steelblue", no = "firebrick")


      ###### DENSITY
      density[[i]] =
        ggpubr::ggdensity(data = mat[,i],
                          title = paste0("**",colnames(mat)[i], "**"),
                          fill = ifelse(AD.test[[i]]$p.value < p.threshold, yes = "gray", no = "indianred"),
                          color = ifelse(AD.test[[i]]$p.value < p.threshold, yes = "gray50", no = "indianred4")) +
        ggpubr::stat_overlay_normal_density(color = ifelse(AD.test[[i]]$p.value < p.threshold, yes = "steelblue", no = "firebrick4"),
                                            linetype = "dashed",
                                            linewidth = 1.5) +
        xlab("Intensity") +
        ylab("Density") +
        ggtext::geom_richtext(data = data.frame(x = -Inf,
                                                y = +Inf,
                                                label = paste0("Anderson-Darling, *P* = ", pval.label)),
                              mapping = aes(x = x,
                                            y = y,
                                            label = label),
                              fill = NA,
                              label.color = ifelse(AD.test[[i]]$p.value < p.threshold, yes = NA, no = "firebrick"),
                              hjust = -0.1,
                              vjust = 1.5,
                              inherit.aes = FALSE) +
        theme(plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
              aspect.ratio = 1,
              axis.line = element_blank(),
              panel.background = element_rect(fill = NA, linewidth = 1, color = "black"))


    }

    names(AD.test) = colnames(mat)
    names(qqplot) = colnames(mat)
    names(density) = colnames(mat)


    ### Print message if required
    normality = sapply(X = AD.test, FUN = function(x){x$p.value < p.threshold}, USE.NAMES = T)

    if (verbose == T) {
      if (all(normality)) {
        message("All samples display a normal distribution.")

      } else {
        message(paste0("The following samples do not display a normal distribution: ",
                       paste0(names(normality)[isFALSE(normality)], collapse = ", "), "."))
      }
    }


    ### Create output object
    DEprot.normality =
      new(Class = "DEprot.normality",
          norm.statement = all(normality),
          norm.AD.tests = AD.test,
          qqplots = qqplot,
          densities = density,
          p.threshold = p.threshold)


    return(DEprot.normality)
  } #END
