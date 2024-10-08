#' @title impute.counts
#'
#' @description Function that allows for the imputation of missing values using the \href{https://www.rdocumentation.org/packages/missForest/}{missForest} algorithm.
#'
#' @param DEprot.object A \code{DEprot object}, as generated by \link{load.counts}.
#' @param max.iterations Max number of iterations for the missForest algorithm. Default: \code{100}.
#' @param variable.wise.OOBerror Logical value to define whether the OOB error is returned for each variable separately. Default: \code{TRUE}.
#' @param use.normalized.data Logical value indicating whether the imputation should be performed based on the rationalized data. Default: \code{TRUE}.
#' @param overwrite.imputation Logical value to indicate whether, in the case already available, the table of imputed counts should be overwritten. Default: \code{FALSE}.
#' @param cores Number of cores used to run the missForest algorithm. If \code{cores} is 1 (or lower), the imputation will be run in parallel. Two modes are possible and can be defined by the parameter \code{parallel.mode}. Default: \code{1}.
#' @param parallel.mode Define the mode to use for the parallelization, ignored when \code{cores} is more than 1. One among: 'variables', 'forests'. Default: \code{"variables"}. See also the documentation of the \href{https://www.rdocumentation.org/packages/missForest/versions/1.5/topics/missForest}{missForest function}.
#' @param verbose Logical valued indicating whether processing messages should be printed. Default: \code{FALSE}.
#'
#' @seealso \href{https://www.rdocumentation.org/packages/missForest/}{missForest package}.
#'
#' @return A \code{DEprot} object. The boxplot showing the distribution of the protein intensity is remade and added to the slot (\code{boxplot.imputed}). A list with parameters and other info about the imputation is added as well in the \code{imputation} slot.
#'
#' @export impute.counts


impute.counts =
  function(DEprot.object,
           max.iterations = 100,
           variable.wise.OOBerror = TRUE,
           use.normalized.data = TRUE,
           overwrite.imputation = FALSE,
           cores = 1,
           parallel.mode = "variables",
           verbose = FALSE) {

    ### Libraries
    require(dplyr)
    require(ggplot2)
    require(doParallel)
    require(doRNG)

    ### check object
    if (!("DEprot" %in% class(DEprot.object))) {
      warning("The input must be an object of class 'DEprot'.")
      return(DEprot.object)
    }



    ### Check if imputation already available
    if (DEprot.object@imputed == T) {
      if (overwrite.imputation == F) {
        warning(paste0("The 'DEprot' object contains already an imputed table.\n",
                       "If you wish to overwrite the imputation, set the parameter 'overwrite.imputation = TRUE'."))
        return(DEprot.object)
      }
    }


    ### Check if normalized data are available
    if (use.normalized.data == T) {
      if (is.null(DEprot.object@norm.counts)) {
        warning(paste0("You asked to use normalized data for the imputation, but normalized data are not available.\n",
                       "To perform imputation on raw data, set 'use.normalized.data = FALSE'."))
        return(DEprot.object)
      } else {
        cnt = DEprot.object@norm.counts
      }
    } else {
      cnt = DEprot.object@raw.counts
    }



    ### Run missForest algorithm
    # Parallelize if required
    start.time = Sys.time()

    cores = ifelse(test = cores > 1, yes = min(c(cores, nrow(DEprot.object@metadata))), no = 1)
    registerDoParallel(cores = cores)
    #getDoParWorkers()
    registerDoRNG(seed = 1.618)
    DoRNG.check = try(invisible(foreach(i=1:3) %dorng% sqrt(i)))


    if (!("list" %in% class(DoRNG.check)) | cores <= 1) {
      imputed.cnt = missForest::missForest(xmis = t(cnt), maxiter = max.iterations, verbose = verbose, variablewise = variable.wise.OOBerror, parallelize = "no")
    } else {
      if (tolower(parallel.mode) %in% c("variables", "forests")) {
        imputed.cnt = missForest::missForest(xmis = t(cnt), maxiter = max.iterations, verbose = verbose, variablewise = variable.wise.OOBerror, parallelize = tolower(parallel.mode))
      } else {
        warning(paste0("The parallel.mode must be one among: 'variables', 'forests'."))
        return(DEprot.object)
      }
    }

    if (variable.wise.OOBerror == TRUE) {
      names(imputed.cnt$OOBerror) = colnames(t(cnt))
    }

    end.time = Sys.time()
    time.taken = round(end.time - start.time,2)


    ### Replot distributions
    # melt counts table
    melt.cnt =
      suppressMessages(reshape2::melt(as.data.frame(t(imputed.cnt$ximp)))) %>%
      dplyr::mutate(variable = factor(variable, levels = colnames(t(imputed.cnt$ximp))))

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
                  fill = "forestgreen",
                  color = NA) +
      geom_boxplot(data = melt.cnt,
                   mapping = aes(x = variable,
                                 y = value,
                                 group = variable),
                   fill = "white",
                   color = "darkgreen",
                   width = 0.15,
                   outlier.color = "black",
                   outlier.stroke = NA,
                   outlier.size = 2,
                   outlier.alpha = 0.25) +
      geom_line(data = data.frame(cnt.stats),
                mapping = aes(x = variable,
                              y = max,
                              group = 1),
                color = "indianred",
                linetype = 2,
                inherit.aes = F) +
      geom_line(data = data.frame(cnt.stats),
                mapping = aes(x = variable,
                              y = min,
                              group = 1),
                color = "steelblue",
                linetype = 2,
                inherit.aes = F) +
      ylab(ifelse(is.null(DEprot.object@log.base),
                  yes = "Intensity",
                  no = paste0(ifelse(DEprot.object@log.base == exp(1),
                                     yes = "ln", no = paste0("log~",DEprot.object@log.base,"~")),
                              "(Intensity)"))) +
      ggtitle(label = "**Imputed data**", subtitle = "*missForest*") +
      xlab("Sample") +
      theme_classic() +
      theme(axis.text.y = element_text(color = "black"),
            axis.text.x = element_text(color = "black", hjust = 1, angle = 30),
            axis.title = ggtext::element_markdown(color = "black"),
            axis.ticks.y = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5),
            aspect.ratio = 10/ncol(t(imputed.cnt$ximp)))


    ### Update object with new counts, imputation method and boxplot
    DEprot.object@imputed = T
    DEprot.object@imputation = list(method = "missForest",
                                    max.iterations = max.iterations,
                                    OOBerror = imputed.cnt$OOBerror,
                                    parallelization.mode = ifelse(cores <=1, yes = "none", no = parallel.mode),
                                    cores = cores,
                                    processing.time = paste(gsub("Time difference of ", "",as.character(time.taken)), attributes(time.taken)$units))

    DEprot.object@imputed.counts = t(imputed.cnt$ximp)
    DEprot.object@boxplot.imputed = boxplot


    ### Return updated object
    return(DEprot.object)

  } # END function
