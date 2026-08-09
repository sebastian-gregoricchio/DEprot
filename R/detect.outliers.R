#' @title detect.outliers
#'
#' @description Automatic, threshold-based flagging of low-quality samples. Three orthogonal per-sample
#' quality metrics are combined: (i) the median correlation of each sample against all the others (or against
#' the replicates of its own group), (ii) a robust squared Mahalanobis distance computed in the space of the
#' first principal components, and (iii) the fraction of missing values. Each metric is compared to its own
#' threshold and a sample is called an outlier when at least \code{min.flags} metrics are triggered.
#' The vector of flagged samples can be passed directly to \link{filter.samples}.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param which.data String indicating which type of counts should be used for the correlation and the PCA. One among: 'raw', 'normalized', 'norm', 'randomized', 'random', 'imputed', 'imp'. Default: \code{"imputed"}.
#' @param sample.subset String vector indicating the column names (samples) to keep in the counts table (the 'column.id' in the metadata table). Default: \code{NULL} (no subsetting).
#' @param group.column String indicating the column of the metadata table defining the sample groups. When provided, the median correlation of a sample is computed only against the other samples of the same group (replicate correlation). Default: \code{NULL} (all the samples are used).
#' @param correlation.method String indicating the correlation method to use. Possible options: 'pearson', 'spearman', 'kendall'. Default: \code{"pearson"}.
#' @param missingness.data String indicating which counts should be used to quantify the missing values. One among: 'auto', 'raw', 'normalized', 'norm', 'randomized', 'random'. With \code{"auto"} the first available table among raw, normalized and randomized counts is used. Default: \code{"auto"}.
#' @param n.PCs Integer indicating the number of principal components used to compute the Mahalanobis distances. Default: \code{3}.
#' @param center.data Logical value indicating whether the data should be centered for the PCA. Default: \code{TRUE}.
#' @param scale.data Logical value indicating whether the data should be scaled for the PCA. Default: \code{TRUE}.
#' @param correlation.z.th Numeric (negative) value indicating the robust Z-score below which the median correlation of a sample is considered too low. Default: \code{-2.5}.
#' @param correlation.min Numeric value indicating an absolute floor for the median correlation: samples below this value are flagged independently of their Z-score. Default: \code{NULL} (not applied).
#' @param mahalanobis.padj.th Numeric value indicating the adjusted p-value threshold applied to the chi-square test of the Mahalanobis distances. Default: \code{0.05}.
#' @param missingness.z.th Numeric (positive) value indicating the robust Z-score above which the missing rate of a sample is considered too high. Default: \code{2.5}.
#' @param missingness.max Numeric value between 0-1 indicating an absolute cap for the missing rate: samples above this value are flagged independently of their Z-score. Default: \code{NULL} (not applied).
#' @param padj.method String indicating the multiple-testing correction applied to the chi-square p-values. Any method supported by \link[stats]{p.adjust}. Default: \code{"BH"}.
#' @param min.flags Integer indicating the minimum number of metrics that must be triggered to call a sample an outlier. Automatically capped to the number of metrics that could be computed. Default: \code{2}.
#' @param verbose Logical value indicating whether messages should be printed. Default: \code{TRUE}.
#'
#' @details The three metrics are deliberately independent: the correlation captures samples whose global
#' expression profile does not resemble the others, the Mahalanobis distance captures samples lying far from
#' the bulk of the data in the reduced PC space, and the missing rate captures samples with a poor identification
#' depth. Because the principal components returned by \link{perform.PCA} are orthogonal by construction, the
#' Mahalanobis distance is computed with a diagonal covariance matrix estimated robustly (median and MAD of
#' each PC): this avoids the singularity that a full covariance matrix would generate whenever the number of
#' samples is close to, or smaller than, the number of components. The resulting squared distances are compared
#' to a chi-square distribution with as many degrees of freedom as the components effectively used.\cr
#' The correlations are computed with \code{use = "pairwise.complete.obs"}, so that unimputed tables do not
#' collapse to the few proteins quantified in all the samples. All the Z-scores are robust
#' (\code{(x - median(x)) / mad(x)}) and are always computed across all the samples, also when
#' \code{group.column} is provided.
#'
#' @return A \code{DEprot.outliers} object, containing the per-sample metric table (\code{metrics}), the vector
#' of the flagged samples (\code{outliers}) and the diagnostic plots (\code{plot}, \code{plot.list}).
#'
#' @name detect.outliers
#'
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#' @importFrom patchwork wrap_plots
#' @importFrom methods new slot
#' @importFrom stats cor mad median p.adjust pchisq sd
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' outliers <- detect.outliers(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                             which.data = "imputed",
#'                             group.column = "condition")
#'
#' outliers
#' plot(outliers)
#'
#' @export detect.outliers


detect.outliers =
  function(DEprot.object,
           which.data = "imputed",
           sample.subset = NULL,
           group.column = NULL,
           correlation.method = "pearson",
           missingness.data = "auto",
           n.PCs = 3,
           center.data = TRUE,
           scale.data = TRUE,
           correlation.z.th = -2.5,
           correlation.min = NULL,
           mahalanobis.padj.th = 0.05,
           missingness.z.th = 2.5,
           missingness.max = NULL,
           padj.method = "BH",
           min.flags = 2,
           verbose = TRUE) {

    ### Functions
    say = function(...) {if (isTRUE(verbose)) {message(...)}}

    robust.z = # robust Z-score; returns NAs when the sample is not dispersed enough to define a scale
      function(x) {
        center = median(x, na.rm = TRUE)
        scale = mad(x, na.rm = TRUE)
        if (is.na(scale) | scale == 0) {scale = stats::sd(x, na.rm = TRUE)}
        if (is.na(scale) | scale == 0) {return(rep(NA, length(x)))}
        return((x - center) / scale)
      }

    #############################################


    ### check object
    if (!("DEprot" %in% class(DEprot.object)) & !("DEprot.analyses" %in% class(DEprot.object))) {
      stop("The input must be an object of class 'DEprot'.")
    }

    ### check the correlation method
    if (!(tolower(correlation.method) %in% c("pearson", "spearman", "kendall"))) {
      stop("The 'correlation.method' must be one among: 'pearson', 'spearman', 'kendall'.")
    }
    correlation.method = tolower(correlation.method)


    ### Check and extract table
    if (tolower(which.data) == "raw") {
      if (!is.null(DEprot.object@raw.counts)) {
        mat = DEprot.object@raw.counts
        data.used = "raw"
      } else {
        stop(paste0("Use of RAW counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
      }
    } else if (tolower(which.data) %in% c("norm", "normalized", "normal")) {
      if (!is.null(DEprot.object@norm.counts)) {
        mat = DEprot.object@norm.counts
        data.used = "normalized"
      } else {
        stop(paste0("Use of NORMALIZED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
      }
    } else if (tolower(which.data) %in% c("imputed", "imp", "impute")) {
      if (!is.null(DEprot.object@imputed.counts)) {
        mat = DEprot.object@imputed.counts
        data.used = "imputed"
      } else {
        stop(paste0("Use of IMPUTED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
      }
    } else if (tolower(which.data) %in% c("randomized", "random")) {
      if (!is.null(DEprot.object@random.counts)) {
        mat = DEprot.object@random.counts
        data.used = "randomized"
      } else {
        stop(paste0("Use of RANDOMIZED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
      }
    } else {
      stop(paste0("The 'which.data' value is not recognized.\n",
                  "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
    }


    ### subset table and metadata
    if (!is.null(sample.subset)) {
      if (!all(sample.subset %in% colnames(mat))) {
        stop(paste0("The following samples are not available in the counts table: ",
                    paste0(sample.subset[!(sample.subset %in% colnames(mat))], collapse = ", "), "."))
      }
      mat = mat[,which(colnames(mat) %in% sample.subset), drop = FALSE]
    }

    meta = dplyr::filter(DEprot.object@metadata, column.id %in% colnames(mat))
    meta = meta[match(colnames(mat), meta$column.id),]
    samples = colnames(mat)

    if (length(samples) < 4) {
      stop("At least 4 samples are required to compute the outlier statistics.")
    }

    ### check the group column
    if (!is.null(group.column)) {
      if (!(group.column %in% colnames(meta))) {
        stop(paste0("The column '", group.column, "' is not available in the metadata table."))
      }
      groups = as.character(meta[,group.column])
    } else {
      groups = rep("all.samples", length(samples))
    }

    mat[is.nan(mat)] = NA



    ##########################
    ### METRIC 1: correlation
    ##########################
    corr.mat = cor(mat, method = correlation.method, use = "pairwise.complete.obs")
    diag(corr.mat) = NA

    median.correlation =
      sapply(X = samples,
             FUN = function(s) {
               partners = samples[groups == groups[which(samples == s)] & samples != s]
               if (length(partners) == 0) {partners = samples[samples != s]}
               return(median(corr.mat[s, partners], na.rm = TRUE))
             },
             USE.NAMES = TRUE)

    correlation.z = robust.z(median.correlation)

    flag.correlation = (!is.na(correlation.z) & correlation.z <= correlation.z.th)
    if (!is.null(correlation.min)) {
      flag.correlation = flag.correlation | (median.correlation < correlation.min)
    }

    if (all(is.na(correlation.z)) & is.null(correlation.min)) {
      say("All the samples display the same median correlation: the correlation metric will not be used.")
      correlation.available = FALSE
    } else {
      correlation.available = TRUE
    }



    ###########################
    ### METRIC 2: PC-distances
    ###########################
    PCA = DEprot::perform.PCA(DEprot.object = DEprot.object,
                              sample.subset = samples,
                              which.data = data.used,
                              center.data = center.data,
                              scale.data = scale.data)

    PC.tb = PCA@PCs[match(samples, PCA@PCs$column.id),]
    PC.available = colnames(PC.tb)[grepl("^PC[0-9]+$", colnames(PC.tb))]
    n.PCs.used = min(c(round(n.PCs,0), length(PC.available), length(samples)-1), na.rm = TRUE)

    if (n.PCs.used < round(n.PCs,0)) {
      say(paste0("Only ", n.PCs.used, " principal components could be used (", round(n.PCs,0), " were required)."))
    }

    PC.scores = as.matrix(PC.tb[,paste0("PC", 1:n.PCs.used), drop = FALSE])

    # robust standardization of each PC (the PCs are orthogonal: a diagonal covariance is sufficient)
    PC.z = apply(X = PC.scores, MARGIN = 2, FUN = robust.z)
    PC.z = matrix(PC.z, nrow = nrow(PC.scores), dimnames = list(samples, colnames(PC.scores)))
    PC.z = PC.z[,colSums(is.na(PC.z)) == 0, drop = FALSE]
    df.chisq = ncol(PC.z)

    if (df.chisq == 0) {
      say("None of the principal components displays a measurable dispersion: the Mahalanobis metric will not be used.")
      mahalanobis.distance = rep(NA, length(samples))
      mahalanobis.pvalue = rep(NA, length(samples))
      mahalanobis.padj = rep(NA, length(samples))
      flag.mahalanobis = rep(FALSE, length(samples))
      mahalanobis.available = FALSE
    } else {
      mahalanobis.distance = rowSums(PC.z^2)
      mahalanobis.pvalue = pchisq(q = mahalanobis.distance, df = df.chisq, lower.tail = FALSE)
      mahalanobis.padj = p.adjust(mahalanobis.pvalue, method = padj.method)
      flag.mahalanobis = (mahalanobis.padj < mahalanobis.padj.th)
      mahalanobis.available = TRUE
    }



    ##########################
    ### METRIC 3: missingness
    ##########################
    missingness.priority = c(raw = "raw.counts", normalized = "norm.counts", randomized = "random.counts")

    missingness.available.slots = missingness.priority[sapply(missingness.priority, function(s){!is.null(methods::slot(DEprot.object, s))})]

    if (tolower(missingness.data) == "auto") {
      missingness.slot = missingness.available.slots[1]
    } else if (tolower(missingness.data) == "raw") {
      missingness.slot = missingness.available.slots["raw"]
    } else if (tolower(missingness.data) %in% c("norm", "normalized", "normal")) {
      missingness.slot = missingness.available.slots["normalized"]
    } else if (tolower(missingness.data) %in% c("randomized", "random")) {
      missingness.slot = missingness.available.slots["randomized"]
    } else {
      stop(paste0("The 'missingness.data' value is not recognized.\n",
                  "       Please indicated a count type among 'auto', 'raw', 'normalized', 'randomized'."))
    }

    if (length(missingness.slot) == 0 | all(is.na(missingness.slot))) {
      say(ifelse(tolower(missingness.data) == "auto",
                 yes = "No unimputed counts table is available: the missingness metric will not be used.",
                 no = paste0("The ", tolower(missingness.data), " counts are not available: the missingness metric will not be used.")))
      missing.rate = rep(NA, length(samples))
      missingness.z = rep(NA, length(samples))
      flag.missingness = rep(FALSE, length(samples))
      missingness.available = FALSE
      missingness.data.used = NA
    } else {
      missingness.data.used = names(missingness.slot)
      miss.mat = methods::slot(DEprot.object, unname(missingness.slot))

      if (!all(samples %in% colnames(miss.mat))) {
        stop(paste0("The ", missingness.data.used, " counts do not contain all the samples analyzed. ",
                    "Use 'missingness.data' to indicate another table."))
      }

      miss.mat = miss.mat[,samples, drop = FALSE]
      miss.mat[is.nan(miss.mat)] = NA
      missing.rate = colSums(is.na(miss.mat)) / nrow(miss.mat)
      missingness.z = robust.z(missing.rate)

      flag.missingness = (!is.na(missingness.z) & missingness.z >= missingness.z.th)
      if (!is.null(missingness.max)) {
        flag.missingness = flag.missingness | (missing.rate > missingness.max)
      }

      if (all(missing.rate == 0)) {
        say(paste0("The ", missingness.data.used, " counts do not contain missing values: the missingness metric will not be used."))
        missingness.available = FALSE
      } else if (all(is.na(missingness.z)) & is.null(missingness.max)) {
        say("All the samples display the same missing rate: the missingness metric will not be used.")
        missingness.available = FALSE
      } else {
        missingness.available = TRUE
      }
    }



    ####################
    ### Combine metrics
    ####################
    metrics.available = c(correlation = correlation.available,
                          mahalanobis = mahalanobis.available,
                          missingness = missingness.available)

    if (sum(metrics.available) == 0) {
      stop("None of the three metrics could be computed on this data set.")
    }

    min.flags.used = min(c(round(min.flags,0), sum(metrics.available)), na.rm = TRUE)
    if (min.flags.used != round(min.flags,0)) {
      say(paste0("Only ", sum(metrics.available), " metric(s) could be computed: 'min.flags' has been set to ", min.flags.used, "."))
    }

    if (correlation.available == FALSE) {flag.correlation = rep(FALSE, length(samples))}
    if (mahalanobis.available == FALSE) {flag.mahalanobis = rep(FALSE, length(samples))}
    if (missingness.available == FALSE) {flag.missingness = rep(FALSE, length(samples))}

    n.flags = as.numeric(flag.correlation) + as.numeric(flag.mahalanobis) + as.numeric(flag.missingness)

    metrics =
      data.frame(column.id = samples,
                 group = groups,
                 median.correlation = as.numeric(median.correlation),
                 correlation.z = as.numeric(correlation.z),
                 flag.correlation = flag.correlation,
                 mahalanobis.distance = as.numeric(mahalanobis.distance),
                 mahalanobis.pvalue = as.numeric(mahalanobis.pvalue),
                 mahalanobis.padj = as.numeric(mahalanobis.padj),
                 flag.mahalanobis = flag.mahalanobis,
                 missing.rate = as.numeric(missing.rate),
                 missingness.z = as.numeric(missingness.z),
                 flag.missingness = flag.missingness,
                 n.flags = n.flags,
                 outlier = (n.flags >= min.flags.used),
                 row.names = NULL,
                 stringsAsFactors = FALSE) %>%
      dplyr::left_join(meta, by = "column.id") %>%
      dplyr::arrange(dplyr::desc(n.flags), dplyr::desc(mahalanobis.distance))

    outliers = metrics$column.id[metrics$outlier]

    say(paste0(length(outliers), " sample(s) flagged as outlier out of ", length(samples),
               ifelse(length(outliers) > 0,
                      yes = paste0(": ", paste0(outliers, collapse = ", "), "."),
                      no = ".")))



    ##################
    ### Metric plots
    ##################
    plot.tb =
      metrics %>%
      dplyr::mutate(column.id = factor(column.id, levels = metrics$column.id),
                    status = factor(ifelse(outlier == TRUE, yes = "outlier", no = "retained"),
                                    levels = c("retained", "outlier")))

    metric.plot =
      function(y.values, flags, y.label, title, thresholds = NULL) {
        tb = dplyr::mutate(plot.tb, metric.value = y.values, flag = flags)
        tb = tb[!is.na(tb$metric.value),]

        p =
          ggplot(data = tb,
                 aes(x = column.id,
                     y = metric.value)) +
          geom_segment(mapping = aes(xend = column.id, yend = -Inf),
                       color = "gray70",
                       linewidth = 0.3)

        if (!is.null(thresholds)) {
          p = p + geom_hline(yintercept = thresholds[!is.na(thresholds)],
                             linetype = 2,
                             linewidth = 0.3,
                             color = "indianred")
        }

        p +
          geom_point(mapping = aes(fill = flag),
                     shape = 21,
                     stroke = NA,
                     size = 2.5) +
          scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "indianred"),
                            guide = "none") +
          ylab(y.label) +
          xlab(NULL) +
          ggtitle(label = title) +
          theme_classic() +
          theme(axis.text = element_text(color = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                panel.background = element_blank(),
                plot.background = element_blank(),
                panel.grid.major.y = element_line(color = "gray", linewidth = 0.1),
                plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
                aspect.ratio = 0.35)
      }

    plot.list = list()

    if (correlation.available == TRUE) {
      plot.list[["correlation"]] =
        metric.plot(y.values = plot.tb$median.correlation,
                    flags = plot.tb$flag.correlation,
                    y.label = paste0("Median ", correlation.method, " correlation"),
                    title = paste0("**Median inter-sample correlation**",
                                   ifelse(!is.null(group.column), yes = "<br>*within group*", no = "")),
                    thresholds = c(median(plot.tb$median.correlation, na.rm = TRUE) +
                                     correlation.z.th * mad(plot.tb$median.correlation, na.rm = TRUE),
                                   correlation.min))
    }

    if (mahalanobis.available == TRUE) {
      plot.list[["mahalanobis"]] =
        metric.plot(y.values = -log10(plot.tb$mahalanobis.padj),
                    flags = plot.tb$flag.mahalanobis,
                    y.label = paste0("-log~10~(", padj.method, " p-value)"),
                    title = paste0("**PC-space Mahalanobis distance**<br>*", df.chisq, " PCs*"),
                    thresholds = -log10(mahalanobis.padj.th))

      plot.list[["mahalanobis"]] = plot.list[["mahalanobis"]] + theme(axis.title.y = ggtext::element_markdown())
    }

    if (missingness.available == TRUE) {
      plot.list[["missingness"]] =
        metric.plot(y.values = plot.tb$missing.rate * 100,
                    flags = plot.tb$flag.missingness,
                    y.label = "Missing values (%)",
                    title = paste0("**Missing rate**<br>*", missingness.data.used, " counts*"),
                    thresholds = c((median(plot.tb$missing.rate, na.rm = TRUE) +
                                      missingness.z.th * mad(plot.tb$missing.rate, na.rm = TRUE)) * 100,
                                   missingness.max * 100))
    }

    combined.plot = patchwork::wrap_plots(plot.list, ncol = 1)



    ################
    ### Export data
    ################
    outlier.object =
      methods::new(Class = "DEprot.outliers",
                   metrics = metrics,
                   outliers = outliers,
                   sample.subset = samples,
                   data.used = data.used,
                   correlation.method = correlation.method,
                   correlation.matrix = corr.mat,
                   PCA = PCA,
                   missingness.data.used = missingness.data.used,
                   group.column = group.column,
                   metrics.available = metrics.available,
                   parameters = list(n.PCs = df.chisq,
                                     correlation.z.th = correlation.z.th,
                                     correlation.min = correlation.min,
                                     mahalanobis.padj.th = mahalanobis.padj.th,
                                     padj.method = padj.method,
                                     missingness.z.th = missingness.z.th,
                                     missingness.max = missingness.max,
                                     min.flags = min.flags.used),
                   plot = combined.plot,
                   plot.list = plot.list)

    return(outlier.object)
  } # END function
