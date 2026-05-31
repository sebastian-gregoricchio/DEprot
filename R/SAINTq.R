#' @title SAINTq
#'
#' @description
#' Native R implementation of the \strong{SAINTq} interaction score
#' (Teo et al., \emph{Proteomics} 2016, doi:10.1002/pmic.201500499) for
#' protein-level LFQ data stored in a \code{DEprot} object. No external
#' \code{saintq} binary is required.
#'
#' At the protein level SAINTq is identical to the intensity module of
#' SAINTexpress / SAINT-MS1 (Teo et al., \emph{J Proteomics} 2014; Choi et al.,
#' \emph{J Proteome Res} 2012). This function reproduces that published model:
#' a \emph{semi-supervised two-component mixture} in which, for every prey, the
#' background ("false") component is estimated from the negative-control runs and
#' the "true" component is placed at a fixed fold-change above that background
#' (the \eqn{\theta_T = 5\,\theta_F} rule of SAINTexpress). For each bait
#' replicate the posterior probability of a true interaction is obtained by
#' Bayes' theorem; \code{AvgP} is the average over replicates, \code{MaxP} the
#' maximum, the fold change is mean(bait)/mean(control), and \code{BFDR} is the
#' Bayesian FDR computed directly from the probabilities.
#'
#' @param DEprot.object A \code{DEprot} (or \code{DEprot.analyses}) object.
#' @param metadata.column String indicating the ID of the column of the metadata that holds the group label of each sample.
#' @param control String indicating the group name in \code{metadata.column} that represent the negative controls (e.g. IgG, empty vector).
#' @param bait String (or character vector) indicating the group name(s) in \code{metadata.column} to score as bait. Each distinct value is scored separately against the control. Default: \code{NULL}, scores every 'non-control' group.
#' @param fold Numeric value > 0. Fold change separating the true component from the background, i.e. \code{mean_True = mean_False + log2(fold)}. Default: \code{5} (default in SAINTexpress).
#' @param prior.pi Either \code{NULL} (estimate the prior proportion of true interactions by Expectation-Maximization, default) or a fixed value between 0 and 1. Default: \code{NULL}.
#' @param sd.scale Numeric value > 0. Standard deviation of the true component relative to the background (\code{sigma_T = sd.scale * sigma_F}). Default: \code{1}.
#' @param min.sd  Numeric value >= 0. Lower bound on the per-prey background SD, to avoid over-confident scores from preys with near-constant controls. Default: \code{0.1}.
#' @param control.background Optional numeric (log2 scale). Background level used for preys never detected in the controls. Default: \code{NULL}, uses the lowest observed control value.
#' @param best.n.rep Optional integer. If set, \code{AvgP} averages only the top \code{best.n.rep} replicate probabilities of each bait (equivalent to SAINTexpress's "best R replicates"). \code{NULL} uses all replicates.
#' @param which.data String indicating which data should be used. One among "imputed", "randomized", "normalized", "raw". Default: \code{"imputed"}.
#' @param viridis.palette Character string indicating the color map option to use (passed to \code{viridis}). Eight options are available: "magma" (or "A"), "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D"), "cividis" (or "E"), "rocket" (or "F"), "mako" (or "G"), "turbo" (or "H"). Default: \code{"magma"}.
#' @param viridis.direction Sets the order of colors in the scale. If \code{1} colors are as output by \code{viridis_pal}, if -1 the order of colors is reversed. Default: \code{-1}.
#' @param viridis.begin The (corrected) hue in [0,1] at which the color map begins. Default: \code{0.15}.
#' @param viridis.end The (corrected) hue in [0,1] at which the color map ends. Default: \code{0.9}.
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return An object of class \code{DEprot.SAINTq}.
#'
#' @section Faithfulness:
#' This is an independent re-implementation of the \emph{published} model, not a
#' re-wrapping of the C++ program. It returns scores with the same definitions
#' and behavior, but exact numbers can differ from the official \code{SAINTq}
#' binary because of low-level differences (parameter estimation, missing-value
#' and control-compression heuristics). For publication-critical work, validate a
#' subset against the official tool.
#'
#' @references
#' Teo G. \emph{et al.} (2016) SAINTq. \emph{Proteomics} 16:2238-2245.
#' Teo G. \emph{et al.} (2014) SAINTexpress. \emph{J Proteomics} 100:37-43.
#' Choi H. \emph{et al.} (2012) SAINT-MS1. \emph{J Proteome Res} 11:2619-2624.
#'
#' @import ggplot2
#' @importFrom stats dnorm
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom ggpubr theme_pubr
#' @importFrom ggtext element_markdown
#' @importFrom viridis scale_color_viridis
#'
#'
#' @examples
#' \dontrun{
#'   sq <- SAINTq(DEprot.object = dpo,
#'                metadata.column = "Target",
#'                control = "IgG",
#'                bait = "BRD4")
#'
#'
#'   sq_BRD4_HDAC1 <- SAINTq(DEprot.object = dpo,
#'                           metadata.column = "Target",
#'                           control = "IgG",
#'                           bait = c("BRD4", "HDAC1"))
#'
#'
#'   sq_all <- SAINTq(DEprot.object = dpo,
#'                    metadata.column = "Target",
#'                    control = "IgG",
#'                    bait = NULL)  #all non-IgG samples in 'Target'
#' }
#'
#' @export


SAINTq =
  function(DEprot.object,
           metadata.column,
           control,
           bait = NULL,
           fold = 5,
           prior.pi = NULL,
           sd.scale = 1,
           min.sd = 0.1,
           control.background = NULL,
           best.n.rep = NULL,
           which.data = "imputed",
           viridis.palette = "magma",
           viridis.direction = -1,
           viridis.begin = 0.15,
           viridis.end = 0.9,
           verbose = TRUE) {

    ## ---------------------  check parameters  ---------------------
    ### Check parameters
    if (missing(DEprot.object)) {stop("Provide a `DEprot.object` of class 'DEprot' or 'DEprot.analyses'.", call. = FALSE)}
    if (missing(metadata.column)) {stop("Provide a character indicating the ID of the `metadata.column` that holds the group label of each sample.", call. = FALSE)}
    if (missing(control) || length(control) > 1) {stop("Provide one 'control' group.", call. = FALSE)}

    ### check object and parameters
    if (!("DEprot" %in% class(DEprot.object))) {
      if (!("DEprot.analyses" %in% class(DEprot.object))) {
        stop("The `DEprot.object` must be an object of class 'DEprot' or 'DEprot.analyses'.")
      }
    }


    ### Check and extract table
    if (tolower(which.data) == "raw") {
      if (!is.null(DEprot.object@raw.counts)) {
        mat = DEprot.object@raw.counts
        data.used = "raw"
      } else {
        stop(paste0("Use of RAW counts was required, but not available.\n",
                    "Please indicated a count type among 'raw', 'normalized', 'randomized, 'imputed', using the option 'which.data'."))
      }
    } else if (tolower(which.data) %in% c("norm", "normalized", "normalised", "normal")) {
      if (!is.null(DEprot.object@norm.counts)) {
        mat = DEprot.object@norm.counts
        data.used = "normalized"
      } else {
        stop(paste0("Use of NORMALIZED counts was required, but not available.\n",
                    "Please indicated a count type among 'raw', 'normalized', 'randomized, 'imputed', using the option 'which.data'."))
      }
    } else if (tolower(which.data) %in% c("imputed", "imp", "impute")) {
      if (!is.null(DEprot.object@imputed.counts)) {
        mat = DEprot.object@imputed.counts
        data.used = "imputed"
      } else {
        stop(paste0("Use of IMPUTED counts was required, but not available.\n",
                    "Please indicated a count type among 'raw', 'normalized', 'randomized, 'imputed', using the option 'which.data'."))
      }
    } else if (tolower(which.data) %in% c("randomized", "randomised", "random")) {
      if (!is.null(DEprot.object@random.counts)) {
        mat = DEprot.object@random.counts
        data.used = "randomized"
      } else {
        stop(paste0("Use of RANDOMIZED counts was required, but not available.\n",
                    "Please indicated a count type among 'raw', 'normalized', 'randomized, 'imputed', using the option 'which.data'."))
      }
    } else {
      stop(paste0("The 'which.data' value is not recognized.\n",
                  "Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
    }


    ## Convert table to log2
    if (!is.numeric(DEprot.object@log.base)) {
      warning("The log.base is not numeric, linear counts are assumed. Counts matrix will be converted to log2+1 values to analyze the data.")
      mat.log2 = log2(mat + 1)
    } else if (as.numeric(DEprot.object@log.base) != 2) {
      warning("The log.base is not 2, counts will be converted to log2 values to analyze the data.")
      mat.log2 = log2(DEprot.object@log.base^mat)
    } else {
      mat.log2 = mat
    }

    # ------------------------------------------------------------------------------------


    ## Map samples to groups
    meta = as.data.frame(DEprot.object@metadata, stringsAsFactors = FALSE)

    if (!metadata.column %in% colnames(meta)) {
      stop("'metadata.column' = '", metadata.column, "' is not in @metadata. Available: ",
           paste(colnames(meta), collapse = ", "), call. = FALSE)
    }


    # Get control samples
    ctrl.cols = which(colnames(mat.log2) %in% meta[which(meta[,metadata.column] == control), "column.id"])
    if (length(ctrl.cols) < 1) stop("No samples match the 'control'.", call. = FALSE)

    # Get bait samples
    if (is.null(bait)) {bait = unique(meta[which(meta[,metadata.column] != control), metadata.column])}
    bait = intersect(bait, unique(meta[,metadata.column]))

    if (length(bait) < 1) {stop("None of the requested 'bait' groups were found.", call. = FALSE)}



    ## per-prey background ('false') component from the controls -------------------------
    ctrl.mat = mat.log2[, ctrl.cols, drop = FALSE]

    n.na.ctrl = rowSums(!is.na(ctrl.mat))

    mean.control = rowMeans(ctrl.mat, na.rm = TRUE)
    mean.control[!is.finite(mean.control)] = NA

    sd.control = apply(ctrl.mat, 1, function(x) {stats::sd(x, na.rm = TRUE)})
    sd.control[n.na.ctrl < 2] = NA

    sigma.global = stats::median(sd.control[is.finite(sd.control) & sd.control > 0], na.rm = TRUE)
    if (!is.finite(sigma.global) || sigma.global <= 0) {
      sigma.global = stats::sd(as.vector(mat.log2), na.rm = TRUE)
      if (!is.finite(sigma.global) || sigma.global <= 0) sigma.global = 1
    }

    bg = ifelse(test = is.null(control.background),
                yes = min(ctrl.mat, na.rm = TRUE),
                no = control.background)
    if (!is.finite(bg)) {bg = min(mat.log2, na.rm = TRUE)}
    if (!is.finite(bg)) {bg = 0}

    mean.control.false.use = ifelse(n.na.ctrl >= 1 & is.finite(mean.control), mean.control, bg)
    sd.control.false.use = ifelse(is.finite(sd.control) & sd.control > 0, sd.control, sigma.global)
    sd.control.false.use = pmax(sd.control.false.use, min.sd)



    ## True component (fixed offset) -----------------------------------------------------
    delta = log2(fold)
    mean.control.true.use = mean.control.false.use + delta
    sd.control.true.use = pmax(sd.control.false.use * sd.scale, min.sd)



    ## Estimate prior proportion of true interactions (EM) -------------------------------
    bait.cols = which(colnames(mat.log2) %in% meta[which(meta[,metadata.column] %in% bait), "column.id"])
    Xall = as.vector(mat.log2[, bait.cols, drop = FALSE])
    ridx = rep(seq_len(nrow(mat.log2)), times = length(bait.cols))
    keep = !is.na(Xall)
    Xall = Xall[keep]; ridx = ridx[keep]
    fF.all = stats::dnorm(x = Xall, mean = mean.control.false.use[ridx], sd = sd.control.false.use[ridx])
    fT.all = stats::dnorm(x = Xall, mean = mean.control.true.use[ridx], sd = sd.control.true.use[ridx])

    if (is.null(prior.pi)) {
      pi.t = 0.1
      for (it in seq_len(200L)) {
        den = pi.t * fT.all + (1 - pi.t) * fF.all
        g = ifelse(den > 0, pi.t * fT.all / den, 0)
        new = mean(g)
        if (abs(new - pi.t) < 1e-8) { pi.t = new; break }
        pi.t = new
      }
      pi.t = min(max(pi.t, 1e-6), 1 - 1e-6)
    } else {
      if (prior.pi <= 0 || prior.pi >= 1) stop("'prior.pi' must be in a range of (0, 1).", call. = FALSE)
      pi.t = prior.pi
    }

    if (isTRUE(verbose)) {
      message(sprintf("Prior P(true) = %.4f; true/false fold = %g (Delta log2 = %.3f).",
                      pi.t, fold, delta))}



    ## Score every bait against the pooled controls --------------------------------------
    prob.fun = function(x, i) {
      fF = stats::dnorm(x, mean.control.false.use[i], sd.control.false.use[i])
      fT = stats::dnorm(x, mean.control.true.use[i], sd.control.true.use[i])
      den = pi.t * fT + (1 - pi.t) * fF
      ifelse(den > 0, pi.t * fT / den, 0)
    }

    res.list =
      lapply(unique(bait),
             function(b) {
               bc = which(colnames(mat.log2) %in% meta[which(meta[,metadata.column] == b), "column.id"])
               Xb = mat.log2[, bc, drop = FALSE]
               detected = rowSums(!is.na(Xb)) > 0

               if (!any(detected)) {return(NULL)}
               idx = which(detected)

               prob = matrix(NA, length(idx), length(bc))

               for (r in seq_along(bc)) {
                 prob[, r] = prob.fun(Xb[idx, r], idx)}

               avgP =
                 vapply(seq_along(idx),
                        function(k) {
                          p = prob[k, ]; p = p[!is.na(p)]
                          if (!length(p)) {return(0)}
                          if (!is.null(best.n.rep)) {p = sort(p, decreasing = TRUE)[seq_len(min(best.n.rep, length(p)))]}
                          return(mean(p))},
                        numeric(1))

               maxP =
                 vapply(seq_along(idx),
                        function(k) {
                          p = prob[k, ]
                          p = p[!is.na(p)]
                          return(ifelse(test = length(p) != 0, yes = max(p), no = 0))
                        },
                        numeric(1))

               avg.bait = rowMeans(Xb[idx, , drop = FALSE], na.rm = TRUE)
               log2fc = avg.bait - mean.control.false.use[idx]

               return(data.frame(Bait = b,
                                 Control = control,
                                 Prey = rownames(mat.log2)[idx],
                                 n.rep = rowSums(!is.na(Xb[idx, , drop = FALSE])),
                                 AvgP = avgP, MaxP = maxP,
                                 log2.FoldChange_bait.vs.control = log2fc,
                                 avg.bait = avg.bait,
                                 avg.ctrl = mean.control.false.use[idx],
                                 stringsAsFactors = FALSE,
                                 row.names = NULL))
             }) # end scoring bait


    ## Combine results
    saint = do.call(rbind, res.list)

    if (is.null(saint) || !nrow(saint)) {
      warning("No bait-prey pairs could be scored.")
      return(invisible(saint))
    }

    ## Compute Bayesian FDR from the probabilities
    ord = order(saint$AvgP, decreasing = TRUE)
    bfdr = numeric(nrow(saint))
    bfdr[ord] = cumsum(1 - saint$AvgP[ord]) / seq_len(nrow(saint))
    saint$bFDR = bfdr
    saint = saint[ord, , drop = FALSE]
    rownames(saint) = NULL



    # Preparing the output ---------------------------------------------------------------
    ## Re-split back the table in multiple tables and plot volcanoes
    saint.tb.list = list()
    saint.plot.list = list()

    for (i in seq_along(unique(saint$Bait))) {
      saint.tb.list[[i]] = saint[saint$Bait == unique(saint$Bait)[i],]

      saint.plot.list[[i]] =
        ggplot(data = saint.tb.list[[i]],
             aes(x = log2.FoldChange_bait.vs.control,
                 y = -log10(bFDR),
                 size = AvgP,
                 color = AvgP)) +
        geom_point(stroke = NA, alpha = 0.5) +
        ggtitle(paste0("**",unique(saint$Bait)[i],"** *vs* **",control,"**")) +
        xlab(paste0("log<sub>2</sub>(Fold Change<sub>","bait","</sup>/<sub>","ctrl","</sub></sub>)")) +
        ylab("-log<sub>10</sub>(bFDR)") +
        viridis::scale_color_viridis(option = viridis.palette, direction = viridis.direction, begin = viridis.begin, end = viridis.end) +
        ggpubr::theme_pubr(legend = "right") +
        theme(plot.title = ggtext::element_markdown(hjust = 0.5),
              axis.title.x = ggtext::element_markdown(),
              axis.title.y = ggtext::element_markdown(),
              aspect.ratio = 1)
    }

    names(saint.tb.list) = unique(saint$Bait)
    names(saint.plot.list) = unique(saint$Bait)


    parameters = list(prior.pi = pi.t,
                      fold = fold,
                      delta.log2 = delta,
                      sd.scale = sd.scale,
                      min.sd = min.sd,
                      sigma.global = sigma.global,
                      background = bg,
                      which.data = which.data,
                      control = control,
                      baits = bait,
                      best.n.rep = best.n.rep)


    ## Build output object
    saint.out =
      new(Class = "DEprot.SAINTq",
          scores = saint.tb.list,
          volcanoes = saint.plot.list,
          parameters = parameters)


    if (isTRUE(verbose)) {message("Scored ", nrow(saint), " bait-prey pairs across ", length(bait), " bait(s).")}

    return(saint.out)
  }
