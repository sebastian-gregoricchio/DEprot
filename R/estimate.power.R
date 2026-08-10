#' @title estimate.power
#'
#' @description Estimates the average power and the number of replicates required by a future experiment,
#' starting from the effect sizes and the variances observed in a differential analysis. The multiplicity
#' is handled by controlling the False Discovery Rate following the approach of Jung (2005) and its
#' closed-form simplification by Liu & Hwang (2007). The calculation is self-contained: it only uses the
#' fold changes and the per-group standard deviations stored in the results table, and requires no
#' external power/sample-size package.
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast Number indicating the position of the contrast to use. Default: \code{1}.
#' @param sample.size.range Numeric vector of length 2 indicating the minimum and maximum number of samples \emph{per group} to test. Default: \code{c(2,30)}.
#' @param target.power Number between 0 and 1 indicating the average power to reach. Default: \code{0.8}.
#' @param fdr Number between 0 and 1 indicating the False Discovery Rate to control. Default: \code{NULL}, the \code{padj.th} stored in the object is used.
#' @param pi0 Number between 0 and 1 indicating the proportion of non-differential proteins. Default: \code{NULL}, estimated from the p-value distribution.
#' @param pi0.lambda Number between 0 and 1 indicating the tuning parameter used to estimate \code{pi0} (Storey's estimator: \code{mean(p > lambda)/(1 - lambda)}). Default: \code{0.5}.
#' @param effect.size Either the string 'empirical', to average the power over the effect sizes observed for the responsive proteins, or a numeric value indicating a single standardized effect size (Cohen's d) to plan for. Default: \code{"empirical"}.
#' @param desired.FC Numeric value indicating a desired fold change (linear scale) to plan for. When provided, the effect size is computed as \code{log2(desired.FC)/sd}, where \code{sd} is the quantile of the observed pooled standard deviations defined by \code{sd.quantile}. It overrides \code{effect.size}. Default: \code{NULL}.
#' @param sd.quantile Number between 0 and 1 indicating the quantile of the observed pooled standard deviations to use together with \code{desired.FC}. Default: \code{0.5} (median).
#' @param min.effect.size Numeric value indicating the minimum absolute effect size to keep among the responsive proteins. Default: \code{NULL}, no filtering.
#' @param hedges.correction Logical value indicating whether the small-sample bias of the standardized effect sizes should be corrected (Hedges' g). Default: \code{TRUE}.
#' @param approximation String indicating the distribution used for the power calculation. One among: 't' (non-central t, exact) and 'normal' (normal approximation, as in the original Jung's formulation). Default: \code{"t"}.
#' @param max.iterations Numeric value indicating the maximum number of iterations allowed to solve the average power at each sample size. Default: \code{100}.
#' @param tolerance Numeric value indicating the tolerance used to solve the average power at each sample size. Default: \code{1e-6}.
#' @param line.color String indicating any R-supported color to use for the curves. Default: \code{"black"}.
#' @param power.threshold.line.color String indicating any R-supported color to use for the target power line. Default: \code{"firebrick"}.
#'
#' @return An object of class \code{DEprot.power}.
#'
#' @details For a given number of replicates \code{n} per group, the average power over the \code{m1} responsive
#' proteins and the per-test significance level \code{alpha} are linked by the FDR:
#' \code{alpha = FDR * m1 * (1-beta) / ((1 - FDR) * m0)}, with \code{m1 = m * (1 - pi0)} and \code{m0 = m - m1}.
#' Since the two quantities depend on each other, at each sample size the couple is obtained by root finding
#' on \code{power(alpha(x)) - x}, taking the largest root because \code{x = 0} is always a trivial solution.
#' The realized FDR reported in the results is therefore equal to the level requested by construction.
#' The power of the single protein is obtained from a two-sided two-sample t-test with non-centrality parameter
#' \code{d * sqrt(n/2)} (\code{d * sqrt(n)} and \code{df = n-1} for paired designs).
#'
#' The standardized effect sizes are computed as \code{log2(FoldChange)} divided by the pooled standard deviation
#' reconstructed from the \code{sd.<group>} columns, and the per-protein group sizes are recovered from the
#' \code{sd}/\code{sem} pair, so that proteins with missing values are treated correctly. For paired analyses
#' the effect size is a Cohen's dz, obtained from \code{lfcSE}; when \code{lfcSE} is not available
#' (\code{stat.test = "wilcoxon"}) the pooled standard deviation is used instead, which ignores the
#' within-pair correlation and makes the estimate conservative.
#'
#' The function works on the output of all the differential functions of DEprot (\code{diff.analyses},
#' \code{diff.analyses.limma}, \code{diff.analyses.prolfqua} and \code{diff.analyses.proDA}), since it only needs
#' the fold changes and the per-group dispersions. For \code{diff.analyses.proDA} the standard deviations are
#' computed on the observed values only, so for the proteins affected by dropout they are estimated on the upper
#' part of the intensity distribution and are under-estimated: restricting the estimation to the proteins
#' quantified in all the samples (\code{min.effect.size} is not enough for this) gives a more conservative answer.
#'
#' Two further limitations are worth keeping in mind. Effect sizes estimated from a pilot experiment are inflated for the
#' proteins that were selected \emph{because} they came out significant (winner's curse), so the empirical mode is
#' optimistic when \code{pi0} is high; planning on \code{desired.FC} avoids the problem. Moreover, when the
#' differential analyses have been run on imputed counts the variance is compressed and every estimate below
#' becomes optimistic: a warning is raised in that case.
#'
#' @references Jung SH (2005), Sample size for FDR-control in microarray data analysis. \emph{Bioinformatics} 21(14):3097-3104.
#' @references Liu P, Hwang JTG (2007), Quick calculation for sample size while controlling false discovery rate with application to microarray analysis. \emph{Bioinformatics} 23(6):739-746.
#'
#' @author Sebastian Gregoricchio
#'
#' @import ggplot2
#' @import dplyr
#' @import ggtext
#' @importFrom stats qt pt qnorm pnorm quantile density
#'
#' @examples
#' pwr <- estimate.power(DEprot::test.toolbox$diff.exp.limma, contrast = 1)
#' pwr
#' pwr@power.table
#'
#' @export estimate.power


estimate.power =
  function(DEprot.analyses.object,
           contrast = 1,
           sample.size.range = c(2, 30),
           target.power = 0.8,
           fdr = NULL,
           pi0 = NULL,
           pi0.lambda = 0.5,
           effect.size = "empirical",
           desired.FC = NULL,
           sd.quantile = 0.5,
           min.effect.size = NULL,
           hedges.correction = TRUE,
           approximation = "t",
           max.iterations = 100,
           tolerance = 1e-6,
           line.color = "black",
           power.threshold.line.color = "firebrick") {

    # ## Libraries
    # require(dplyr)
    # require(ggplot2)

    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      stop("The input must be an object of class 'DEprot.analyses'.")
    }

    ### check and collect contrast
    if (!is.numeric(contrast)) {
      stop("The 'contrast' must be a numeric value.")
    } else if (contrast > length(DEprot.analyses.object@analyses.result.list)) {
      stop("The 'contrast' indicated is not available.")
    }

    results = DEprot.analyses.object@analyses.result.list[[contrast]]$results
    contrasts.info = DEprot.analyses.object@contrasts[[contrast]]
    params = DEprot.analyses.object@differential.analyses.params

    ### check other parameters
    if (!(tolower(approximation) %in% c("t", "normal", "norm"))) {
      stop("The `approximation` must be one among: 't', 'normal'.")
    }

    if (target.power <= 0 | target.power >= 1) {
      stop("The `target.power` must be a number between 0 and 1.")
    }

    if (is.null(fdr)) {fdr = params$padj.th}
    if (fdr <= 0 | fdr >= 1) {stop("The `fdr` must be a number between 0 and 1.")}

    sample.size.range = sort(round(sample.size.range))
    if (sample.size.range[1] < 2) {sample.size.range[1] = 2}


    ### warn when the analyses have been run on imputed counts
    if (grepl("imput", tolower(paste0(params$counts.used, collapse = " ")))) {
      warning(paste0("The differential analyses have been performed on imputed counts ('", params$counts.used, "'): ",
                     "the imputation compresses the variance and all the estimates below are optimistic.\n",
                     "Consider re-running the differential analyses on the normalized counts, restricted to the proteins quantified in all the samples of the contrast."), call. = FALSE)
    }




    ################################################################
    ### Collect the columns required (their name depends on the contrast)
    get.column = function(tb, prefix, group) {
      col = paste0(prefix, ".", group)
      if (!(col %in% colnames(tb))) {col = grep(paste0("^", prefix, "\\."), colnames(tb), value = TRUE)[ifelse(group == contrasts.info$var.1, 1, 2)]}
      if (length(col) == 0) {col = NA_character_}
      return(col)
    }

    sd.1.col = get.column(results, "sd", contrasts.info$var.1)
    sd.2.col = get.column(results, "sd", contrasts.info$var.2)
    sem.1.col = get.column(results, "sem", contrasts.info$var.1)
    sem.2.col = get.column(results, "sem", contrasts.info$var.2)
    fc.col = grep("^log2\\.Fold_", colnames(results), value = TRUE)[1]

    if (any(is.na(c(sd.1.col, sd.2.col, sem.1.col, sem.2.col, fc.col)))) {
      stop("The results table does not contain the `sd.<group>`/`sem.<group>` columns required: the object has probably been generated with an older version of DEprot.")
    }


    ### Per-protein group sizes: taken from the `n.detected.<group>` columns when available
    ### (diff.analyses.proDA), otherwise recovered from the sd/sem pair (sem = sd/sqrt(n))
    n.det.1.col = get.column(results, "n.detected", contrasts.info$var.1)
    n.det.2.col = get.column(results, "n.detected", contrasts.info$var.2)

    if (!is.na(n.det.1.col) & !is.na(n.det.2.col)) {
      n.1 = results[,n.det.1.col]
      n.2 = results[,n.det.2.col]
    } else {
      n.1 = round((results[,sd.1.col] / results[,sem.1.col])^2)
      n.2 = round((results[,sd.2.col] / results[,sem.2.col])^2)
    }

    n.1[!is.finite(n.1)] = length(contrasts.info$group.1)
    n.2[!is.finite(n.2)] = length(contrasts.info$group.2)

    paired.test = isTRUE(contrasts.info$paired.test)




    ################################################################
    ### Standardized effect sizes
    sd.pooled = sqrt(((n.1 - 1) * results[,sd.1.col]^2 + (n.2 - 1) * results[,sd.2.col]^2) / (n.1 + n.2 - 2))

    if (paired.test == TRUE) {
      if ("lfcSE" %in% colnames(results) & any(is.finite(results$lfcSE))) {
        n.pairs = pmin(n.1, n.2)
        sd.ref = results$lfcSE * sqrt(n.pairs)
        df.pilot = n.pairs - 1
      } else {
        warning("The analyses are paired but `lfcSE` is not available (Wilcoxon test): the pooled standard deviation is used instead, which ignores the within-pair correlation and makes the estimate conservative.", call. = FALSE)
        sd.ref = sd.pooled
        df.pilot = n.1 + n.2 - 2
      }
    } else {
      sd.ref = sd.pooled
      df.pilot = n.1 + n.2 - 2
    }

    d.observed = results[,fc.col] / sd.ref

    ## Hedges' correction of the small-sample bias
    if (hedges.correction == TRUE) {
      J = 1 - 3 / (4 * df.pilot - 1)
      J[!is.finite(J) | J <= 0] = 1
      d.observed = d.observed * J
    }

    finite.d = is.finite(d.observed) & is.finite(results$p.value)

    if (sum(finite.d) == 0) {
      stop("None of the proteins has a finite effect size: the power cannot be estimated.")
    }

    m = sum(finite.d)




    ################################################################
    ### Proportion of non-differential proteins
    if (is.null(pi0)) {
      pi0 = mean(results$p.value[finite.d] > pi0.lambda, na.rm = TRUE) / (1 - pi0.lambda)
      pi0 = min(1, max(pi0, 1/m))
    } else if (pi0 <= 0 | pi0 > 1) {
      stop("The `pi0` must be a number between 0 and 1.")
    }

    m1 = max(1, round(m * (1 - pi0)))
    m0 = m - m1

    if (m0 <= 0) {
      stop("The estimated proportion of null proteins is zero: provide `pi0` manually.")
    }




    ################################################################
    ### Effect sizes to plan for
    if (!is.null(desired.FC)) {
      ## MSstats-like mode: a single target fold change over a quantile of the observed dispersion
      sd.target = stats::quantile(sd.ref[finite.d], probs = sd.quantile, na.rm = TRUE)
      d.target = abs(log2(desired.FC)) / as.numeric(sd.target)
      effect.size.mode = paste0("desired.FC = ", desired.FC, " (sd quantile: ", sd.quantile, ")")

    } else if (is.numeric(effect.size)) {
      d.target = abs(effect.size)
      effect.size.mode = paste0("fixed d = ", round(d.target, 3))

    } else if (identical(effect.size, "empirical")) {
      ## the m1 proteins with the smallest p-values are taken as the responsive ones
      responsive.idx = order(results$p.value[finite.d])[1:m1]
      d.target = abs(d.observed[finite.d][responsive.idx])
      d.target = d.target[is.finite(d.target)]

      if (!is.null(min.effect.size)) {
        d.target = d.target[d.target >= min.effect.size]
        if (length(d.target) == 0) {stop("No responsive protein is left after the `min.effect.size` filtering.")}
      }

      effect.size.mode = "empirical"

    } else {
      stop("The `effect.size` must be either 'empirical' or a numeric value.")
    }




    ################################################################
    ### Power function (two-sided test)
    compute.average.power = function(n, alpha) {
      if (paired.test == TRUE) {
        ncp = d.target * sqrt(n)
        df = n - 1
      } else {
        ncp = d.target * sqrt(n / 2)
        df = 2 * n - 2
      }

      if (df < 1) {return(0)}

      if (tolower(approximation) == "t") {
        crit = stats::qt(1 - alpha/2, df = df)
        pwr = stats::pt(crit, df = df, ncp = ncp, lower.tail = FALSE) + stats::pt(-crit, df = df, ncp = ncp, lower.tail = TRUE)
      } else {
        crit = stats::qnorm(1 - alpha/2)
        pwr = stats::pnorm(ncp - crit) + stats::pnorm(-ncp - crit)
      }

      return(mean(pwr, na.rm = TRUE))
    }




    ################################################################
    ### Per-test significance level required to control the FDR at a given average power
    compute.alpha = function(power) {
      return(min((fdr * m1 * power) / ((1 - fdr) * m0), 1))
    }


    ### Solve alpha and average power at a given sample size.
    ### The two define each other, and a plain fixed-point iteration can oscillate at very small
    ### sample sizes: it would then exhaust the iterations and return an alpha that does not
    ### correspond to the power reported. The couple is therefore found by root finding on
    ### `power(alpha(x)) - x`, taking the largest root since x = 0 is always a trivial solution.
    solve.power = function(n) {
      gap = function(x) {compute.average.power(n = n, alpha = compute.alpha(x)) - x}

      grid = c(1e-8, seq(0.005, 1, length.out = 40))
      gap.grid = sapply(grid, gap)

      crossing = which(gap.grid[-length(gap.grid)] > 0 & gap.grid[-1] <= 0)

      if (length(crossing) == 0) {
        return(ifelse(gap.grid[1] <= 0, grid[1], 1))
      }

      k = max(crossing)
      root = try(stats::uniroot(f = gap,
                                interval = c(grid[k], grid[k+1]),
                                tol = tolerance,
                                maxiter = max.iterations)$root,
                 silent = TRUE)

      if (inherits(root, "try-error")) {root = grid[k]}

      return(root)
    }




    ################################################################
    ### Average power at each sample size
    tested.n = seq(sample.size.range[1], sample.size.range[2], 1)

    power.table =
      data.frame(n.per.group = tested.n,
                 alpha = NA_real_,
                 average.power = NA_real_,
                 expected.TP = NA_real_,
                 expected.FP = NA_real_,
                 expected.discoveries = NA_real_,
                 expected.FDR = NA_real_)

    for (k in 1:length(tested.n)) {
      current.power = solve.power(n = tested.n[k])
      alpha = compute.alpha(current.power)

      TP = m1 * current.power
      FP = m0 * alpha

      power.table$alpha[k] = alpha
      power.table$average.power[k] = current.power
      power.table$expected.TP[k] = TP
      power.table$expected.FP[k] = FP
      power.table$expected.discoveries[k] = TP + FP
      power.table$expected.FDR[k] = ifelse((TP + FP) > 0, FP / (TP + FP), NA_real_)
    }


    ### Smallest sample size reaching the target power
    reached = which(power.table$average.power >= target.power)
    n.required = ifelse(length(reached) > 0, power.table$n.per.group[min(reached)], NA_real_)

    n.current = (length(contrasts.info$group.1) + length(contrasts.info$group.2)) / 2




    ################################################################
    ### Plots
    subtitle.label = paste0("*", contrasts.info$metadata.column, ":* **", contrasts.info$var.1, "** *vs* **", contrasts.info$var.2, "**")

    power.plot =
      ggplot(data = power.table,
             aes(x = n.per.group,
                 y = average.power)) +
      geom_hline(yintercept = target.power, linewidth = 0.5, color = power.threshold.line.color) +
      geom_vline(xintercept = n.current, linetype = 3, color = "gray50") +
      geom_line(linewidth = 1, color = line.color) +
      {if (is.finite(n.required)) geom_point(data = power.table[power.table$n.per.group == n.required,], size = 2.5, color = power.threshold.line.color)} +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      scale_x_continuous(expand = c(0.01,0)) +
      xlab("Sample size (per group)") +
      ylab("Average power") +
      ggtitle(label = paste0("Average power at FDR &le; ", fdr, ifelse(is.finite(n.required), paste0(" (*n* = ", n.required, ")"), "")),
              subtitle = subtitle.label) +
      theme_classic() +
      theme(axis.text = element_text(color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.title.x = ggtext::element_markdown(),
            axis.title.y = ggtext::element_markdown(),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5),
            aspect.ratio = 0.6)


    discoveries.plot =
      ggplot(data = power.table,
             aes(x = n.per.group,
                 y = expected.TP)) +
      geom_vline(xintercept = n.current, linetype = 3, color = "gray50") +
      geom_line(linewidth = 1, color = line.color) +
      scale_y_continuous(expand = c(0,0), limits = c(0, m1*1.05)) +
      scale_x_continuous(expand = c(0.01,0)) +
      xlab("Sample size (per group)") +
      ylab("Expected true discoveries") +
      ggtitle(label = paste0("*m* = ", m, ", *m*<sub>1</sub> = ", m1, " (&pi;<sub>0</sub> = ", round(pi0, 3), ")"),
              subtitle = subtitle.label) +
      theme_classic() +
      theme(axis.text = element_text(color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.title.x = ggtext::element_markdown(),
            axis.title.y = ggtext::element_markdown(),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5),
            aspect.ratio = 0.6)


    effect.size.plot =
      ggplot(data = data.frame(effect.size = d.target),
             aes(x = effect.size)) +
      geom_density(fill = line.color, color = line.color, alpha = 0.25) +
      geom_vline(xintercept = median(d.target, na.rm = TRUE), linetype = 2, color = "gray30") +
      scale_x_continuous(expand = c(0.01,0)) +
      scale_y_continuous(expand = c(0,0)) +
      xlab("Standardized effect size (|*d*|)") +
      ylab("Density") +
      ggtitle(label = paste0("Effect sizes | mode: ", effect.size.mode),
              subtitle = paste0("median |*d*| = ", round(median(d.target, na.rm = TRUE), 3))) +
      theme_classic() +
      theme(axis.text = element_text(color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.title.x = ggtext::element_markdown(),
            axis.title.y = ggtext::element_markdown(),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5),
            aspect.ratio = 0.6)




    ################################################################
    ### Export DEprot.power object
    DEprot.power.object =
      new(Class = "DEprot.power",
          power.table = power.table,
          n.required = n.required,
          effect.size = d.target,
          pi0 = pi0,
          power.plot = power.plot,
          discoveries.plot = discoveries.plot,
          effect.size.plot = effect.size.plot,
          params = list(contrast = paste0(contrasts.info$var.1, ".vs.", contrasts.info$var.2),
                        n.current = n.current,
                        paired.test = paired.test,
                        stat.test = params$stat.test,
                        counts.used = params$counts.used,
                        fdr = fdr,
                        target.power = target.power,
                        pi0.lambda = pi0.lambda,
                        effect.size.mode = effect.size.mode,
                        hedges.correction = hedges.correction,
                        approximation = tolower(approximation),
                        m = m,
                        m1 = m1,
                        m0 = m0))

    return(DEprot.power.object)

  } # END function
