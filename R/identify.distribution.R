#' @title identify.distribution
#'
#' @description Identifies the best-fitting theoretical distribution (normal, t, chi-squared, or F)
#' for a given vector (such as test statistics) using maximum likelihood estimation.
#' Comparison is based on AIC (Akaike Information Criterion) and BIC (Bayesian Information Criterion) values,
#' and a Kolmogorov-Smirnov test is performed on the best-fitting distribution.
#'
#' @param values Numeric vector of test statistic values (e.g., from a Wilcoxon or t-test). \code{NA} values are removed internally.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{best}: }{String indicating the best-fitting distribution based on AIC.}
#'   \item{\code{comparison}: }{Data frame with AIC and BIC values for all successfully fitted distributions.}
#'   \item{\code{fits}: }{Named list of \code{fitdist} objects for each fitted distribution.}
#'   \item{\code{ks_test_best}: }{Kolmogorov-Smirnov test result for the best-fitting distribution.}
#' }
#'
#' @import fitdistrplus
#'
#' @name identify.distribution
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' result <- identify.distribution(values = DEprot::test.toolbox$diff.exp.limma@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$results$statistic)
#'
#' result$best
#' result$comparison
#'
#' @export identify.distribution



identify.distribution =
  function(values) {

    stats = as.numeric(na.omit(values))

    # Define scaled/shifted t
    dt.scaled = function(x, df, mean = 0, sd = 1) dt((x - mean) / sd, df) / sd
    pt.scaled = function(q, df, mean = 0, sd = 1) pt((q - mean) / sd, df)
    qt.scaled = function(p, df, mean = 0, sd = 1) qt(p, df) * sd + mean

    # Define shifted/scaled F
    df.scaled = function(x, df1, df2, mean = 0, sd = 1) df((x - mean) / sd, df1, df2) / sd
    pf.scaled = function(q, df1, df2, mean = 0, sd = 1) pf((q - mean) / sd, df1, df2)
    qf.scaled = function(p, df1, df2, mean = 0, sd = 1) qf(p, df1, df2) * sd + mean

    # Define shifted/scaled chi-squared
    dchisq.scaled = function(x, df, mean = 0, sd = 1) dchisq((x - mean) / sd, df) / sd
    pchisq.scaled = function(q, df, mean = 0, sd = 1) pchisq((q - mean) / sd, df)
    qchisq.scaled = function(p, df, mean = 0, sd = 1) qchisq(p, df) * sd + mean

    results = list()

    # Normal
    results$norm = tryCatch(fitdist(data = stats,
                                      distr = "norm"),
                              error = function(e) NULL)

    # t
    results$t = tryCatch(fitdist(data = stats,
                                 distr = "t.scaled",
                                 start = list(df = 5, mean = median(stats), sd = sd(stats)),
                                 lower = c(df = 1, mean = -Inf, sd = 0.01)),
                         error = function(e) NULL)

    # Chi-squared
    results$chisq = tryCatch(fitdist(data = stats,
                                     distr = "chisq.scaled",
                                     start = list(df = 5, mean = min(stats) - 0.1, sd = sd(stats)),
                                     lower = c(df = 1, mean = -Inf, sd = 0.01)),
                             error = function(e) NULL)

    # F
    results$f = tryCatch(fitdist(data = stats,
                                 distr = "f.scaled",
                                 start = list(df1 = 5, df2 = 10, mean = min(stats) - 0.1, sd = sd(stats)),
                                 lower = c(df1 = 1, df2 = 1, mean = -Inf, sd = 0.01)),
                           error = function(e) NULL)

    # Collect AIC/BIC
    valid = Filter(Negate(is.null), results)

    if (length(valid) == 0){
      warning("No distribution could be fitted.")
      return(list(best = NA,
                  comparison = NA,
                  fits = NA,
                  ks_test_best = NA))
    }

    comparison = data.frame(distribution = names(valid),
                            AIC = sapply(valid, function(f) f$aic),
                            BIC = sapply(valid, function(f) f$bic),
                            row.names = NULL)

    comparison = comparison[order(comparison$AIC), ]
    best = comparison$distribution[1]

    # KS test for best fit
    best_fit = valid[[best]]
    ks = do.call(ks.test,
                 c(list(x = stats, y = paste0("p", gsub("\\.", ".", best_fit$distname))),
                   as.list(best_fit$estimate)))

    return(list(best = best,
                comparison = comparison,
                fits = valid,
                ks_test_best = ks))

  } # END function



