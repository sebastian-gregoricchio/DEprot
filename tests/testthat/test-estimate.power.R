## Objects used by the tests: the differential analyses are run once and reused.
## The imputed object is used on purpose, so that the warning on the variance compression
## is raised and can be tested; it is suppressed everywhere else.

dpo.ttest =
  suppressWarnings(suppressMessages(diff.analyses(DEprot.object = DEprot::test.toolbox$dpo.imp,
                                                  contrast.list = list(c("condition", "FBS", "6h.DMSO"),
                                                                       c("condition", "6h.10nM.E2", "6h.DMSO")),
                                                  linear.FC.th = 1.2)))


dpo.paired =
  suppressWarnings(suppressMessages(diff.analyses(DEprot.object = DEprot::test.toolbox$dpo.imp,
                                                  contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                                                  paired.test = TRUE,
                                                  replicate.column = "replicate",
                                                  linear.FC.th = 1.2)))


dpo.wilcox.paired =
  suppressWarnings(suppressMessages(diff.analyses(DEprot.object = DEprot::test.toolbox$dpo.imp,
                                                  contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                                                  stat.test = "wilcoxon",
                                                  paired.test = TRUE,
                                                  replicate.column = "replicate",
                                                  linear.FC.th = 1.2)))


power.ttest = suppressWarnings(estimate.power(DEprot.analyses.object = dpo.ttest, contrast = 1))



######################################    ERRORS    ######################################

test_that("errors are returned if the object is not of class DEprot.analyses", {
  expect_error(estimate.power(DEprot::sample.config))
  expect_error(estimate.power(DEprot::test.toolbox$dpo.imp))
})



test_that("errors are returned if the contrast indicated is not available", {
  expect_error(estimate.power(dpo.ttest, contrast = 100))
  expect_error(estimate.power(dpo.ttest, contrast = "1"))
})



test_that("errors are returned if the parameters are out of range", {
  expect_error(suppressWarnings(estimate.power(dpo.ttest, target.power = 1)))
  expect_error(suppressWarnings(estimate.power(dpo.ttest, target.power = 0)))
  expect_error(suppressWarnings(estimate.power(dpo.ttest, fdr = 0)))
  expect_error(suppressWarnings(estimate.power(dpo.ttest, fdr = 1)))
  expect_error(suppressWarnings(estimate.power(dpo.ttest, pi0 = 1.5)))
  expect_error(suppressWarnings(estimate.power(dpo.ttest, pi0 = 0)))
  expect_error(suppressWarnings(estimate.power(dpo.ttest, approximation = "chisq")))
  expect_error(suppressWarnings(estimate.power(dpo.ttest, effect.size = "median")))
})



test_that("an error is returned if all the proteins are considered responsive", {
  ## pi0 close to zero leaves no null protein, and the FDR cannot be converted into an alpha
  expect_error(suppressWarnings(estimate.power(dpo.ttest, pi0 = 1e-9)))
})



test_that("an error is returned if no protein is left after the effect size filtering", {
  expect_error(suppressWarnings(estimate.power(dpo.ttest, min.effect.size = 1e6)))
})



test_that("an error is returned if the dispersion columns are missing", {
  dpo.stripped = dpo.ttest
  res = dpo.stripped@analyses.result.list[[1]]$results
  dpo.stripped@analyses.result.list[[1]]$results = res[,!grepl("^sd\\.|^sem\\.", colnames(res))]

  expect_error(suppressWarnings(estimate.power(dpo.stripped, contrast = 1)))
})



test_that("an error is returned if none of the effect sizes is finite", {
  dpo.na = dpo.ttest
  sd.columns = grep("^sd\\.", colnames(dpo.na@analyses.result.list[[1]]$results), value = TRUE)
  dpo.na@analyses.result.list[[1]]$results[,sd.columns] = NA_real_

  expect_error(suppressWarnings(estimate.power(dpo.na, contrast = 1)))
})



##########################################################################################

test_that("a warning is raised when the analyses have been run on imputed counts", {
  warns = testthat::capture_warnings(estimate.power(dpo.ttest, contrast = 1))
  expect_true(any(grepl("imputed", warns)))
})



test_that("a warning is raised for paired analyses without lfcSE (wilcoxon)", {
  ## `lfcSE` is NA for rank-based tests, hence the pooled standard deviation is used instead
  warns = testthat::capture_warnings(estimate.power(dpo.wilcox.paired, contrast = 1))
  expect_true(any(grepl("lfcSE", warns)))
})



test_that("no error is returned for the different effect size modes", {
  expect_no_error(suppressWarnings(estimate.power(dpo.ttest, contrast = 1, effect.size = "empirical")))
  expect_no_error(suppressWarnings(estimate.power(dpo.ttest, contrast = 1, effect.size = 1)))
  expect_no_error(suppressWarnings(estimate.power(dpo.ttest, contrast = 1, desired.FC = 1.5)))
  expect_no_error(suppressWarnings(estimate.power(dpo.ttest, contrast = 1, desired.FC = 1.5, sd.quantile = 0.75)))
  expect_no_error(suppressWarnings(estimate.power(dpo.ttest, contrast = 1, effect.size = "empirical", min.effect.size = 0.1)))
})



test_that("no error is returned for the different approximations and options", {
  expect_no_error(suppressWarnings(estimate.power(dpo.ttest, contrast = 1, approximation = "normal")))
  expect_no_error(suppressWarnings(estimate.power(dpo.ttest, contrast = 1, hedges.correction = FALSE)))
  expect_no_error(suppressWarnings(estimate.power(dpo.ttest, contrast = 1, pi0 = 0.9, pi0.lambda = 0.7)))
  expect_no_error(suppressWarnings(estimate.power(dpo.ttest, contrast = 2, fdr = 0.1)))
  expect_no_error(suppressWarnings(estimate.power(dpo.ttest, contrast = 1, max.iterations = 2, tolerance = 0)))
})



test_that("no error is returned for paired analyses", {
  expect_no_error(suppressWarnings(estimate.power(dpo.paired, contrast = 1)))
})



test_that("the columns are resolved also when their suffix does not match the group name", {
  dpo.renamed = dpo.ttest
  res = dpo.renamed@analyses.result.list[[1]]$results
  colnames(res) = gsub("^sd\\..*$", "sd.x", colnames(res))
  colnames(res)[grep("^sd\\.x$", colnames(res))[2]] = "sd.y"
  dpo.renamed@analyses.result.list[[1]]$results = res

  expect_no_error(suppressWarnings(estimate.power(dpo.renamed, contrast = 1)))
})



##########################################################################################

test_that("the object returned is a DEprot.power object with all its slots filled", {
  expect_s4_class(power.ttest, "DEprot.power")
  expect_s3_class(power.ttest@power.plot, "ggplot")
  expect_s3_class(power.ttest@discoveries.plot, "ggplot")
  expect_s3_class(power.ttest@effect.size.plot, "ggplot")
  expect_true(is.data.frame(power.ttest@power.table))
  expect_true(is.numeric(power.ttest@effect.size))
  expect_true(power.ttest@pi0 > 0 & power.ttest@pi0 <= 1)
  expect_equal(power.ttest@params$m, power.ttest@params$m0 + power.ttest@params$m1)
})



test_that("the power table covers the range requested and starts at 2 samples per group", {
  pwr = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, sample.size.range = c(0, 12)))

  expect_equal(min(pwr@power.table$n.per.group), 2)
  expect_equal(max(pwr@power.table$n.per.group), 12)
  expect_equal(nrow(pwr@power.table), 11)
  expect_true(all(c("n.per.group", "alpha", "average.power",
                    "expected.TP", "expected.FP",
                    "expected.discoveries", "expected.FDR") %in% colnames(pwr@power.table)))
})



test_that("the power is bounded and increases with the number of samples", {
  expect_true(all(power.ttest@power.table$average.power >= 0 & power.ttest@power.table$average.power <= 1))
  expect_true(all(diff(power.ttest@power.table$average.power) >= -1e-8))
  expect_true(all(power.ttest@power.table$alpha > 0 & power.ttest@power.table$alpha <= 1))
})



test_that("the realized FDR corresponds to the level requested", {
  ## by construction FP/(TP+FP) must equal the FDR at every sample size
  pwr = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, fdr = 0.05))

  expect_equal(pwr@power.table$expected.FDR,
               rep(0.05, nrow(pwr@power.table)),
               tolerance = 1e-6)
})



test_that("the expected discoveries are consistent with the power and the alpha", {
  tb = power.ttest@power.table

  expect_equal(tb$expected.TP, power.ttest@params$m1 * tb$average.power, tolerance = 1e-8)
  expect_equal(tb$expected.FP, power.ttest@params$m0 * tb$alpha, tolerance = 1e-8)
  expect_equal(tb$expected.discoveries, tb$expected.TP + tb$expected.FP, tolerance = 1e-8)
  expect_true(all(tb$expected.TP <= power.ttest@params$m1 + 1e-8))
})



test_that("the required sample size is the first one reaching the target power", {
  pwr = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, target.power = 0.8))

  if (is.finite(pwr@n.required)) {
    reached = pwr@power.table$n.per.group[pwr@power.table$average.power >= 0.8]
    expect_equal(pwr@n.required, min(reached))
  } else {
    expect_true(all(pwr@power.table$average.power < 0.8))
  }
})



test_that("the required sample size is NA when the target is not reached", {
  pwr = suppressWarnings(estimate.power(dpo.ttest,
                                        contrast = 1,
                                        sample.size.range = c(2, 3),
                                        target.power = 0.999))

  expect_true(is.na(pwr@n.required))
  expect_s3_class(pwr@power.plot, "ggplot")
})



test_that("a null effect size gives no power at any sample size", {
  pwr = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, effect.size = 0))

  expect_true(all(pwr@power.table$average.power < 1e-3))
  expect_true(is.na(pwr@n.required))
})



test_that("a larger effect size requires fewer samples", {
  small = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, effect.size = 0.8, sample.size.range = c(2, 60)))
  large = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, effect.size = 2.5, sample.size.range = c(2, 60)))

  expect_true(large@n.required <= small@n.required)
})



test_that("a stricter FDR requires more samples", {
  loose = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, fdr = 0.2, sample.size.range = c(2, 60)))
  strict = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, fdr = 0.001, sample.size.range = c(2, 60)))

  expect_true(all(strict@power.table$alpha < loose@power.table$alpha))
  expect_true(all(strict@power.table$average.power <= loose@power.table$average.power + 1e-8))
  expect_true(is.na(strict@n.required) | isTRUE(strict@n.required >= loose@n.required))
})



test_that("the effect sizes correspond to the fold changes over the pooled standard deviation", {
  res = dpo.ttest@analyses.result.list[[1]]$results
  vars = dpo.ttest@contrasts[[1]]

  n.1 = round((res[,paste0("sd.", vars$var.1)] / res[,paste0("sem.", vars$var.1)])^2)
  n.2 = round((res[,paste0("sd.", vars$var.2)] / res[,paste0("sem.", vars$var.2)])^2)

  ## proteins with a null variance in one group give 0/0: the function falls back on the group sizes
  n.1[!is.finite(n.1)] = length(vars$group.1)
  n.2[!is.finite(n.2)] = length(vars$group.2)

  sd.pooled = sqrt(((n.1 - 1) * res[,paste0("sd.", vars$var.1)]^2 + (n.2 - 1) * res[,paste0("sd.", vars$var.2)]^2) / (n.1 + n.2 - 2))
  d = res[,paste0("log2.Fold_", vars$var.1, ".vs.", vars$var.2)] / sd.pooled

  finite.d = is.finite(d) & is.finite(res$p.value)
  d = d[finite.d]
  p = res$p.value[finite.d]

  pwr = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, hedges.correction = FALSE))

  expect_equal(sort(pwr@effect.size),
               sort(abs(d[order(p)][1:length(pwr@effect.size)])),
               tolerance = 1e-8)
})



test_that("the Hedges correction shrinks the effect sizes", {
  corrected = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, hedges.correction = TRUE))
  uncorrected = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, hedges.correction = FALSE))

  expect_true(median(corrected@effect.size) < median(uncorrected@effect.size))
})



test_that("the paired analyses use the standard error of the differences", {
  res = dpo.paired@analyses.result.list[[1]]$results
  vars = dpo.paired@contrasts[[1]]

  n.1 = round((res[,paste0("sd.", vars$var.1)] / res[,paste0("sem.", vars$var.1)])^2)
  n.2 = round((res[,paste0("sd.", vars$var.2)] / res[,paste0("sem.", vars$var.2)])^2)

  ## proteins with a null variance in one group give 0/0: the function falls back on the group sizes
  n.1[!is.finite(n.1)] = length(vars$group.1)
  n.2[!is.finite(n.2)] = length(vars$group.2)

  n.pairs = pmin(n.1, n.2)

  d = res[,paste0("log2.Fold_", vars$var.1, ".vs.", vars$var.2)] / (res$lfcSE * sqrt(n.pairs))

  finite.d = is.finite(d) & is.finite(res$p.value)
  d = d[finite.d]
  p = res$p.value[finite.d]

  pwr = suppressWarnings(estimate.power(dpo.paired, contrast = 1, hedges.correction = FALSE))

  expect_equal(sort(pwr@effect.size),
               sort(abs(d[order(p)][1:length(pwr@effect.size)])),
               tolerance = 1e-8)
})



test_that("the parameters used are stored in the object", {
  pwr = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, fdr = 0.1, target.power = 0.9, desired.FC = 2))

  expect_equal(pwr@params$fdr, 0.1)
  expect_equal(pwr@params$target.power, 0.9)
  expect_equal(pwr@params$paired.test, FALSE)
  expect_equal(pwr@params$stat.test, "t.test")
  expect_true(grepl("desired.FC", pwr@params$effect.size.mode))
  expect_equal(length(pwr@effect.size), 1)
})



test_that("the pi0 provided is used instead of the estimated one", {
  pwr = suppressWarnings(estimate.power(dpo.ttest, contrast = 1, pi0 = 0.75))

  expect_equal(pwr@pi0, 0.75)
  expect_equal(pwr@params$m1, round(pwr@params$m * 0.25))
})



##########################################################################################

test_that("the show and plot methods run without errors", {
  expect_output(show(power.ttest))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(power.ttest))
  expect_no_error(plot(power.ttest, nrow = 1))
})



test_that("the function runs on the limma object shipped with the package", {
  res = DEprot::test.toolbox$diff.exp.limma@analyses.result.list[[1]]$results
  skip_if_not(any(grepl("^sd\\.", colnames(res))), "the shipped object predates the dispersion columns")

  expect_no_error(suppressWarnings(estimate.power(DEprot::test.toolbox$diff.exp.limma, contrast = 1)))
})



test_that("the function runs on a proDA object and uses the detected sample sizes", {
  skip_if_not_installed("proDA")

  dpo.proDA =
    suppressWarnings(suppressMessages(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.norm,
                                                          contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                                                          check.missingness = FALSE,
                                                          linear.FC.th = 1.2)))

  pwr = suppressWarnings(estimate.power(dpo.proDA, contrast = 1))

  expect_s4_class(pwr, "DEprot.power")
  expect_equal(pwr@params$stat.test, "proDA")
  expect_equal(pwr@params$counts.used, "normalized")
  expect_true(all(is.finite(pwr@power.table$average.power)))
})
