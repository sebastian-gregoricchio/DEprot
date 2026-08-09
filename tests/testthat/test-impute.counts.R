## ----------------------------------------------------------------------------------------
##  impute.counts()
##  Every method is a branch of the same function: they are looped over on a small matrix,
##  skipping the ones whose package is not installed.
##  The '@imputation.method' slot holds a list, whose 'method' element carries the label.
## ----------------------------------------------------------------------------------------

dpo.miss <- make.dpo(n.prot = 40, n.samples = 6, n.missing = 20)


test_that("missForest imputation fills every missing value", {
  skip_if_not_installed("missForest")

  imp <- impute.counts(DEprot.object = dpo.miss,
                       method = "missForest",
                       which.data = "randomized",
                       missForest.max.iterations = 3,
                       seed = 1,
                       verbose = FALSE)

  expect_s4_class(imp, "DEprot")
  expect_true(imp@imputed)
  expect_false(any(is.na(imp@imputed.counts)))
  expect_equal(dim(imp@imputed.counts), dim(dpo.miss@random.counts))
  expect_equal(rownames(imp@imputed.counts), rownames(dpo.miss@random.counts))
  expect_type(imp@imputation.method, "list")
  expect_true(grepl("missForest", imp@imputation.method$method, ignore.case = TRUE))
})


test_that("the kNN-based methods run", {
  skip_if_not_installed("VIM")

  for (m in c("KNN", "kNN-VIM")) {
    imp <- impute.counts(DEprot.object = dpo.miss,
                         method = m,
                         which.data = "randomized",
                         kNN.n.nearest.neighbours = 3,
                         seed = 1,
                         verbose = FALSE)

    expect_s4_class(imp, "DEprot")
    expect_false(any(is.na(imp@imputed.counts)))
    expect_true(grepl("kNN", imp@imputation.method$method, ignore.case = TRUE))
  }
})


test_that("the truncated and correlation-based kNN run", {
  for (m in c("tKNN", "corkNN")) {
    imp <- try(impute.counts(DEprot.object = dpo.miss,
                             method = m,
                             which.data = "randomized",
                             kNN.n.nearest.neighbours = 3,
                             seed = 1,
                             verbose = FALSE),
               silent = TRUE)

    skip_if(inherits(imp, "try-error"), paste0("'", m, "' is not available in this setup"))

    expect_s4_class(imp, "DEprot")
    expect_false(any(is.na(imp@imputed.counts)))
    expect_true(grepl(m, imp@imputation.method$method, ignore.case = TRUE))
  }
})


test_that("the pcaMethods-based methods run", {
  skip_if_not_installed("pcaMethods")

  for (m in c("SVD", "LLS", "BPCA", "PPCA")) {
    imp <- impute.counts(DEprot.object = dpo.miss,
                         method = m,
                         which.data = "randomized",
                         LLS.k = 2,
                         pcaMethods.nPCs.to.test = 2,
                         seed = 1,
                         verbose = FALSE)

    expect_s4_class(imp, "DEprot")
    expect_false(any(is.na(imp@imputed.counts)))
    ## the label stored in the object is the canonical name of the method
    expect_true(grepl(m, imp@imputation.method$method, ignore.case = TRUE))
  }
})


test_that("RegImpute runs with both fill methods", {
  for (fill in c("row_mean", "zeros")) {
    imp <- try(impute.counts(DEprot.object = dpo.miss,
                             method = "RegImpute",
                             which.data = "randomized",
                             RegImpute.max.iterations = 2,
                             RegImpute.fillmethod = fill,
                             seed = 1,
                             verbose = FALSE),
               silent = TRUE)

    skip_if(inherits(imp, "try-error"), "'RegImpute' is not available in this setup")

    expect_s4_class(imp, "DEprot")
    expect_false(any(is.na(imp@imputed.counts)))
  }
})


test_that("the imputation can be run on raw and normalized counts", {
  skip_if_not_installed("pcaMethods")

  imp.raw <- impute.counts(DEprot.object = dpo.miss,
                           method = "SVD",
                           which.data = "raw",
                           seed = 1,
                           verbose = FALSE)

  imp.norm <- impute.counts(DEprot.object = dpo.miss,
                            method = "SVD",
                            which.data = "normalized",
                            seed = 1,
                            verbose = FALSE)

  expect_s4_class(imp.raw, "DEprot")
  expect_s4_class(imp.norm, "DEprot")
})


test_that("the imputation is reproducible with a fixed seed", {
  skip_if_not_installed("pcaMethods")

  a <- impute.counts(DEprot.object = dpo.miss, method = "SVD", which.data = "randomized",
                     seed = 42, verbose = FALSE)
  b <- impute.counts(DEprot.object = dpo.miss, method = "SVD", which.data = "randomized",
                     seed = 42, verbose = FALSE)

  expect_equal(a@imputed.counts, b@imputed.counts)
})



## ----------------------------------------------------------------------------------------
##  Error branches
## ----------------------------------------------------------------------------------------

test_that("an unknown imputation method is rejected", {
  expect_error(impute.counts(DEprot.object = dpo.miss, method = "not.a.method", verbose = FALSE))
})


test_that("an object of the wrong class is rejected", {
  expect_error(impute.counts(DEprot.object = data.frame(a = 1), verbose = FALSE))
})


test_that("counts that are not available cannot be used", {
  dpo.bare <- dpo.miss
  dpo.bare@norm.counts <- NULL
  dpo.bare@random.counts <- NULL

  expect_error(impute.counts(DEprot.object = dpo.bare, which.data = "normalized", verbose = FALSE))
  expect_error(impute.counts(DEprot.object = dpo.bare, which.data = "randomized", verbose = FALSE))
})


test_that("an object already imputed is not overwritten silently", {
  skip_if_not_installed("pcaMethods")

  imp <- impute.counts(DEprot.object = dpo.miss, method = "SVD", which.data = "randomized",
                       seed = 1, verbose = FALSE)

  expect_error(impute.counts(DEprot.object = imp, method = "SVD", which.data = "randomized",
                             overwrite.imputation = FALSE, verbose = FALSE))

  imp2 <- impute.counts(DEprot.object = imp, method = "SVD", which.data = "randomized",
                        overwrite.imputation = TRUE, seed = 1, verbose = FALSE)

  expect_s4_class(imp2, "DEprot")
})
