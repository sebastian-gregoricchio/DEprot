######################################    ERRORS    ######################################
test_that("errors are returned if the table is not numeric", {
  cnt = DEprot::unimputed.counts
  cnt[8,1:10] = "A"
  expect_error(load.counts2(counts = iris,
                            metadata = DEprot::sample.config[-5,],
                            data.type = "norm",
                            log.base = 2))

})



test_that("errors are returned if the counts and metadata do not match", {
  expect_error(load.counts2(counts = DEprot::unimputed.counts,
                            metadata = DEprot::sample.config[1:6,],
                            data.type = "norm",
                            log.base = 2))

})



test_that("errors are returned if the column.id idnidcated is not present in the metadata", {
  expect_error(load.counts2(counts = DEprot::unimputed.counts,
                            metadata = DEprot::sample.config[1:6,],
                            data.type = "norm",
                            log.base = 2,
                            column.id = "not a column"))

})


##########################################################################################

test_that("the function load.counts.2 is working with normalized data", {
  expect_no_error(load.counts2(counts = DEprot::unimputed.counts,
                               metadata = DEprot::sample.config,
                               data.type = "norm",
                               log.base = 2))

})


test_that("the function load.counts.2 is working with raw data", {
  expect_no_error(load.counts2(counts = DEprot::unimputed.counts,
                               metadata = DEprot::sample.config,
                               data.type = "raw",
                               log.base = 2))

})


test_that("the function load.counts.2 is working with imputed data", {
  expect_no_error(load.counts2(counts = DEprot::unimputed.counts,
                               metadata = DEprot::sample.config,
                               data.type = "imputed",
                               log.base = 2))

})


test_that("the function load.counts.2 is working with randomized data", {
  expect_no_error(load.counts2(counts = DEprot::unimputed.counts,
                               metadata = DEprot::sample.config,
                               data.type = "randomized",
                               log.base = 2))

})



test_that("the function load.counts.2 is working with linear data", {
  expect_no_error(suppressMessages(load.counts2(counts = DEprot::unimputed.counts,
                                                metadata = DEprot::sample.config,
                                                data.type = "randomized",
                                                log.base = 1)))

})


##########################################################################################
##                                    protein.info                                      ##
##########################################################################################

test_that("the protein.info slot is NULL when no annotation is provided", {
  dpo = load.counts2(counts = DEprot::unimputed.counts,
                     metadata = DEprot::sample.config,
                     data.type = "raw",
                     log.base = 2)
  expect_null(dpo@protein.info)
})


test_that("an annotation provided at loading is stored and aligned to the counts", {
  prot.ids = rownames(DEprot::unimputed.counts)
  info = data.frame(gene.name = paste0("GENE", seq_along(prot.ids)),
                    row.names = prot.ids,
                    stringsAsFactors = FALSE)

  dpo = suppressMessages(load.counts2(counts = DEprot::unimputed.counts,
                                      metadata = DEprot::sample.config,
                                      data.type = "raw",
                                      log.base = 2,
                                      protein.info = info))

  expect_true(is.data.frame(dpo@protein.info))
  expect_equal(rownames(dpo@protein.info), rownames(dpo@raw.counts))

  # load.counts2 discards the rows containing only NAs: the annotation follows the
  # proteins effectively retained in the object, not the ones of the input table
  expect_lte(nrow(dpo@protein.info), length(prot.ids))
  expect_equal(dpo@protein.info$gene.name, info[rownames(dpo@raw.counts), "gene.name"])
})



test_that("the identifiers of the annotation can be given in a dedicated column", {
  prot.ids = rownames(DEprot::unimputed.counts)
  info = data.frame(uniprot = prot.ids,
                    gene.name = paste0("GENE", seq_along(prot.ids)),
                    stringsAsFactors = FALSE)

  dpo = load.counts2(counts = DEprot::unimputed.counts,
                     metadata = DEprot::sample.config,
                     data.type = "raw",
                     log.base = 2,
                     protein.info = info,
                     protein.info.id.column = "uniprot")

  expect_equal(rownames(dpo@protein.info), rownames(dpo@raw.counts))
  expect_false("uniprot" %in% colnames(dpo@protein.info))
})


test_that("errors are returned if the annotation contains duplicated protein IDs", {
  prot.ids = rownames(DEprot::unimputed.counts)
  info = data.frame(gene.name = paste0("GENE", seq_along(prot.ids)),
                    row.names = prot.ids,
                    stringsAsFactors = FALSE)
  info$prot.id = c(prot.ids[1], prot.ids[-length(prot.ids)])   # first ID duplicated
  rownames(info) = NULL

  expect_error(load.counts2(counts = DEprot::unimputed.counts,
                            metadata = DEprot::sample.config,
                            data.type = "raw",
                            log.base = 2,
                            protein.info = info))
})
