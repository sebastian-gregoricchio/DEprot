test_that("get.results works for DEprot.analyses objects", {
  expect_no_error(get.results(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, 1))
})


test_that("get.results does not work for object not of class DEprot.analyses", {
  expect_error(get.results(DEprot.analyses.object = DEprot::test.toolbox$dpo.imp, 1))
})

test_that("get.results does not work if the contrast requested is not available", {
  expect_error(get.results(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, 100))
})


##########################################################################################
##                                    protein.info                                      ##
##########################################################################################

## Objects and a minimal annotation built on the protein IDs of the counts
dpa.noinfo <- DEprot::test.toolbox$diff.exp.limma
prot.ids <- rownames(dpa.noinfo@imputed.counts)

info.3col <- data.frame(gene.name = paste0("GENE", seq_along(prot.ids)),
                        n.peptides = seq_along(prot.ids),
                        description = paste("Protein", seq_along(prot.ids)),
                        row.names = prot.ids,
                        stringsAsFactors = FALSE)

dpa.info <- add.protein.info(dpa.noinfo, info.3col)


test_that("by default no annotation column is appended", {
  res <- get.results(dpa.info, 1)
  expect_false(any(colnames(info.3col) %in% colnames(res)))
})


test_that("the keyword 'all' appends the whole annotation", {
  res <- get.results(dpa.info, 1, protein.info.columns = "all")
  expect_true(all(colnames(info.3col) %in% colnames(res)))
  expect_equal(nrow(res), nrow(get.results(dpa.info, 1)))
})


test_that("the keywords are case-insensitive and NULL is read as 'none'", {
  expect_true("gene.name" %in% colnames(get.results(dpa.info, 1, "All")))
  expect_false("gene.name" %in% colnames(get.results(dpa.info, 1, "NONE")))
  expect_false("gene.name" %in% colnames(get.results(dpa.info, 1, NULL)))
})


test_that("a subset of the annotation columns can be requested", {
  res <- get.results(dpa.info, 1, protein.info.columns = c("gene.name", "n.peptides"))
  expect_true(all(c("gene.name", "n.peptides") %in% colnames(res)))
  expect_false("description" %in% colnames(res))
})


test_that("the annotation columns can be prefixed", {
  res <- get.results(dpa.info, 1, protein.info.columns = "all", protein.info.prefix = "info.")
  expect_true(all(paste0("info.", colnames(info.3col)) %in% colnames(res)))
  expect_false("gene.name" %in% colnames(res))
})


test_that("the annotation is bound by protein ID and not by position", {
  res <- get.results(dpa.info, 1, protein.info.columns = "all")
  expect_equal(res$gene.name, info.3col[as.character(res$prot.id), "gene.name"])
})


test_that("unannotated proteins are returned as NAs", {
  partial <- info.3col[1:10, , drop = FALSE]
  obj <- suppressMessages(add.protein.info(dpa.noinfo, partial))
  res <- get.results(obj, 1, protein.info.columns = "all")
  expect_equal(nrow(res), length(prot.ids))
  expect_true(any(is.na(res$gene.name)))
})


test_that("colliding column names are made unique without touching the results", {
  info <- info.3col
  colnames(info)[1] <- "padj"          # already present in the results table
  obj <- add.protein.info(dpa.noinfo, info)
  expect_warning(res <- get.results(obj, 1, protein.info.columns = "all"), "same name")
  expect_true("padj" %in% colnames(res))
  expect_true("padj.1" %in% colnames(res))
  # the original results column is preserved
  expect_equal(res$padj, get.results(dpa.info, 1)$padj)
})


test_that("requesting non-existing annotation columns returns an error", {
  expect_error(get.results(dpa.info, 1, protein.info.columns = "not.a.column"),
               "not available")
})


test_that("a non-character protein.info.columns returns an error", {
  expect_error(get.results(dpa.info, 1, protein.info.columns = 3), "string vector")
})


test_that("'all' is silent when the object carries no annotation", {
  expect_no_warning(res <- get.results(dpa.noinfo, 1, protein.info.columns = "all"))
  expect_equal(colnames(res), colnames(get.results(dpa.noinfo, 1)))
})


test_that("explicit columns warn when the object carries no annotation", {
  expect_warning(get.results(dpa.noinfo, 1, protein.info.columns = "gene.name"),
                 "No 'protein.info' table")
})


test_that("an annotation column named as a keyword is reported as ambiguous", {
  info <- info.3col
  colnames(info)[1] <- "all"
  obj <- add.protein.info(dpa.noinfo, info)
  expect_warning(get.results(obj, 1, protein.info.columns = "all"), "keyword")
})


test_that("the contrast attribute is preserved when the annotation is appended", {
  res <- get.results(dpa.info, 1, protein.info.columns = "all")
  expect_equal(attributes(res)$contrast, names(dpa.info@analyses.result.list)[1])
})
