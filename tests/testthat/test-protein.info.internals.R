## ----------------------------------------------------------------------------------------
##  internal_functions.protein.info.R
##  '.check.protein.info' validates and re-orders an annotation table, '.append.protein.info'
##  binds it to a results table. Both are reached through 'add.protein.info' / 'get.results',
##  but most of their branches are error and collision handling that a normal call never hits.
## ----------------------------------------------------------------------------------------

info.table <-
  function(ids, extra = TRUE) {
    tb <- data.frame(gene.name = toupper(ids), stringsAsFactors = FALSE)
    if (isTRUE(extra)) {tb$description <- paste("protein", seq_along(ids))}
    rownames(tb) <- ids
    return(tb)
  }

prot.ids <- rownames(tb.dpo.imp@imputed.counts)



## ----------------------------------------------------------------------------------------
##  .get.protein.info / .get.protein.ids
## ----------------------------------------------------------------------------------------

test_that("the annotation and the protein IDs are read from any DEprot class", {
  expect_null(DEprot:::.get.protein.info(tb.dpo.imp))
  expect_equal(DEprot:::.get.protein.ids(tb.dpo.imp), prot.ids)
  expect_equal(DEprot:::.get.protein.ids(tb.limma), rownames(tb.limma@imputed.counts))
})


test_that("an object without any counts has no protein IDs", {
  bare <- tb.dpo.raw
  bare@raw.counts <- NULL
  bare@norm.counts <- NULL
  bare@random.counts <- NULL
  bare@imputed.counts <- NULL

  expect_null(DEprot:::.get.protein.ids(bare))
})



## ----------------------------------------------------------------------------------------
##  .check.protein.info
## ----------------------------------------------------------------------------------------

check.info <-
  function(protein.info, ids = prot.ids, ...) {
    DEprot:::.check.protein.info(protein.info = protein.info,
                                 protein.ids = ids,
                                 arg.name = "protein.info",
                                 ...)
  }


test_that("an empty annotation is read as no annotation at all", {
  expect_null(check.info(NULL))
  expect_null(check.info(list()))
  expect_null(check.info(NA))
})


test_that("a matrix is accepted and converted to a data.frame", {
  m <- matrix(toupper(prot.ids), ncol = 1,
              dimnames = list(prot.ids, "gene.name"))

  out <- check.info(m, verbose = FALSE)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), length(prot.ids))
  expect_equal(rownames(out), prot.ids)
})


test_that("an object that is neither a table nor a matrix is rejected", {
  expect_error(check.info(1:10))
  expect_error(check.info("not a table"))
})


test_that("a table without any column is read as no annotation at all", {
  ## a data.frame with zero columns has length 0, hence it is caught by the early exit
  ## rather than by the column check further down
  empty <- data.frame(row.names = prot.ids)
  expect_null(check.info(empty))
})


test_that("the protein IDs can be given in a dedicated column", {
  tb <- data.frame(my.id = prot.ids,
                   gene.name = toupper(prot.ids),
                   stringsAsFactors = FALSE)

  out <- check.info(tb, id.column = "my.id", verbose = FALSE)

  expect_equal(rownames(out), prot.ids)
  ## the ID column is consumed, not kept as an annotation column
  expect_false("my.id" %in% colnames(out))
})


test_that("a 'prot.id' column is used when the table has no row names", {
  tb <- data.frame(prot.id = prot.ids,
                   gene.name = toupper(prot.ids),
                   stringsAsFactors = FALSE)

  out <- check.info(tb, verbose = FALSE)

  expect_equal(rownames(out), prot.ids)
})


test_that("an ID column that does not exist is reported", {
  expect_error(check.info(info.table(prot.ids), id.column = "not.a.column"))
})


test_that("a table with no identifiable protein IDs is rejected", {
  tb <- data.frame(gene.name = toupper(prot.ids), stringsAsFactors = FALSE)
  expect_error(check.info(tb))
})


test_that("a table carrying only the IDs has nothing to annotate", {
  tb <- data.frame(prot.id = prot.ids, stringsAsFactors = FALSE)
  expect_error(check.info(tb))
})


test_that("missing and duplicated protein IDs are rejected", {
  na.ids <- info.table(prot.ids)
  rownames(na.ids) <- NULL
  na.ids$prot.id <- c(NA, prot.ids[-1])
  expect_error(check.info(na.ids))

  empty.ids <- na.ids
  empty.ids$prot.id <- c("", prot.ids[-1])
  expect_error(check.info(empty.ids))

  dup <- data.frame(prot.id = c(prot.ids[1], prot.ids),
                    gene.name = c("A", toupper(prot.ids)),
                    stringsAsFactors = FALSE)
  expect_error(check.info(dup))
})


test_that("the annotation is re-ordered on the proteins of the object", {
  shuffled <- info.table(rev(prot.ids))

  out <- check.info(shuffled, verbose = FALSE)

  expect_equal(rownames(out), prot.ids)
  expect_equal(out$gene.name, toupper(prot.ids))
})


test_that("the proteins missing from the annotation are filled with NA", {
  partial <- info.table(prot.ids[1:5])

  out <- suppressMessages(check.info(partial, verbose = TRUE))

  expect_equal(nrow(out), length(prot.ids))
  expect_true(all(is.na(out$gene.name[-(1:5)])))
})


test_that("the annotations of proteins absent from the object are discarded", {
  extended <- info.table(c(prot.ids, paste0("extra.", 1:3)))

  expect_message(check.info(extended, verbose = TRUE), "discarded")

  out <- suppressMessages(check.info(extended, verbose = TRUE))
  expect_equal(nrow(out), length(prot.ids))
  expect_false(any(grepl("extra", rownames(out))))
})


test_that("an annotation matching no protein raises a warning", {
  unrelated <- info.table(paste0("other.", 1:5))

  expect_warning(check.info(unrelated, verbose = FALSE))

  ## the table is kept, aligned on the proteins of the object: every value is therefore NA
  out <- suppressWarnings(check.info(unrelated, verbose = FALSE))
  expect_equal(nrow(out), length(prot.ids))
  expect_true(all(is.na(out$gene.name)))
})


test_that("the messages can be silenced", {
  expect_silent(invisible(check.info(info.table(prot.ids[1:5]), verbose = FALSE)))
})



## ----------------------------------------------------------------------------------------
##  .append.protein.info
## ----------------------------------------------------------------------------------------

annotated <- suppressMessages(add.protein.info(DEprot.object = tb.limma,
                                               protein.info = info.table(rownames(tb.limma@imputed.counts))))

results.tb <- get.results(DEprot.analyses.object = tb.limma, contrast = 1)


append.info <-
  function(data = results.tb, object = annotated, ...) {
    DEprot:::.append.protein.info(data = data, DEprot.object = object, ...)
  }


test_that("'none' returns the table untouched", {
  expect_identical(append.info(protein.info.columns = "none"), results.tb)
  ## a NULL is read as the 'none' keyword
  expect_identical(append.info(protein.info.columns = NULL), results.tb)
})


test_that("'all' appends the whole annotation", {
  out <- append.info(protein.info.columns = "all")

  expect_true(all(c("gene.name", "description") %in% colnames(out)))
  expect_equal(nrow(out), nrow(results.tb))
  ## the annotation follows the protein IDs, not the row order
  expect_equal(out$gene.name, toupper(out$prot.id))
})


test_that("a subset of columns can be appended", {
  out <- append.info(protein.info.columns = "gene.name")

  expect_true("gene.name" %in% colnames(out))
  expect_false("description" %in% colnames(out))
})


test_that("columns that do not exist in the annotation are reported", {
  expect_error(append.info(protein.info.columns = "not.a.column"))
})


test_that("'protein.info.columns' must be a character vector", {
  expect_error(append.info(protein.info.columns = 1))
  expect_error(append.info(protein.info.columns = character(0)))
})


test_that("a prefix can be added to the annotation columns", {
  out <- append.info(protein.info.columns = "all", protein.info.prefix = "info.")

  expect_true("info.gene.name" %in% colnames(out))
  expect_false("gene.name" %in% colnames(out))
})


test_that("the colliding column names are renamed with a warning", {
  ## an annotation column named as one of the results columns would silently shadow it
  colliding <- info.table(rownames(tb.limma@imputed.counts))
  colnames(colliding)[1] <- "padj"

  object <- suppressMessages(add.protein.info(DEprot.object = tb.limma,
                                              protein.info = colliding,
                                              overwrite = TRUE))

  expect_warning(append.info(object = object, protein.info.columns = "all"))

  out <- suppressWarnings(append.info(object = object, protein.info.columns = "all"))
  expect_true("padj.1" %in% colnames(out))
  ## the results column keeps its own values
  expect_equal(out$padj, results.tb$padj)
})


test_that("nothing is appended when the object carries no annotation", {
  expect_warning(append.info(object = tb.limma, protein.info.columns = "all"))
  expect_identical(suppressWarnings(append.info(object = tb.limma, protein.info.columns = "all")),
                   results.tb)
})


test_that("nothing is appended without an object to read the annotation from", {
  expect_warning(DEprot:::.append.protein.info(data = results.tb,
                                               DEprot.object = NULL,
                                               protein.info.columns = "all"))
})


test_that("get.results forwards the arguments to the appending helper", {
  out <- get.results(DEprot.analyses.object = annotated,
                     contrast = 1,
                     protein.info.columns = "all",
                     protein.info.prefix = "annot.")

  expect_true("annot.gene.name" %in% colnames(out))
  expect_error(get.results(DEprot.analyses.object = annotated,
                           contrast = 1,
                           protein.info.columns = "not.a.column"))
})
