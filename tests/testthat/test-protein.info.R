## Real example objects.
## Notice that they have been generated before the introduction of the 'protein.info'
## slot, hence they also serve to test the backward compatibility.
dpo <- DEprot::test.toolbox$dpo.imp
dpa <- DEprot::test.toolbox$diff.exp.limma
dpo@protein.info <- NULL
dpa@protein.info <- NULL

prot.ids <- rownames(dpo@imputed.counts)


## Minimal annotation: a single column built on the protein IDs of the counts
info.1col <- data.frame(gene.name = paste0("GENE", seq_along(prot.ids)),
                        row.names = prot.ids,
                        stringsAsFactors = FALSE)

## Annotation with more than one column, used for the column selection
info.3col <- data.frame(gene.name = paste0("GENE", seq_along(prot.ids)),
                        n.peptides = seq_along(prot.ids),
                        description = paste("Protein", seq_along(prot.ids)),
                        row.names = prot.ids,
                        stringsAsFactors = FALSE)


# ---------------------------------------------------------------------
# .check.protein.info: identifiers
# ---------------------------------------------------------------------

test_that("the protein IDs are read from the row names", {
  res <- DEprot:::.check.protein.info(info.1col, prot.ids)
  expect_true(is.data.frame(res))
  expect_equal(rownames(res), prot.ids)
  expect_equal(res$gene.name, info.1col$gene.name)
})


test_that("the protein IDs are read from a 'prot.id' column when row names are absent", {
  info <- data.frame(prot.id = prot.ids,
                     gene.name = info.1col$gene.name,
                     stringsAsFactors = FALSE)
  res <- DEprot:::.check.protein.info(info, prot.ids)
  expect_equal(rownames(res), prot.ids)
  expect_false("prot.id" %in% colnames(res))   # the ID column is not duplicated
  expect_equal(colnames(res), "gene.name")
})


test_that("the protein IDs are read from a custom column", {
  info <- data.frame(uniprot = prot.ids,
                     gene.name = info.1col$gene.name,
                     stringsAsFactors = FALSE)
  res <- DEprot:::.check.protein.info(info, prot.ids, id.column = "uniprot")
  expect_equal(rownames(res), prot.ids)
  expect_false("uniprot" %in% colnames(res))
})


test_that("a matrix is accepted and converted to a data.frame", {
  mat <- as.matrix(info.1col)
  res <- DEprot:::.check.protein.info(mat, prot.ids)
  expect_true(is.data.frame(res))
  expect_equal(rownames(res), prot.ids)
})


test_that("NULL returns NULL", {
  expect_null(DEprot:::.check.protein.info(NULL, prot.ids))
})


# ---------------------------------------------------------------------
# .check.protein.info: alignment
# ---------------------------------------------------------------------

test_that("an unsorted annotation is re-ordered on the counts", {
  shuffled <- info.1col[sample(seq_along(prot.ids)), , drop = FALSE]
  res <- DEprot:::.check.protein.info(shuffled, prot.ids)
  expect_equal(rownames(res), prot.ids)
  expect_equal(res$gene.name, info.1col$gene.name)   # values follow their own protein
})


test_that("proteins missing from the annotation are filled with NAs", {
  partial <- info.1col[1:10, , drop = FALSE]
  res <- suppressMessages(DEprot:::.check.protein.info(partial, prot.ids))
  expect_equal(nrow(res), length(prot.ids))
  expect_equal(rownames(res), prot.ids)
  expect_equal(sum(is.na(res$gene.name)), length(prot.ids) - 10)
  expect_false(any(is.na(res$gene.name[1:10])))
})


test_that("annotations of proteins absent from the counts are discarded", {
  extra <- rbind(info.1col,
                 data.frame(gene.name = "GENE.EXTRA", row.names = "protein.not.there"))
  res <- suppressMessages(DEprot:::.check.protein.info(extra, prot.ids))
  expect_equal(nrow(res), length(prot.ids))
  expect_false("protein.not.there" %in% rownames(res))
})


test_that("a message reports the unannotated and the discarded proteins", {
  expect_message(DEprot:::.check.protein.info(info.1col[1:10, , drop = FALSE], prot.ids),
                 "not annotated")
  extra <- rbind(info.1col,
                 data.frame(gene.name = "X", row.names = "protein.not.there"))
  expect_message(DEprot:::.check.protein.info(extra, prot.ids), "discarded")
})


test_that("a warning is returned when no protein ID matches", {
  info <- data.frame(gene.name = c("A", "B"),
                     row.names = c("wrong.1", "wrong.2"),
                     stringsAsFactors = FALSE)
  expect_warning(res <- DEprot:::.check.protein.info(info, prot.ids), "None of the protein IDs")
  expect_true(all(is.na(res$gene.name)))
})


# ---------------------------------------------------------------------
# .check.protein.info: errors
# ---------------------------------------------------------------------

test_that("duplicated protein IDs return an error", {
  info <- data.frame(prot.id = c(prot.ids[1], prot.ids[1]),
                     gene.name = c("A", "B"),
                     stringsAsFactors = FALSE)
  expect_error(DEprot:::.check.protein.info(info, prot.ids), "duplicated")
})


test_that("missing or empty protein IDs return an error", {
  info <- data.frame(prot.id = c(prot.ids[1], NA),
                     gene.name = c("A", "B"),
                     stringsAsFactors = FALSE)
  expect_error(DEprot:::.check.protein.info(info, prot.ids), "missing values")
})


test_that("a table without identifiers returns an error", {
  info <- data.frame(gene.name = info.1col$gene.name, stringsAsFactors = FALSE)
  expect_error(DEprot:::.check.protein.info(info, prot.ids), "no protein IDs")
})


test_that("a table containing only the identifiers returns an error", {
  info <- data.frame(prot.id = prot.ids, stringsAsFactors = FALSE)
  expect_error(DEprot:::.check.protein.info(info, prot.ids), "no annotation column")
})


test_that("a non-tabular input returns an error", {
  expect_error(DEprot:::.check.protein.info(seq_along(prot.ids), prot.ids), "data.frame")
})


test_that("a non-existing id.column returns an error", {
  expect_error(DEprot:::.check.protein.info(info.1col, prot.ids, id.column = "not.a.column"),
               "not present")
})


# ---------------------------------------------------------------------
# add.protein.info
# ---------------------------------------------------------------------

test_that("add.protein.info attaches the annotation to a DEprot object", {
  res <- add.protein.info(DEprot.object = dpo, protein.info = info.1col)
  expect_s4_class(res, "DEprot")
  expect_true(is.data.frame(res@protein.info))
  expect_equal(rownames(res@protein.info), prot.ids)
  expect_equal(colnames(res@protein.info), "gene.name")
})


test_that("add.protein.info works on DEprot.analyses objects as well", {
  res <- add.protein.info(DEprot.object = dpa, protein.info = info.1col)
  expect_s4_class(res, "DEprot.analyses")
  expect_equal(rownames(res@protein.info), prot.ids)
})


test_that("add.protein.info aligns an unsorted and incomplete annotation", {
  partial <- info.1col[sample(1:10), , drop = FALSE]
  res <- suppressMessages(add.protein.info(dpo, partial))
  expect_equal(rownames(res@protein.info), prot.ids)
  expect_equal(sum(!is.na(res@protein.info$gene.name)), 10)
})


test_that("add.protein.info does not overwrite an existing annotation by default", {
  res <- add.protein.info(dpo, info.1col)
  expect_error(add.protein.info(res, info.3col), "overwrite")
})


test_that("add.protein.info replaces the annotation when overwrite is TRUE", {
  res <- add.protein.info(dpo, info.1col)
  res <- add.protein.info(res, info.3col, overwrite = TRUE)
  expect_equal(colnames(res@protein.info), colnames(info.3col))
})


test_that("add.protein.info removes the annotation when protein.info is NULL", {
  dpo.test = dpo
  dpo.test@protein.info = info.1col

  info = NULL
  dpo.info = add.protein.info(DEprot.object =dpo.test,
                              protein.info = info)

  res = suppressMessages(add.protein.info(DEprot.object = dpo.info,
                                          protein.info = NULL))

  expect_null(res@protein.info)
  expect_message(add.protein.info(DEprot.object = dpo.info, protein.info = NULL))
  expect_message(add.protein.info(DEprot.object = dpo, protein.info = NULL))
})


test_that("add.protein.info accepts a custom identifier column", {
  info <- data.frame(uniprot = prot.ids,
                     gene.name = info.1col$gene.name,
                     stringsAsFactors = FALSE)
  res <- add.protein.info(dpo, info, id.column = "uniprot")
  expect_equal(rownames(res@protein.info), prot.ids)
})


test_that("add.protein.info returns errors on wrong inputs", {
  expect_error(add.protein.info(iris, info.1col), "DEprot")
  expect_error(add.protein.info(dpo), "protein.info")
})


# ---------------------------------------------------------------------
# get.protein.info
# ---------------------------------------------------------------------

test_that("get.protein.info returns the annotation table", {
  res <- add.protein.info(dpo, info.1col)
  info <- get.protein.info(res)
  expect_true(is.data.frame(info))
  expect_equal(rownames(info), prot.ids)
  expect_equal(info$gene.name, info.1col$gene.name)
})


test_that("get.protein.info works on DEprot.analyses objects as well", {
  res <- add.protein.info(dpa, info.1col)
  expect_equal(rownames(get.protein.info(res)), prot.ids)
})


test_that("get.protein.info returns NULL and a message when no annotation is stored", {
  expect_message(info <- get.protein.info(dpo), "No 'protein.info' table")
  expect_null(info)
})


test_that("get.protein.info does not work for objects of other classes", {
  expect_error(get.protein.info(iris), "DEprot")
  expect_error(get.protein.info(DEprot::test.toolbox$geneset), "DEprot")
})


# ---------------------------------------------------------------------
# Backward compatibility (objects without the slot)
# ---------------------------------------------------------------------

test_that("objects predating the slot are handled as objects without annotation", {
  # test.toolbox has been generated before the introduction of protein.info
  expect_null(DEprot:::.get.protein.info(dpo))
  expect_no_error(add.protein.info(dpo, info.1col))
})


# ---------------------------------------------------------------------
# Propagation of the annotation
# ---------------------------------------------------------------------

test_that("filter.proteins keeps the annotation aligned to the counts (keep)", {
  res <- add.protein.info(dpo, info.1col)
  kept <- prot.ids[1:10]
  res <- filter.proteins(res, proteins = kept, mode = "keep")
  expect_equal(rownames(res@protein.info), rownames(res@imputed.counts))
  expect_equal(rownames(res@protein.info), kept)
})


test_that("filter.proteins keeps the annotation aligned to the counts (remove)", {
  res <- add.protein.info(dpo, info.1col)
  removed <- prot.ids[1:10]
  res <- filter.proteins(res, proteins = removed, mode = "remove")
  expect_equal(rownames(res@protein.info), rownames(res@imputed.counts))
  expect_false(any(removed %in% rownames(res@protein.info)))
})


test_that("remove.undetected.proteins keeps the annotation aligned to the counts", {
  res <- add.protein.info(dpo, info.1col)
  res <- remove.undetected.proteins(res, min.n.samples = 3, which.data = "imputed")
  expect_equal(rownames(res@protein.info), rownames(res@imputed.counts))
})


test_that("the annotation survives the differential analyses", {
  res <- add.protein.info(dpo, info.1col)
  res <- diff.analyses.limma(DEprot.object = res,
                             contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                             include.rep.model = TRUE,
                             replicate.column = "replicate",
                             linear.FC.th = 1.2)
  expect_s4_class(res, "DEprot.analyses")
  expect_equal(rownames(res@protein.info), prot.ids)
})
