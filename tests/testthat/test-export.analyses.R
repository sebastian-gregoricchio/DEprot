test_that("export.analyse writes files correctly", {
  tmp_dir <- tempfile("test_output_")
  dir.create(tmp_dir)

  suppressMessages(suppressWarnings(DEprot::export.analyses(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast.subset = 1, output.folder = tmp_dir, verbose = TRUE)))

  files <- list.files(tmp_dir)
  expect_true(length(files) > 0)

  unlink(tmp_dir, recursive = TRUE)  # cleanup
})
