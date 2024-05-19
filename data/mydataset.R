unimputed.counts = readRDS(url("https://sebastian-gregoricchio.github.io/Rseb/inst/extdata/qPCR_results_rep1.rds"))
sample.config = readRDS(url("https://sebastian-gregoricchio.github.io/Rseb/inst/extdata/qPCR_results_rep2.rds"))

################################################################################
# Generate data files
usethis::use_data(unimputed.counts,
                  sample.config,
                  overwrite = T)
