unimputed.counts = readRDS(url("https://sebastian-gregoricchio.github.io/Rseb/inst/extdata/unimputed.counts.rds"))
sample.config = readRDS(url("https://sebastian-gregoricchio.github.io/Rseb/inst/extdata/sample.config.rds"))

################################################################################
# Generate data files
usethis::use_data(unimputed.counts,
                  sample.config,
                  overwrite = T)
