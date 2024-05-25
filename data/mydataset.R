unimputed.counts = readRDS(url("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/unimputed.counts.rds"))

sample.config = data.frame(column.id = c("Sample_A", "Sample_B", "Sample_C",
                                         "Sample_D", "Sample_E", "Sample_F", "Sample_G", "Sample_H", "Sample_I",
                                         "Sample_J", "Sample_K", "Sample_L"),
                           sample.id = c("BCa_FBS_rep1",
                                         "BCa_6h.DMSO_rep1", "BCa_6h.10nM.E2_rep1", "BCa_FBS_rep2", "BCa_6h.DMSO_rep2",
                                         "BCa_6h.10nM.E2_rep2", "BCa_FBS_rep3", "BCa_6h.DMSO_rep3", "BCa_6h.10nM.E2_rep3",
                                         "BCa_FBS_rep4", "BCa_6h.DMSO_rep4", "BCa_6h.10nM.E2_rep4"),
                           cell = c("BCa",
                                    "BCa", "BCa", "BCa", "BCa", "BCa", "BCa", "BCa", "BCa", "BCa",
                                    "BCa", "BCa"),
                           condition = c("FBS", "6h.DMSO", "6h.10nM.E2",
                                         "FBS", "6h.DMSO", "6h.10nM.E2", "FBS", "6h.DMSO", "6h.10nM.E2",
                                         "FBS", "6h.DMSO", "6h.10nM.E2"),
                           combined.id = c("BCa_FBS", "BCa_6h.DMSO",
                                           "BCa_6h.10nM.E2", "BCa_FBS", "BCa_6h.DMSO", "BCa_6h.10nM.E2",
                                           "BCa_FBS", "BCa_6h.DMSO", "BCa_6h.10nM.E2", "BCa_FBS", "BCa_6h.DMSO",
                                           "BCa_6h.10nM.E2"),
                           replicate = c(rep(c("rep1", "rep2", "rep3", "rep4"), each = 3)))

#dpo.QN = readRDS(url("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/dpo.QN.rds"))
#dpo.imputed = readRDS(url("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/dpo.imputed.rds"))
#dpo.DE.results = readRDS(url("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/dpo.DE.results.rds"))

################################################################################
# Generate data files
usethis::use_data(unimputed.counts,
                  sample.config,
                  #dpo.QN,
                  #dpo.imputed,
                  #dpo.DE.results,
                  overwrite = T)
