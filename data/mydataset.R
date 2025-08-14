# unimputed.counts = readRDS(url("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/unimputed.counts.rds"))
# dpo.imputed.counts = readRDS(url("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/imputed.counts.rds"))
#
# corum_v4.1 = data.table::fread("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/corum_v4.1.tsv", data.table = F)
#
# corum_v5.0 = data.table::fread("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/corum_v5.0.tsv", data.table = F)
#
# sample.config = data.frame(column.id = c("Sample_A", "Sample_B", "Sample_C",
#                                          "Sample_D", "Sample_E", "Sample_F", "Sample_G", "Sample_H", "Sample_I",
#                                          "Sample_J", "Sample_K", "Sample_L"),
#                            sample.id = c("BCa_FBS_rep1",
#                                          "BCa_6h.DMSO_rep1", "BCa_6h.10nM.E2_rep1", "BCa_FBS_rep2", "BCa_6h.DMSO_rep2",
#                                          "BCa_6h.10nM.E2_rep2", "BCa_FBS_rep3", "BCa_6h.DMSO_rep3", "BCa_6h.10nM.E2_rep3",
#                                          "BCa_FBS_rep4", "BCa_6h.DMSO_rep4", "BCa_6h.10nM.E2_rep4"),
#                            cell = c("BCa",
#                                     "BCa", "BCa", "BCa", "BCa", "BCa", "BCa", "BCa", "BCa", "BCa",
#                                     "BCa", "BCa"),
#                            condition = c("FBS", "6h.DMSO", "6h.10nM.E2",
#                                          "FBS", "6h.DMSO", "6h.10nM.E2", "FBS", "6h.DMSO", "6h.10nM.E2",
#                                          "FBS", "6h.DMSO", "6h.10nM.E2"),
#                            combined.id = c("BCa_FBS", "BCa_6h.DMSO",
#                                            "BCa_6h.10nM.E2", "BCa_FBS", "BCa_6h.DMSO", "BCa_6h.10nM.E2",
#                                            "BCa_FBS", "BCa_6h.DMSO", "BCa_6h.10nM.E2", "BCa_FBS", "BCa_6h.DMSO",
#                                            "BCa_6h.10nM.E2"),
#                            replicate = c(rep(c("rep1", "rep2", "rep3", "rep4"), each = 3)))
#
# #dpo.QN = readRDS(url("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/dpo.QN.rds"))
# #dpo.imputed = readRDS(url("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/dpo.imputed.rds"))
# #dpo.DE.results = readRDS(url("https://sebastian-gregoricchio.github.io/DEprot/inst/extdata/dpo.DE.results.rds"))
#
#
#
# ###################
# ## TEST TOOLBOX
# test.unimputed.lfq = DEprot::unimputed.counts[1:51,]
# test.dpo.raw = DEprot::load.counts2(counts = test.unimputed.lfq, metadata = DEprot::sample.config, data.type = "raw", log.base = 2)
# test.dpo.raw = DEprot::remove.undetected.proteins(test.dpo.raw, which.data = "raw")
# test.dpo.norm = DEprot::normalize.counts(test.dpo.raw)
# test.dpo.imp.miss = DEprot::impute.counts(test.dpo.norm)
# test.diff.exp.limma = DEprot::diff.analyses.limma(DEprot.object = test.dpo.imp.miss,
#                                                   contrast.list = list(c("condition", "FBS", "6h.DMSO"),
#                                                                        c("condition", "6h.10nM.E2", "6h.DMSO")),
#                                                   include.rep.model = TRUE,
#                                                   replicate.column = "replicate")
#
# set.seed(1234)
# test.geneset = do.call(rbind,
#                        list(data.frame(gs_name = "set1",
#                                        gene_symbol = unique(sample(rownames(test.unimputed.lfq),size = 23))),
#                             data.frame(gs_name = "set2",
#                                        gene_symbol = unique(sample(rownames(test.unimputed.lfq),size = 15))),
#                             data.frame(gs_name = "set3",
#                                        gene_symbol = unique(sample(rownames(test.unimputed.lfq),size = 10)))
#                        ))
#
#
#
# test.gsea.results = DEprot::geneset.enrichment(DEprot.analyses.object = test.diff.exp.limma,
#                                                contrast = 1,
#                                                TERM2GENE = test.geneset, enrichment.type = "gsea",
#                                                pvalueCutoff = 1,
#                                                qvalueCutoff = 1)
#
# test.ora.results = DEprot::geneset.enrichment(DEprot.analyses.object = test.diff.exp.limma,
#                                               contrast = 1,
#                                               TERM2GENE = test.geneset,
#                                               enrichment.type = "ORA",
#                                               pvalueCutoff = 1,
#                                               qvalueCutoff = 1,
#                                               diff.status.category = "6h.DMSO")
#
#
#
# test.toolbox = list(unimputed.lfq = test.unimputed.lfq,
#                     dpo.raw = test.dpo.raw,
#                     dpo.norm = test.dpo.norm,
#                     dpo.imp = test.dpo.imp.miss,
#                     diff.exp.limma = test.diff.exp.limma,
#                     geneset = test.geneset,
#                     gsea.results = test.gsea.results,
#                     ora.results = test.ora.results)
#
#



# ################################################################################
# # Generate data files
# usethis::use_data(unimputed.counts,
#                   dpo.imputed.counts,
#                   sample.config,
#                   corum_v4.1,
#                   corum_v5.0,
#                   #dpo.QN,
#                   #dpo.imputed,
#                   #dpo.DE.results,
#                   test.toolbox
#                   overwrite = TRUE)
