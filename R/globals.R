# Declares the names used in non-standard evaluation (dplyr / ggplot2 aes) so that
# R CMD check does not flag them as "no visible binding for global variable".
# This is a plain top-level call (not roxygen), so it does not affect NAMESPACE.

utils::globalVariables(c(
  # --- PCoA plotting internals ---
  "axis.id", "Proportion.of.Variance", "Cumulative.Proportion", "Broken.stick",
  "loading.x.scaled", "loading.y.scaled", "variable", "shape", "color",
  "observed", "reconstructed.2D",

  # --- pre-existing NSE column names ---
  ".id", "AvgP", "FDR", "Group", "ID", "Intensity.log2", "NES", "P.Value", "PC",
  "PCs", "RMSE", "Sample", "Var1", "Var2", "X..i..", "adj.P.Val", "bFDR",
  "basemean.log2", "column.id", "contrast", "correlation", "correlation.coeff",
  "dataset", "diff.status", "expected.values", "fill", "fraction", "gene",
  "geneList", "geneSet_size", "group", "i", "imputation.method", "lfq.A", "lfq.B",
  "log2.FoldChange_bait.vs.control", "log2.Fold_group1.vs.group2",
  "log2.mean.group1", "log2.mean.group2", "log2FC", "mean.score", "n.NAs.in.row",
  "n.tot.NAs", "n.total", "p.value", "padj", "parameter", "percentage",
  "presence", "processing.time", "prot.id", "pval", "rank.cor", "rank.fc",
  "rank.stat", "ranking.score", "runningScore", "score", "setSize", "shared",
  "statistic", "status", "value", "x", "xmax", "xmin", "y", "ymax", "ymin"
))
