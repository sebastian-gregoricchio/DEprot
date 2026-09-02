#' @import methods
NULL

#' @title DEprot class
#'
#' @slot metadata The data.frame corresponding to the metadata table describing the samples. Class: \code{"ANY"}.
#' @slot protein.info A data.frame (or \code{NULL}) containing extra information about the proteins (e.g., gene symbol, description, number of peptides). It has one row per protein, row-by-row aligned with the counts tables, and the same row names. It is kept synchronized with the counts by all the functions that add or remove proteins. Class: \code{"ANY"}.
#' @slot raw.counts Numeric matrix (rows: proteins, columns: samples) of the raw counts. Class: \code{"ANY"}.
#' @slot norm.counts Numeric matrix (rows: proteins, columns: samples) of the normalized counts. Class: \code{"ANY"}.
#' @slot random.counts Numeric matrix (rows: proteins, columns: samples) of the randomized counts. Class: \code{"ANY"}.
#' @slot imputed.counts Numeric matrix (rows: proteins, columns: samples) of the imputed counts. Class: \code{"ANY"}.
#' @slot log.base Numeric value indicating the base of the logarithm expressing the values in the loaded data. Class: \code{"ANY"}.
#' @slot log.transformed Logical value indicating whether the data are log-transformed or not. Class: \code{"logical"}.
#' @slot normalized Logical value indicating whether the data are normalized. Class: \code{"logical"}.
#' @slot normalization.method String (or any other class) value indicating the normalization method. Class: \code{"ANY"}. Class: \code{"ANY"}.
#' @slot randomized Logical value indicating whether the bottom distribution randomization has been applied. Class: \code{"logical"}.
#' @slot randomization.method List indicating the parameters used for the randomization (bottom distribution). Class" \code{"ANY"}.
#' @slot imputed Logical value indicating whether the data are imputed. Class: \code{"logical"}.
#' @slot imputation.method String (or any other class) value indicating the imputation method. Class: \code{"ANY"}.
#' @slot boxplot.raw Ggplot object showing the distribution of the raw values per sample. Class: \code{"ANY"}.
#' @slot boxplot.norm Ggplot object showing the distribution of the normalized values per sample. Class: \code{"ANY"}.
#' @slot boxplot.random Ggplot object showing the distribution of the randomized values per sample. Class: \code{"ANY"}.
#' @slot boxplot.imputed Ggplot object showing the distribution of the imputed values per sample. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot",
         slots = list(metadata = "ANY",
                      protein.info = "ANY",
                      raw.counts = "ANY",
                      norm.counts = "ANY",
                      random.counts = "ANY",
                      imputed.counts = "ANY",
                      log.base = "numeric",
                      log.transformed = "logical",
                      randomized = "logical",
                      randomization.method = "ANY",
                      normalized = "logical",
                      normalization.method = "ANY",
                      imputed = "logical",
                      imputation.method = "ANY",
                      boxplot.raw = "ANY",
                      boxplot.norm = "ANY",
                      boxplot.random = "ANY",
                      boxplot.imputed = "ANY"))






#' @title DEprot.analyses class
#'
#' @description Extension of the \linkS4class{DEprot} class: it stores the same data and
#'   pre-processing information, plus the results of the differential analyses.
#'   All the slots of \linkS4class{DEprot} are inherited and behave identically.
#'
#' @slot analyses.result.list List containing the differential results for each contrast. Class: \code{"ANY"}. The list contains the following elements:
#'      \describe{
#'        \item{\code{results}: }{a data.frame containing the results of the analyses; includes average expression of each group, basemean, foldchange, pvalue and p.adj, test statistic, degrees of freedom and differential.status}
#'        \item{\code{n.diff}: }{a summary table showing the number of proteins in each differential expression status (up/down/unresponsive, null)}
#'        \item{\code{PCA.data}: }{output of \link{perform.PCA} for the subset of samples analyzed in a specific contrast}
#'        \item{\code{PCA.plots}: }{combination of 3 plots: scatter PC1-vs-PC2, scatter PC2-vs-PC3, and cumulative bar plot}
#'        \item{\code{correlations}: }{combination of Pearson and Spearman correlation heatmaps (obtained by \link{plot.correlation.heatmap}) for the subset of samples analyzed in a specific contrast}
#'        \item{\code{volcano}: }{volcano plot showing the log2(FoldChange) x -log10(p.adjusted) of differential expression results; it can be regenerated using \link{plot.volcano}}
#'        \item{\code{MA.plot}: }{MA-plot showing the log2(basemean) x log2(FoldChange) of differential expression results; it can be regenerated using \link{plot.MA}}}
#' @slot contrasts List of contrasts. each contrast is a vector indicating, in the order: metadata.table.column - groupA - groupB; (groupA / group B). Class: \code{"ANY"}.
#' @slot differential.analyses.params List of parameters used to run the differential analyses (fold change thresholds, p-value threshold, p-adjustment method, etc.). Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.analyses",
         contains = "DEprot",
         slots = list(analyses.result.list = "ANY",
                      contrasts = "ANY",
                      differential.analyses.params = "ANY"))





#' @title DEprot.PCA class
#'
#' @slot PCA.metadata metadata of the samples used in the PCA (subset of the original \code{DEprot@@metadata}. Class: \code{"ANY"}.
#' @slot sample.subset vector containing the list of samples analyzed. Class: \code{"ANY"}.
#' @slot data.used vector indicating the type of counts used (imputed, normalized, raw). Class: \code{"ANY"}.
#' @slot prcomp object of class \code{prcomp} (or output from \code{pcaMethods::pca}, method = "nipals") corresponding to the full PCA output. Class: \code{"ANY"}.
#' @slot PCs data.frame combining the PC scores and the metadata table, useful for replotting. Class: \code{"ANY"}.
#' @slot importance statistical summary table for the PCA analyses per each PC. Class: \code{"ANY"}.
#' @slot cumulative.PC.plot ggplot object corresponding to out put of \code{plot.PC.cumulative} for this object. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.PCA",
         slots = list(PCA.metadata = "ANY",
                      sample.subset = "ANY",
                      data.used = "ANY",
                      prcomp = "ANY",
                      PCs = "ANY",
                      importance = "ANY",
                      cumulative.PC.plot = "ANY"))






#' @title DEprot.correlation class
#'
#' @slot heatmap ggplot object corresponding to the correlation heatmap. Class: \code{"ANY"}.
#' @slot corr.metadata metadata of the samples used in the correlation (subset of the original \code{DEprot@@metadata}. Class: \code{"ANY"}.
#' @slot sample.subset vector containing the list of samples analyzed. Class: \code{"ANY"}.
#' @slot data.used vector indicating the type of counts used (imputed, normalized, raw). Class: \code{"ANY"}.
#' @slot corr.matrix the correlation matrix on which the heatmap is base on. Class: \code{"ANY"}.
#' @slot distance object of class \code{dist} corresponding to the output of \code{as.dist(1 - correlation.matrix)}. Class: \code{"ANY"}.
#' @slot cluster \code{hclust} object generated by \code{hclust(d = as.dist(1 - correlation.matrix), method = clustering.method)}. Class: \code{"ANY"}.
#' @slot method String indicating the method used for the correlation (e.g., 'pearson', 'spearman', 'kendal'). Class: \code{"ANY"}.

#'
#' @export

setClass(Class = "DEprot.correlation",
         slots = list(heatmap = "ANY",
                      corr.metadata = "ANY",
                      sample.subset = "ANY",
                      data.used = "ANY",
                      corr.matrix = "ANY",
                      distance = "ANY",
                      cluster = "ANY",
                      method = "ANY"))





#' @title DEprot.PCoA class
#'
#' @slot PCoA.metadata metadata of the samples used in the PCoA (subset of the original \code{DEprot@@metadata}). Class: \code{"ANY"}.
#' @slot sample.subset vector containing the list of samples analyzed. Class: \code{"ANY"}.
#' @slot data.used vector indicating the type of counts used (imputed, normalized, raw), inherited from the \code{DEprot.correlation} object. Class: \code{"ANY"}.
#' @slot correlation.method string indicating the correlation method used to build the dissimilarities ('pearson', 'spearman', 'kendall'). Class: \code{"ANY"}.
#' @slot correlation.matrix the sample-by-sample correlation matrix the PCoA is based on. Class: \code{"ANY"}.
#' @slot distance.input object of class \code{dist} as stored in the input \code{DEprot.correlation} object, i.e. \code{as.dist(1 - correlation.matrix)}. Class: \code{"ANY"}.
#' @slot distance object of class \code{dist} effectively used for the ordination, after transformation and/or non-euclidean correction. Class: \code{"ANY"}.
#' @slot distance.method string describing how the dissimilarities used were obtained. Class: \code{"ANY"}.
#' @slot cmdscale full output of \code{stats::cmdscale} (\code{points}, \code{eig}, \code{GOF}, \code{ac}). Class: \code{"ANY"}.
#' @slot PCos data.frame combining the principal coordinates and the metadata table, useful for replotting. Class: \code{"ANY"}.
#' @slot eigenvalues numeric vector of all the eigenvalues returned by the ordination. Class: \code{"ANY"}.
#' @slot importance statistical summary table of the ordination per each principal coordinate (eigenvalue, proportion/percentage of variance, cumulative proportion, broken-stick expectation). Class: \code{"ANY"}.
#' @slot euclidean.diagnostics list summarizing the amount and the magnitude of negative eigenvalues, and the goodness-of-fit of the embedding. Class: \code{"ANY"}.
#' @slot axis.loadings data.frame of protein-vs-axis correlations (only if a \code{DEprot} object was provided to \link{perform.PCoA}, otherwise \code{NULL}). Class: \code{"ANY"}.
#' @slot scatter.plot ggplot object corresponding to the scatter of the two principal coordinates requested by the user. Class: \code{"ANY"}.
#' @slot scatter.plot.123 patchwork object combining the PCo1-vs-PCo2 and PCo3-vs-PCo2 scatters. Class: \code{"ANY"}.
#' @slot cumulative.PCo.plot ggplot object showing the proportion of variance and its cumulative curve for each principal coordinate. Class: \code{"ANY"}.
#' @slot shepard.plot ggplot object comparing the input dissimilarities with the distances reconstructed in the reduced space (embedding quality diagnostic). Class: \code{"ANY"}.
#' @slot parameters list of the parameters used to compute the ordination. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.PCoA",
         slots = list(PCoA.metadata = "ANY",
                      sample.subset = "ANY",
                      data.used = "ANY",
                      correlation.method = "ANY",
                      correlation.matrix = "ANY",
                      distance.input = "ANY",
                      distance = "ANY",
                      distance.method = "ANY",
                      cmdscale = "ANY",
                      PCos = "ANY",
                      eigenvalues = "ANY",
                      importance = "ANY",
                      euclidean.diagnostics = "ANY",
                      axis.loadings = "ANY",
                      scatter.plot = "ANY",
                      scatter.plot.123 = "ANY",
                      cumulative.PCo.plot = "ANY",
                      shepard.plot = "ANY",
                      parameters = "ANY"))






#' @title DEprot.upset class
#'
#' @slot upset Ggplot object corresponding to upset plot displaying the overlaps of differential proteins between contrasts. Class: \code{"ANY"}.
#' @slot obs.matrix Logical matrix indicating all the proteins that are differentially expressed at least in a contrast (rows). Columns indicate a specific contrast. The logical values indicate whether a protein is found differential in a specific contrast (column). Therefore, this table can be used to extract the proteins included in a specific overlap. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.upset",
         slots = list(upset = "ANY",
                      obs.matrix = "ANY"))






#' @title DEprot.contrast.heatmap class
#'
#' @slot heatmap Ggplot object corresponding to any heatmap generate by either \link{heatmap.contrasts}. Class: \code{"ANY"}.
#' @slot cluster The \code{hclust} object of the rows (proteins). Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.contrast.heatmap",
         slots = list(heatmap = "ANY",
                      cluster = "ANY"))


#' @title DEprot.counts.heatmap class
#'
#' @slot heatmap Ggplot object corresponding to any heatmap generate by either \link{heatmap.counts}. Class: \code{"ANY"}.
#' @slot row.cluster The \code{hclust} object of the rows (proteins). Class: \code{"ANY"}.
#' @slot column.cluster The \code{hclust} object of the columns (samples). Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.counts.heatmap",
         slots = list(heatmap = "ANY",
                      row.cluster = "ANY",
                      column.cluster = "ANY"))



#' @title DEprot.enrichResult class
#'
#' @slot enrichment.discovery the direct output from \code{clusterProfiler::GSEA} or \code{clusterProfiler::enricher} (GSEA and ORA, respectively). Class: \code{"ANY"}.
#' @slot protein.network a string plot showing protein networks (\code{clusterProfiler::cnetplot}). Class: \code{"ANY"}.
#' @slot pathway.network a list with clusters and string plot showing pathway/set networks (\code{aPEAR::enrichmentNetwork}). Class: \code{"ANY"}.
#' @slot NES.plot (GSEA only) a bar plot showing the NES scores for each significantly enriched geneSet. Class: \code{"ANY"}.
#' @slot dotplot_gene.ratio a dotplot showing the geneRatios for each significantly enriched geneSet. Class: \code{"ANY"}.
#' @slot dotplot_fold.enrichment (ORA only) a dotplot showing the foldEnrichment for each significantly enriched geneSet. Class: \code{"ANY"}.
#' @slot parameters a list containing the parameters used to run the analyses. Class: \code{"ANY"}.
#' @slot affinity.propagation results of representative elements for each geneset.
#'
#' @export

setClass(Class = "DEprot.enrichResult",
         slots = list(enrichment.discovery = "ANY",
                      protein.network = "ANY",
                      pathway.network = "ANY",
                      NES.plot = "ANY",
                      dotplot_gene.ratio = "ANY",
                      dotplot_fold.enrichment = "ANY",
                      parameters = "ANY",
                      affinity.propagation = "ANY"))




#' @title DEprot.pvalues class
#'
#' @slot pvalue.distribution ggplot object depicting the histogram of the distribution of the p-values of the differential expression test of a specific contrast. Class: \code{"ANY"}.
#' @slot padjusted.distribution ggplot object depicting the histogram of the distribution of the adjusted p-values of the differential expression test of a specific contrast. Class: \code{"ANY"}.
#' @slot pvalue.rank ggplot object depicting the curve of the ranked p-values. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.pvalues",
         slots = list(pvalue.distribution = "ANY",
                      padjusted.distribution = "ANY",
                      pvalue.rank = "ANY"))




#' @title DEprot.normality class
#'
#' @slot norm.statement Logical value indicating whether the samples are normally distributed (\code{TRUE}). Class: \code{"ANY"}.
#' @slot norm.AD.tests List of Anderson-Darling normality test results for each sample. Class: \code{"ANY"}.
#' @slot qqplots List of ggplots objects depicting the Q-Q plots for the Anderson-Darling normality test. Class: \code{"ANY"}.
#' @slot densities List of ggplots objects depicting the destiny distribution of the intensities overlapped to a theoretical normal distribution. Class: \code{"ANY"}.
#' @slot p.threshold Numeric value indicating the P-value threshold to define whether a sample passed the normality test. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.normality",
         slots = list(norm.statement = "ANY",
                      norm.AD.tests = "ANY",
                      qqplots = "ANY",
                      densities = "ANY",
                      p.threshold = "ANY"))





#' @title DEprot.RMSE class
#'
#' @slot original.DEprot.object Object of class \code{DEprot} used to compute the RMSE. Class: \code{"ANY"}.
#' @slot percentage.test Percentage of the total proteins that should be used to perform the comparisons. Class: \code{"ANY"}.
#' @slot seed Seed used for the randomization. Class: \code{"ANY"}.
#' @slot fraction.missing.values Fraction of missing values in the original table. Class: \code{"ANY"}.
#' @slot test.dataset Subset of the original table used for the comparisons. Class: \code{"ANY"}.
#' @slot imputed.objects List of the output of \link{impute.counts} (class \code{DEprot}) using the different imputation methods. Class: \code{"ANY"}.
#' @slot RMSE.tables List of data.frames, one per tested imputation method, containing the coordinates of the imputed values and including the following columns (Class: \code{"ANY"}.):
#'      \describe{
#'        \item{\code{row.id}: }{id of the row (protein)}
#'        \item{\code{col.id}: }{id of the column (sample)}
#'        \item{\code{expected.values}: }{value measured in the experiment}
#'        \item{\code{imputation.method}: }{id of the method used for the imputation}
#'        \item{\code{imputed.values}: }{value imputed by \code{DEprot}}
#'        \item{\code{residuals}: }{difference of the the values, \code{imputed - expected}}
#'        \item{\code{sq.residuals}: }{the squared value of the residuals}}
#' @slot RMSE.scores A table indicating the methods and the corresponding RMSE score. Class: \code{"ANY"}.
#' @slot correlation.plots A list of ggplot objects with the correlation between observed/expected and imputed values. Class: \code{"ANY"}.
#' @slot density.residuals A list of ggplot objects depicting the distribution of the residuals. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.RMSE",
         slots = list(original.DEprot.object = "ANY",
                      percentage.test = "ANY",
                      seed = "ANY",
                      fraction.missing.values = "ANY",
                      test.dataset = "ANY",
                      imputed.objects = "ANY",
                      RMSE.tables = "ANY",
                      RMSE.scores = "ANY",
                      correlation.plots = "ANY",
                      density.residuals = "ANY"))




#' @title DEprot.SAINTq class
#'
#' @description
#' S4 object storing the results of \code{\link{SAINTq}}: the SAINTq
#' interaction-scoring tables, the corresponding volcano plots, and the
#' parameters used to compute the scores.
#'
#' @slot scores A named \code{list} (one element per bait) of \code{data.frame}s
#'   containing the SAINTq scoring tables. Each table reports one row per
#'   bait-prey pair with the following columns:
#'   \itemize{
#'     \item \code{Bait}: name of the bait group (a value of the \code{metadata.column}) scored in this table; constant within each per-bait element of the list.
#'     \item \code{Control}: name of the control group (a value of the \code{metadata.column}).
#'     \item \code{Prey}: identifier of the prey protein (the row name of the LFQ count matrix).
#'     \item \code{n.rep}: number of bait replicates in which the prey was quantified (i.e. non-missing values).
#'     \item \code{AvgP}: \emph{average probability}, the main SAINTq score.
#'       It is the posterior probability of a true (bait-specific) interaction
#'       averaged across the bait replicates. Ranges in \code{[0, 1]}; higher
#'       values indicate higher confidence.
#'     \item \code{MaxP}: the maximum per-replicate posterior probability across the bait replicates (the score of the single best replicate).
#'     \item \code{log2.FoldChange_bait.vs.control}: enrichment of the prey in the bait over the control, on the linear scale, computed as \code{log2(avg.bait) - log2(avg.ctrl)}.
#'     \item \code{avg.bait}: mean \code{log2} intensity of the prey across the bait replicates.
#'     \item \code{avg.ctrl}: mean \code{log2} intensity of the prey across the control run (the background mean \eqn{\mu_F} of the model;
#'       replaced by the global \code{background} value for preys never detected in the control).
#'     \item \code{bFDR}: Bayesian false discovery rate of the interaction,
#'       derived from the \code{AvgP} values as the cumulative mean of \code{(1 - AvgP)}
#'       down the probability-ranked list. Lower values indicate higher confidence.
#'   }
#'   Class: \code{"ANY"}.
#'
#' @slot volcanoes A named \code{list} (one element per bait) of \code{ggplot}
#'   objects. Each is a volcano-style plot of the interactions for that bait,
#'   showing the \code{log2(FoldChange)} on the x-axis against \code{-log10(bFDR)}
#'   on the y-axis, with the size and color of the points encoding the
#'   \code{AvgP} score. Class: \code{"ANY"}.
#'
#' @slot parameters A \code{list} recording the settings and fitted quantities used for the scoring:
#'   \itemize{
#'     \item \code{prior.pi}: prior probability of a true interaction (\eqn{\pi_T});
#'       either estimated by expectation-maximization or fixed by the user.
#'     \item \code{fold}: fold change separating the true component from the
#'       background, such that \code{mu_T = mu_F + log2(fold)}.
#'     \item \code{delta.log2}: \code{log2(fold)}, i.e. the additive shift (in \code{log2} units)
#'       applied to the background mean to define the mean of the true component.
#'     \item \code{sd.scale}: scaling factor of the true-component standard
#'       deviation relative to the background one (\code{sigma_T = sd.scale * sigma_F}).
#'     \item \code{min.sd}: lower bound applied to the per-prey background standard deviation,
#'       preventing over-confident scores from preys with near-constant control intensities.
#'     \item \code{sigma.global}: global background standard deviation (median of the per-prey control standard deviations),
#'       used as a fallback for preys with too few control measurements.
#'     \item \code{background}: \code{log2} background intensity assigned to preys never detected in any control run.
#'     \item \code{which.data}: the DEprot count matrix that was scored, one of \code{"imputed"}, \code{"randomized"}, \code{"normalized"} or \code{"raw"}.
#'     \item \code{control}: name of the control group used as the background.
#'     \item \code{baits}: name(s) of the bait group(s) that were scored.
#'     \item \code{best.n.rep}: number of top-scoring bait replicates used when averaging the per-replicate posterior probabilities
#'       into \code{AvgP} (the "best R replicates" option of SAINTexpress). \code{NULL} means that all replicates of each bait were used.
#'   }
#'   Class: \code{"ANY"}.
#'
#' @seealso \code{\link{SAINTq}}
#'
#' @export

setClass(Class = "DEprot.SAINTq",
         slots = list(scores = "ANY",
                      volcanoes = "ANY",
                      parameters = "ANY"))




#' @title DEprot.missingness class
#'
#' @description
#' S4 object storing the results of \code{\link{missingness.diagnostic}}: the classification of the
#' missing values into MNAR-like (left-censored, to be replaced by \link{randomize.missing.values})
#' and MCAR-like (to be replaced by \link{impute.counts}), together with the corresponding
#' summary tables and diagnostic plots.
#'
#' @slot data.used String indicating the counts table used for the diagnostic ('raw', 'normalized', 'randomized', 'imputed'). Class: \code{"ANY"}.
#' @slot counts.available String vector indicating all the counts tables available in the input object. Class: \code{"ANY"}.
#' @slot metadata Data.frame corresponding to the metadata table of the samples analyzed. Class: \code{"ANY"}.
#' @slot group.column String indicating the metadata column used to define the groups of replicates. Class: \code{"ANY"}.
#' @slot missing.matrix Logical matrix (rows: proteins, columns: samples) indicating the missing values (\code{TRUE}). Class: \code{"ANY"}.
#' @slot imputation.map Character matrix (rows: proteins, columns: samples) indicating, for each value, the strategy
#'   that the \code{DEprot} double-imputation would apply: \code{"detected"} (measured value), \code{"MNAR"}
#'   (missing in at least \code{percentage.missing}\% of the replicates of a group; randomized using the bottom of the
#'   distribution), \code{"MCAR"} (sparse missing value; imputed). Class: \code{"ANY"}.
#' @slot protein.stats Data.frame with one row per protein reporting the number/frequency of missing values (globally and
#'   per group), the average intensity and the assigned \code{missing.class} ('complete', 'MCAR', 'MNAR', 'all.missing'). Class: \code{"ANY"}.
#' @slot sample.stats Data.frame with one row per sample reporting the number and percentage of missing values, split by class. Class: \code{"ANY"}.
#' @slot group.summary Data.frame with one row per group of replicates summarizing the missing values and the number of
#'   proteins per class within the group. Class: \code{"ANY"}.
#' @slot pattern.summary Data.frame with the number and percentage of proteins in each missing-value class. Class: \code{"ANY"}.
#' @slot global.stats List of global metrics: total percentage of missing values, fraction of missing values that are
#'   MNAR-like, estimated \code{LOD50} (intensity at which 50\% of the values are missing), slope and p-value of the
#'   logistic dropout model, p-value of the intensity shift between complete and incomplete proteins, and the intensity
#'   threshold corresponding to the bottom \code{tail.percentage}\% of the distribution. Class: \code{"ANY"}.
#' @slot dropout.model Object of class \code{glm} corresponding to the logistic dropout model
#'   (\code{cbind(n.missing, n.detected) ~ mean.intensity}). Class: \code{"ANY"}.
#' @slot plots List of ggplot objects: \code{detection.density}, \code{dropout.curve}, \code{missingness.heatmap},
#'   \code{missing.per.sample}, \code{detection.frequency}, \code{pattern.barplot}, \code{sample.similarity}, \code{upset}. Class: \code{"ANY"}.
#' @slot jaccard.matrix Numeric matrix of the sample-vs-sample Jaccard similarity of the detection patterns
#'   (\code{intersection / union} of the sets of detected proteins). Class: \code{"ANY"}.
#' @slot jaccard.cluster \code{hclust} object generated by
#'   \code{hclust(d = as.dist(1 - jaccard.matrix), method = cluster.method)}, used to order the samples
#'   and to draw the dendrogram of the \code{sample.similarity} heatmap. Class: \code{"ANY"}.
#' @slot contrast.stats List (one element per contrast, \code{NULL} if no \code{DEprot.analyses} object was provided) containing
#'   the contrast id, the samples of each group, the per-protein classification restricted to the contrast
#'   (including the \code{testable} column, \code{FALSE} for proteins missing in both groups) and the corresponding plots. Class: \code{"ANY"}.
#' @slot parameters List of the parameters used to run the diagnostic (counts used, group column, \code{percentage.missing},
#'   \code{tail.percentage} and the source of these parameters: user-defined or retrieved from the randomization). Class: \code{"ANY"}.
#'
#' @seealso \code{\link{missingness.diagnostic}}, \code{\link{randomize.missing.values}}, \code{\link{impute.counts}}
#'
#' @export

setClass(Class = "DEprot.missingness",
         slots = list(data.used = "ANY",
                      counts.available = "ANY",
                      metadata = "ANY",
                      group.column = "ANY",
                      missing.matrix = "ANY",
                      imputation.map = "ANY",
                      protein.stats = "ANY",
                      sample.stats = "ANY",
                      group.summary = "ANY",
                      pattern.summary = "ANY",
                      global.stats = "ANY",
                      dropout.model = "ANY",
                      plots = "ANY",
                      jaccard.matrix = "ANY",
                      jaccard.cluster = "ANY",
                      contrast.stats = "ANY",
                      parameters = "ANY"))




#' @title DEprot.outliers class
#'
#' @slot metrics data.frame collecting, for each sample, the three quality metrics, their robust Z-scores/p-values, the individual flags, the total number of flags and the final outlier call, merged with the metadata table. Class: \code{"ANY"}.
#' @slot outliers vector containing the samples called as outliers. Class: \code{"ANY"}.
#' @slot sample.subset vector containing the list of samples analyzed. Class: \code{"ANY"}.
#' @slot data.used vector indicating the type of counts used (imputed, normalized, raw). Class: \code{"ANY"}.
#' @slot correlation.method String indicating the method used for the correlation (e.g., 'pearson', 'spearman', 'kendall'). Class: \code{"ANY"}.
#' @slot correlation.matrix the sample-by-sample correlation matrix used to compute the median correlations, with the diagonal set to \code{NA}. Class: \code{"ANY"}.
#' @slot PCA the \code{DEprot.PCA} object generated internally and used to compute the Mahalanobis distances. Class: \code{"ANY"}.
#' @slot missingness.data.used vector indicating the type of counts used to quantify the missing values, \code{NA} when no unimputed table was available. Class: \code{"ANY"}.
#' @slot group.column String indicating the metadata column used to restrict the correlations to the replicates of the same group, \code{NULL} when all the samples were used. Class: \code{"ANY"}.
#' @slot metrics.available logical vector indicating which of the three metrics could effectively be computed. Class: \code{"ANY"}.
#' @slot parameters list of the thresholds used to flag the samples. Class: \code{"ANY"}.
#' @slot plot patchwork object combining the diagnostic plots of all the metrics available. Class: \code{"ANY"}.
#' @slot plot.list list of the individual ggplot objects, one per metric. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.outliers",
         slots = list(metrics = "ANY",
                      outliers = "ANY",
                      sample.subset = "ANY",
                      data.used = "ANY",
                      correlation.method = "ANY",
                      correlation.matrix = "ANY",
                      PCA = "ANY",
                      missingness.data.used = "ANY",
                      group.column = "ANY",
                      metrics.available = "ANY",
                      parameters = "ANY",
                      plot = "ANY",
                      plot.list = "ANY"))



#' @title DEprot.timecourse class
#'
#' @slot results data.frame with one row per protein: statistics, kinetic descriptors, cluster assignment and ranking. Class: \code{"ANY"}.
#' @slot fitted.curves named list of matrices (proteins x grid points), one element per group level ('all' when no group was used). Class: \code{"ANY"}.
#' @slot time.grid numeric vector of the prediction grid positions, expressed on the ORIGINAL time scale. Class: \code{"ANY"}.
#' @slot observed.means named list of matrices (proteins x timepoints) with the mean measured value at each timepoint. Class: \code{"ANY"}.
#' @slot timepoints numeric vector of the observed timepoints. Class: \code{"ANY"}.
#' @slot counts.used matrix of the counts effectively analyzed (only the samples kept). Class: \code{"ANY"}.
#' @slot tc.metadata metadata of the samples used in the time course (subset of the original \code{DEprot@@metadata}). Class: \code{"ANY"}.
#' @slot sample.subset vector containing the list of samples analyzed. Class: \code{"ANY"}.
#' @slot data.used vector indicating the type of counts used (imputed, normalized, raw). Class: \code{"ANY"}.
#' @slot design model matrix used by limma for the fit. Class: \code{"ANY"}.
#' @slot clusters list with the membership matrix, the centroids, the method, the number of clusters and the fuzzifier. \code{NULL} when the clustering was not performed. Class: \code{"ANY"}.
#' @slot params list of all the parameters used for the analysis. Class: \code{"ANY"}.
#' @slot profile.plot ggplot object corresponding to the output of \code{plot.timecourse.profiles} for this object. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.timecourse",
         slots = list(results = "ANY",
                      fitted.curves = "ANY",
                      time.grid = "ANY",
                      observed.means = "ANY",
                      timepoints = "ANY",
                      counts.used = "ANY",
                      tc.metadata = "ANY",
                      sample.subset = "ANY",
                      data.used = "ANY",
                      design = "ANY",
                      clusters = "ANY",
                      params = "ANY",
                      profile.plot = "ANY"))



#' @title DEprot.timecourse.enrichment class
#'
#' @slot results data.frame combining the enrichment results of all the clusters, with a 'cluster' column. Class: \code{"ANY"}.
#' @slot enrichment.per.cluster named list containing the \code{enrichResult} object returned by \code{clusterProfiler::enricher} for each cluster. Class: \code{"ANY"}.
#' @slot dotplot ggplot object showing, for each cluster, the enriched genesets; the number inside each dot corresponds to the count of proteins. Class: \code{"ANY"}.
#' @slot universe vector containing the background gene list used for the tests. Class: \code{"ANY"}.
#' @slot parameters a list containing the parameters used to run the analyses. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.timecourse.enrichment",
         slots = list(results = "ANY",
                      enrichment.per.cluster = "ANY",
                      dotplot = "ANY",
                      universe = "ANY",
                      parameters = "ANY"))





#' @title DEprot.power class
#'
#' @slot power.table Data.frame reporting, for each sample size tested: the per-test significance level required to control the FDR, the average power, the expected number of true and false discoveries and the expected FDR. Class: \code{"ANY"}.
#' @slot n.required Numeric value indicating the smallest number of samples per group reaching the target average power (\code{NA} when the target is not reached within the range tested). Class: \code{"ANY"}.
#' @slot effect.size Numeric vector of the standardized effect sizes (|d|) used for the estimation. Class: \code{"ANY"}.
#' @slot pi0 Numeric value indicating the proportion of non-differential proteins used for the estimation. Class: \code{"ANY"}.
#' @slot power.plot Ggplot object showing the average power as function of the number of samples per group. Class: \code{"ANY"}.
#' @slot discoveries.plot Ggplot object showing the expected number of true discoveries as function of the number of samples per group. Class: \code{"ANY"}.
#' @slot effect.size.plot Ggplot object showing the distribution of the standardized effect sizes used. Class: \code{"ANY"}.
#' @slot params List of the parameters used for the estimation. Class: \code{"ANY"}.
#'
#' @export

setClass(Class = "DEprot.power",
         slots = list(power.table = "ANY",
                      n.required = "ANY",
                      effect.size = "ANY",
                      pi0 = "ANY",
                      power.plot = "ANY",
                      discoveries.plot = "ANY",
                      effect.size.plot = "ANY",
                      params = "ANY"))





# ===================================================================================================================
# ===================================================================================================================
# ===================================================================================================================
#                                     ################# METHODS #################
# ===================================================================================================================
# ===================================================================================================================
# ===================================================================================================================

## ------------------------------------------------------------------------- ##
## Internal helpers for the analysis slots
## ------------------------------------------------------------------------- ##

## analyses.result.list, contrasts and differential.analyses.params exist only in
## DEprot.analyses objects: this returns the value they used to have in a plain
## DEprot object (NULL, or NA for the contrasts), so that the functions accepting
## both classes do not need to branch on the class before reading them
.deprot_analysis_slot = function(object,
                                 slot.name,
                                 default = NULL) {
  if (methods::is(object, "DEprot.analyses")) {
    methods::slot(object, slot.name)
  } else {
    default
  }
}


#' @title DEprot updateObject-method
#' @description Drops the analysis slots from the \code{DEprot} objects built with DEprot <= 2.1.0,
#'   where they were declared but never filled. Run it on the objects restored from an old
#'   \code{.rds} file: \code{dpo = updateObject(dpo)}.
#' @param object Object of class \code{DEprot}.
#' @param ... Not used.
#' @param verbose Logical value indicating whether a message should be printed when the object is updated. Default: \code{FALSE}.
#' @return An object of class \code{DEprot}.
#' @author Sebastian Gregoricchio
#' @export
setGeneric(name = "updateObject",
           def = function(object, ..., verbose = FALSE) {standardGeneric("updateObject")})


#' @rdname updateObject
#' @exportMethod updateObject
setMethod(
  f = "updateObject",
  signature = "DEprot",
  definition = function(object, ..., verbose = FALSE) {
    ## objects written by an older class definition keep their slots in the
    ## attributes: rebuilding from the current definition is enough to drop them
    if (methods::is(object, "DEprot.analyses")) {
      return(object)
    }

    obsolete = c("analyses.result.list", "contrasts", "differential.analyses.params")
    stored = attributes(object)

    if (!any(obsolete %in% names(stored))) {
      return(object)
    }

    if (verbose == TRUE) {
      message("The object was built with DEprot <= 2.1.0: the unused analysis slots are removed.")
    }

    keep = setdiff(methods::slotNames("DEprot"), obsolete)
    args = lapply(keep, function(x) {stored[[x]]})
    names(args) = keep

    do.call(methods::new, c(list(Class = "DEprot"), args))
  }
)


#' @title DEprot show-method
#' @param object Object of class \code{DEprot}
#' @export
setMethod(
  f = "show",
  signature = "DEprot",
  definition = function(object) {

    # slot -> label, in display order
    count.slots <- c(raw.counts     = "raw",
                     norm.counts    = "normalized",
                     random.counts  = "randomized",
                     imputed.counts = "imputed")

    available <- count.slots[!vapply(
      names(count.slots),
      function(s) .deprot_slot_is_empty(methods::slot(object, s)),
      logical(1)
    )]

    # every populated table shares the same dimensions; use the first available
    tb <- if (length(available)) methods::slot(object, names(available)[1]) else NULL

    log.txt <- if (is.na(object@log.base)) "none (linear)" else paste0("log", object@log.base)

    cat("DEprot object:")
    cat("\n           Samples: ", if (is.null(tb)) 0L else ncol(tb))
    cat("\n          Proteins: ", if (is.null(tb)) 0L else nrow(tb))
    cat("\n  Counts available: ", paste(available, collapse = ", "))
    cat("\nLog transformation: ", log.txt)
    cat("\n  Metadata columns: ", paste(colnames(object@metadata), collapse = ", "), "\n")

    info <- .get.protein.info(object)
    if (!is.null(info)) {
      cat("      Protein info: ", paste(colnames(info), collapse = ", "), "\n")
    }
  }
)


#' @title DEprot plot-method
#' @param x Object of class \code{DEprot}.
#' @param y Not used.
#' @param ... Further arguments passed to \code{patchwork::wrap_plots()} (e.g. \code{guides}, \code{widths}, \code{heights}, \code{design}). They must be named.
#' @param ncol,nrow The dimensions of the grid to create - if both are NULL (default) it will use the same logic as \code{facet_wrap()} to set the dimensions.
#' @keywords internal
#' @importFrom patchwork wrap_plots
#' @export
setMethod(
  f = "plot",
  signature = "DEprot",
  definition = function(x, y, ..., ncol = NULL, nrow = NULL) {
    boxplot.slots <- c("boxplot.raw", "boxplot.norm",
                       "boxplot.random", "boxplot.imputed")
    plot.list <- lapply(boxplot.slots, function(s) methods::slot(x, s))
    plot.list <- Filter(Negate(.deprot_slot_is_empty), plot.list)
    if (length(plot.list) == 0L) {
      message("No boxplots are available in this DEprot object.")
      return(invisible(NULL))
    }

    ## wrap_plots() reads its unnamed arguments as further plots to assemble: an unnamed
    ## value passed through '...' would silently become a panel instead of an option
    extra <- list(...)
    if (length(extra) > 0L && (is.null(names(extra)) || any(names(extra) == ""))) {
      stop("Only named arguments of 'patchwork::wrap_plots()' can be passed through '...'.",
           call. = FALSE)
    }

    p <- patchwork::wrap_plots(plot.list, ncol = ncol, nrow = nrow, ...)
    print(p)
    invisible(p)
  }
)



#' @title DEprot.analyses show-method
#' @param object Object of class \code{DEprot.analyses}
#' @export
setMethod(f = "show",
          signature = "DEprot.analyses",
          definition =
            function(object) {

              ## Summarise differential expression results
              recap = data.frame()

              for (i in 1:length(object@analyses.result.list)) {
                recap = rbind(recap,
                              cbind(data.frame(n.contrast = i,
                                               contrast.id = rep(paste0(object@contrasts[[i]]$metadata.column, ": ", object@contrasts[[i]]$var.1, " vs ", object@contrasts[[i]]$var.2), 4),
                                               group.factor = rep(object@contrasts[[i]]$metadata.column, 4),
                                               group1 = rep(object@contrasts[[i]]$var.1, 4),
                                               group2 = rep(object@contrasts[[i]]$var.2, 4),
                                               paired.test = object@contrasts[[i]]$paired.test),
                                    object@analyses.result.list[[i]]$n.diff))

              }

              cat("DEprot.analyses object:")
              cat("\n          Counts used: ", object@differential.analyses.params$counts.used)
              cat("\nFold Change threshold: ", object@differential.analyses.params$linear.FC.th, "(linear)")
              cat("\nFC unresponsive range:  [", object@differential.analyses.params$linear.FC.unresp.range[1],",",object@differential.analyses.params$linear.FC.unresp.range[2],"] (linear)", sep = "")
              cat("\n       padj threshold: ", object@differential.analyses.params$padj.th, "(linear)")
              cat("\n          padj method: ", object@differential.analyses.params$padj.method)
              cat("\n")
              cat("\n")
              cat("\nDifferential results summary:\n")

              print(recap)
            }#end definition
) #end method



#' @title DEprot.analyses summary-method
#' @param object Object of class \code{DEprot.analyses}
#' @export
setMethod(f = "summary",
          signature = "DEprot.analyses",
          definition =
            function(object) {

              ## Summarise differential expression results
              recap = data.frame()

              for (i in 1:length(object@analyses.result.list)) {
                recap = rbind(recap,
                              cbind(data.frame(contrast.id = rep(paste0(object@contrasts[[i]]$metadata.column, ": ", object@contrasts[[i]]$var.1, " vs ", object@contrasts[[i]]$var.2), 4),
                                               group.factor = rep(object@contrasts[[i]]$metadata.column, 4),
                                               group1 = rep(object@contrasts[[i]]$var.1, 4),
                                               group2 = rep(object@contrasts[[i]]$var.2, 4),
                                               paired.test = object@contrasts[[i]]$paired.test),
                                    object@analyses.result.list[[i]]$n.diff))

              }
              recap
            } #end definition
) #end method



#' @title DEprot.analyses plot-method
#' @param x Object of class \code{DEprot.analyses}
#' @param y Not used.
#' @param ... Further arguments passed to \code{patchwork::wrap_plots()} (e.g. \code{guides}, \code{widths}, \code{heights}, \code{design}). They must be named, and they are used only when several contrasts are assembled.
#' @param plot.type String indicating which plots need to be summarized: 'volcano', 'MA', 'correlation', 'PCA'. Default: \code{"volcano"}.
#' @param label.top.n Single integer value (or \code{NULL}) indicating the number of top differentially expressed proteins to label automatically. Differential proteins (up- or down-regulated) are ranked by \code{-log10(padj) * abs(log2FC)} and the top \code{N} are labeled; if fewer than \code{N} differential proteins are available, all of them are labeled. When \code{use.uncorrected.pvalue = TRUE}, the uncorrected p-value is used in the ranking instead of the adjusted one. Any IDs provided through \code{dot.labels} are added to the automatically selected ones, and \code{labels.in.boxes} together with the other label aesthetics apply to them as well. Default: \code{NULL} (no automatic labels).
#' @param ncol,nrow The dimensions of the grid to create - if both are NULL (default) it will use the same logic as \code{facet_wrap()} to set the dimensions.
#' @keywords internal
#' @importFrom patchwork wrap_plots
#' @export
setMethod(f = "plot",
          signature = "DEprot.analyses",
          definition =
            function(x, y, ..., plot.type = "volcano", label.top.n = NULL, ncol = NULL, nrow = NULL) {

              ## wrap_plots() reads its unnamed arguments as further plots to assemble: an
              ## unnamed value passed through '...' would silently become a panel
              extra = list(...)
              if (length(extra) > 0L && (is.null(names(extra)) || any(names(extra) == ""))) {
                stop("Only named arguments of 'patchwork::wrap_plots()' can be passed through '...'.", call. = FALSE)
              }

              if (tolower(plot.type) %in% c("v", "volcano", "volcanos")) {
                plots = lapply(seq_along(x$analyses.result.list),
                               function(i) {
                                 plot.volcano(DEprot.analyses.object = x, contrast = i, label.top.n = label.top.n)
                               }) }
              else if (tolower(plot.type) %in% c("m", "ma", "mas")) {
                plots = lapply(seq_along(x$analyses.result.list),
                               function(i) {
                                 plot.MA(DEprot.analyses.object = x, contrast = i, label.top.n = label.top.n)
                               }) }
              else if (tolower(plot.type) %in% c("c", "cor", "corr", "cors", "corrs", "correlation", "correlations")) {
                plots = lapply(seq_along(x$analyses.result.list),
                               function(i) {
                                 x$analyses.result.list[[i]]$correlations
                               }) }
              else if (tolower(plot.type) %in% c("p", "pc", "pcs", "pca", "pcas")) {
                plots = lapply(seq_along(x$analyses.result.list),
                               function(i) {
                                 x$analyses.result.list[[i]]$PCA.plots
                               }) }
              else {
                stop(paste0("The 'plot.type' value ('", plot.type, "') is not recognized.\n",
                            "       Please indicate one among: 'volcano', 'MA', 'correlation', 'PCA'."), call. = FALSE)
              }

              ## the correlations and the PCA panels are not stored by every differential
              ## function: an empty element would be assembled as a blank panel
              plots = Filter(Negate(.deprot_slot_is_empty), plots)

              if (length(plots) == 0L) {
                message(paste0("No '", plot.type, "' plot is available in this DEprot.analyses object."))
                return(invisible(NULL))
              }

              if (length(plots) > 1) {
                p = patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow, ...)
                print(p)
                invisible(p) }
              else {
                print(plots[[1]])
                invisible(plots[[1]])
              }
            })



#' @title DEprot.PCA show-method
#' @param object Object of class \code{DEprot.PCA}
#' @export
setMethod(f = "show",
          signature = "DEprot.PCA",
          definition =
            function(object) {
              cat("DEprot.PCA object:")
              cat("\n  Samples analyzed: ", object@sample.subset)
              cat("\n         Data used: ", paste(object@data.used, "(log2)\n"))
              cat("\n")
              cat("\n")
              cat("PCs vectors:\n")
              print(object@PCs)})



#' @title DEprot.PCoA show-method
#' @param object Object of class \code{DEprot.PCoA}
#' @export
setMethod(f = "show",
          signature = "DEprot.PCoA",
          definition =
            function(object) {
              cat("DEprot.PCoA object:")
              cat("\n  Samples analyzed: ", object@sample.subset)
              cat("\n         Data used: ", paste(object@data.used, "(log2)"))
              cat("\n       Correlation: ", object@correlation.method)
              cat("\n      Dissimilarity:", object@distance.method)
              cat("\n")

              if (!is.null(object@euclidean.diagnostics)) {
                cat("\n  Negative eigenvalues: ",
                    object@euclidean.diagnostics$n.negative.eigenvalues,
                    paste0("(max |relative| magnitude = ",
                           round(object@euclidean.diagnostics$max.relative.negative.eigenvalue, 4), ")"))
                cat("\n     Goodness-of-fit: ",
                    paste0(round(object@euclidean.diagnostics$GOF[1] * 100, 1), "%"))
                cat("\n")
              }

              cat("\n")
              cat("Principal coordinates:\n")
              print(object@PCos)})



#' @title DEprot.correlation show-method
#' @param object Object of class \code{DEprot.correlation}
#' @export
setMethod(f = "show",
          signature = "DEprot.correlation",
          definition = function(object) {print(object@heatmap)})



#' @title DEprot.upset show-method
#' @param object Object of class \code{DEprot.upset}
#' @export
setMethod(f = "show",
          signature = "DEprot.upset",
          definition = function(object) {print(object@upset)})


#' @title DEprot.contrast.heatmap show-method
#' @param object Object of class \code{DEprot.contrast.heatmap}
#' @export
setMethod(f = "show",
          signature = "DEprot.contrast.heatmap",
          definition = function(object) {print(object@heatmap)})


#' @title DEprot.counts.heatmap show-method
#' @param object Object of class \code{DEprot.counts.heatmap}
#' @export
setMethod(f = "show",
          signature = "DEprot.counts.heatmap",
          definition = function(object) {print(object@heatmap)})



#' @title DEprot.enrichResult show-method
#' @param object Object of class \code{DEprot.enrichResult}
#' @export
setMethod(f = "show",
          signature = "DEprot.enrichResult",
          definition = function(object) {print(object@enrichment.discovery@result)})


#' @title DEprot.pvalues show-method
#' @param object Object of class \code{DEprot.pvalues}
#' @import patchwork
#' @export
setMethod(f = "show",
          signature = "DEprot.pvalues",
          definition =
            function(object) {
              #require(patchwork, quietly = TRUE)
              plot = (object@pvalue.distribution / object@padjusted.distribution) | object@pvalue.rank
              print(plot)
            })


#' @title DEprot.normality show-method
#' @param object Object of class \code{DEprot.normality}
#' @export
setMethod(f = "show",
          signature = "DEprot.normality",
          definition =
            function(object) {
              if (object@norm.statement == TRUE) {
                message("All samples display a normal distribution.")

              } else {
                normality = sapply(X = object@norm.AD.tests, FUN = function(x){x$p.value < object@p.threshold}, USE.NAMES = TRUE)
                message(paste0("The following samples do not display a normal distribution: ",
                               paste0(names(normality)[isFALSE(normality)], collapse = ", "), "."))
              }
            })



#' @title DEprot.normality plot-method
#' @param x Object of class \code{DEprot.normality}
#' @param y Not used.
#' @param ... Further arguments passed to \code{patchwork::wrap_plots()} (e.g. \code{guides}, \code{widths}, \code{heights}, \code{design}). They must be named, and they are used only when several contrasts are assembled.
#' @param n.samples NUmber of samples to display (in order by metadata table). Default: \code{NULL} (all samples are shown).
#' @keywords internal
#' @importFrom patchwork wrap_plots
#' @export
setMethod(f = "plot",
          signature = "DEprot.normality",
          definition =
            function(x, y, ..., n.samples = NULL) {

              if (!is.null(n.samples)) {
                n = min(c(length(x@densities), n.samples))
              } else {
                n = length(x@densities)
              }

              plot = patchwork::wrap_plots(c(x@qqplots[1:n], x@densities[1:n]), byrow = FALSE, ncol = 2, ...)
              print(plot)
              invisible(plot)
            })



#' @title DEprot.RMSE show-method
#' @param object Object of class \code{DEprot.RMSE}
#' @importFrom patchwork wrap_plots
#' @export
setMethod(f = "show",
          signature = "DEprot.RMSE",
          definition =
            function(object) {
              plot = patchwork::wrap_plots(object@correlation.plots)
              print(plot)
              invisible(plot)
            })


#' @title DEprot.RMSE summary-method
#' @param object Object of class \code{DEprot.RMSE}
#' @export
setMethod(f = "summary",
          signature = "DEprot.RMSE",
          definition = function(object) {object@RMSE.scores}
) #end method



#' @title DEprot.SAINTq show-method
#' @param object Object of class \code{DEprot.SAINTq}
#' @importFrom patchwork wrap_plots
#' @export
setMethod(f = "show",
          signature = "DEprot.SAINTq",
          definition =
            function(object) {
              plot = patchwork::wrap_plots(object@volcanoes)
              print(plot)
            })



#' @title DEprot.SAINTq summary-method
#' @param object Object of class \code{DEprot.SAINTq}
#' @export
setMethod(f = "summary",
          signature = "DEprot.SAINTq",
          definition = function(object) {
            tb = do.call(rbind, object@scores)
            rownames(tb) = NULL
            return(tb)}
) #end method





#' @title DEprot.missingness show-method
#' @param object Object of class \code{DEprot.missingness}
#' @export
setMethod(f = "show",
          signature = "DEprot.missingness",
          definition =
            function(object) {

              gs = object@global.stats
              ps = object@pattern.summary

              cat("DEprot.missingness object:")
              cat("\n           Counts used: ", object@data.used,
                  paste0("(available: ", paste(object@counts.available, collapse = ", "), ")"))
              cat("\n              Proteins: ", gs$n.proteins)
              cat("\n               Samples: ", gs$n.samples)
              cat("\n          Group column: ", object@group.column)
              cat("\n     MNAR defined when: ", paste0(">= ", object@parameters$percentage.missing, "% missing values within a group"))
              cat("\n            Parameters: ", object@parameters$parameters.source)
              cat("\n")
              cat("\n        Missing values: ", paste0(gs$n.missing, "/", gs$n.values,
                                                       " (", round(gs$perc.missing, 2), "%)"))
              cat("\n        of which MNAR: ", paste0(gs$n.cells.MNAR, " (", round(gs$perc.missing.MNAR, 1), "% of the missing values)"))
              cat("\n        of which MCAR: ", paste0(gs$n.cells.MCAR, " (", round(100 - gs$perc.missing.MNAR, 1), "% of the missing values)"))
              cat("\n")
              cat("\n  Estimated LOD50: ", ifelse(is.na(gs$LOD50), yes = "not estimable", no = round(gs$LOD50, 3)))
              cat("\n    Dropout slope: ", ifelse(is.na(gs$dropout.slope), yes = "NA",
                                                  no = paste0(round(gs$dropout.slope, 3),
                                                              " (p = ", format.pval(gs$dropout.pvalue, digits = 3), ")")))
              cat("\n  Intensity shift: ", ifelse(is.na(gs$intensity.shift.pvalue), yes = "NA",
                                                  no = paste0("incomplete - complete = ",
                                                              round(gs$median.intensity.incomplete - gs$median.intensity.complete, 3),
                                                              " (p = ", format.pval(gs$intensity.shift.pvalue, digits = 3), ")")))
              cat("\n")

              ## interpretation helper
              mnar.evidence = (!is.na(gs$dropout.slope) && gs$dropout.slope < 0 &&
                                 !is.na(gs$dropout.pvalue) && gs$dropout.pvalue < 0.05)

              cat("\n  Interpretation: ",
                  ifelse(mnar.evidence,
                         yes = "the missingness is intensity-dependent (MNAR/left-censoring is dominant).",
                         no = "no clear intensity-dependence of the missingness (MCAR-like)."))
              cat("\n")
              cat("\n")
              cat("Proteins per missing-value class:\n")
              print(ps)

              if (!is.null(object@contrast.stats)) {
                cat("\nContrast-level diagnostics available for: ",
                    paste(names(object@contrast.stats), collapse = ", "), "\n")
              }
            }#end definition
) #end method



#' @title DEprot.missingness summary-method
#' @param object Object of class \code{DEprot.missingness}
#' @export
setMethod(f = "summary",
          signature = "DEprot.missingness",
          definition =
            function(object) {

              if (is.null(object@contrast.stats)) {
                return(object@group.summary)
              }

              recap =
                do.call(rbind,
                        lapply(X = names(object@contrast.stats),
                               FUN = function(n) {
                                 x = object@contrast.stats[[n]]
                                 cbind(data.frame(contrast.id = n,
                                                  group.factor = x$metadata.column,
                                                  group1 = x$var.1,
                                                  group2 = x$var.2,
                                                  stringsAsFactors = FALSE),
                                       x$summary)
                               }))
              rownames(recap) = NULL
              return(recap)
            }#end definition
) #end method



#' @title DEprot.missingness plot-method
#' @param x Object of class \code{DEprot.missingness}.
#' @param y Not used.
#' @param ... Further arguments passed to \code{patchwork::wrap_plots()} (e.g. \code{guides}, \code{widths}, \code{heights}, \code{design}). They must be named, and they are used only when several contrasts are assembled.
#' @param plot.type String indicating which plot(s) should be returned: 'summary' (default, combination of the four main
#' diagnostics), 'density', 'dropout', 'heatmap', 'samples', 'frequency', 'classes', 'similarity', 'upset', 'contrasts', 'all'.
#' @param contrast Index or name of a contrast (only when \code{plot.type = "contrasts"}). Default: \code{NULL} (all contrasts).
#' @param ncol,nrow The dimensions of the grid to create - if both are NULL (default) it will use the same logic as
#' \code{facet_wrap()} to set the dimensions.
#' @keywords internal
#' @importFrom patchwork wrap_plots
#' @export
setMethod(f = "plot",
          signature = "DEprot.missingness",
          definition =
            function(x, y, ..., plot.type = "summary", contrast = NULL, ncol = NULL, nrow = NULL) {

              p = x@plots

              plot.list =
                switch(EXPR = tolower(plot.type),
                       "summary" = list(p$detection.density, p$dropout.curve, p$pattern.barplot, p$missing.per.sample),
                       "density" = list(p$detection.density),
                       "dropout" = list(p$dropout.curve),
                       "heatmap" = list(p$missingness.heatmap),
                       "samples" = list(p$missing.per.sample),
                       "frequency" = list(p$detection.frequency),
                       "classes" = list(p$pattern.barplot),
                       "similarity" = list(p$sample.similarity),
                       "upset" = list(p$upset),
                       "all" = p,
                       "contrasts" = NULL,
                       stop(paste0("The 'plot.type' is not recognized. Choose among: 'summary', 'density', 'dropout', ",
                                   "'heatmap', 'samples', 'frequency', 'classes', 'similarity', 'upset', 'contrasts', 'all'.")))

              if (tolower(plot.type) == "contrasts") {
                if (is.null(x@contrast.stats)) {
                  message("No contrast-level diagnostic is available in this object.")
                  return(invisible(NULL))
                }

                selected = if (is.null(contrast)) {names(x@contrast.stats)} else {names(x@contrast.stats)[contrast]}
                plot.list = unlist(lapply(x@contrast.stats[selected],
                                          function(i) {list(i$plots$pattern.barplot, i$plots$detection.density)}),
                                   recursive = FALSE)
              }

              plot.list = Filter(Negate(is.null), plot.list)

              if (length(plot.list) == 0) {
                message("No plot is available for the requested 'plot.type'.")
                return(invisible(NULL))
              }

              if (length(plot.list) == 1) {
                print(plot.list[[1]])
                return(invisible(plot.list[[1]]))
              }

              combined = patchwork::wrap_plots(plot.list, ncol = ncol, nrow = nrow, ...)
              print(combined)
              invisible(combined)
            }#end definition
) #end method





#' @title DEprot.outliers show-method
#' @param object Object of class \code{DEprot.outliers}
#' @export
setMethod(f = "show",
          signature = "DEprot.outliers",
          definition =
            function(object) {
              cat("DEprot.outliers object:")
              cat("\n  Samples analyzed: ", length(object@sample.subset))
              cat("\n         Data used: ", paste0(object@data.used, " (log2)"))
              cat("\n Metrics available: ", paste(names(object@metrics.available)[object@metrics.available == TRUE], collapse = ", "))
              cat("\n    Flags required: ", object@parameters$min.flags)
              cat("\n          Outliers: ", ifelse(length(object@outliers) > 0,
                                                   yes = paste(object@outliers, collapse = ", "),
                                                   no = "none"))
              cat("\n")
              cat("\n")
              cat("Sample metrics:\n")
              print(object@metrics[,c("column.id", "group", "median.correlation", "correlation.z",
                                      "mahalanobis.distance", "mahalanobis.padj",
                                      "missing.rate", "missingness.z", "n.flags", "outlier")])})



#' @title DEprot.outliers plot-method
#' @param x Object of class \code{DEprot.outliers}
#' @param y Not used.
#' @param ... Not used.
#' @keywords internal
#' @export
setMethod(f = "plot",
          signature = "DEprot.outliers",
          definition =
            function(x, y, ...) {
              print(x@plot)
              invisible(x@plot)
            })





#' @title DEprot.timecourse show-method
#'
#' @param object An object of class \code{DEprot.timecourse}.
#'
#' @return Prints a summary of the time-course analysis.
#'
#' @export

setMethod(f = "show",
          signature = "DEprot.timecourse",
          definition =
            function(object) {

              p = object@params

              cat("DEprot.timecourse object:\n")
              cat("  Counts used:       ", object@data.used, "\n", sep = "")
              cat("  Time column:       ", p$time.column, " (transform: ", p$time.transform, ")\n", sep = "")
              cat("  Timepoints:        ", paste(object@timepoints, collapse = ", "),
                  "  [n = ", length(object@timepoints), "]\n", sep = "")
              cat("  Model:             ", p$model, " (df = ", p$spline.df, ")\n", sep = "")

              if (!is.null(p$group.column)) {
                cat("  Group column:      ", p$group.column, " (group x time interaction tested)\n", sep = "")
              }

              cat("  Proteins tested:   ", nrow(object@results), "\n", sep = "")
              cat("  Trending proteins: ", sum(object@results$trend.status == "trending"),
                  " (padj < ", p$padj.th, ", amplitude > ", p$log2.amplitude.th, ")\n", sep = "")

              if (!is.null(object@clusters)) {
                cat("  Clusters:          ", object@clusters$k, " (", object@clusters$method,
                    ifelse(is.na(object@clusters$fuzzifier), "", paste0(", m = ", object@clusters$fuzzifier)),
                    ")\n", sep = "")
              }
            })






#' @title DEprot.timecourse summary-method
#'
#' @param object An object of class \code{DEprot.timecourse}.
#' @param ... Not used.
#'
#' @return A data.frame summarizing each cluster: size, dominant shape, median amplitude, median time to half-maximum, median peak time and best-ranked protein. The clusters are returned sorted by their median t.half, meaning in the order in which they respond.
#'
#' @importFrom stats median na.omit
#'
#' @export

setMethod(f = "summary",
          signature = "DEprot.timecourse",
          definition =
            function(object, ...) {

              res = object@results[object@results$trend.status == "trending",]

              ## nothing to summarize if the clustering was not performed
              if (nrow(res) == 0 | is.null(object@clusters)) {return(res)}

              group.levels = object@params$group.levels

              ## with several groups the descriptors are suffixed by the group name
              res$.t.half = .tc.descriptor.column(results = res, descriptor = "t.half", group.levels = group.levels)
              res$.peak.time = .tc.descriptor.column(results = res, descriptor = "peak.time", group.levels = group.levels)

              summary.tb =
                do.call(rbind,
                        lapply(sort(unique(stats::na.omit(res$cluster))),
                               function(k) {
                                 sub = res[which(res$cluster == k),]
                                 shapes = sort(table(sub$trend.shape), decreasing = TRUE)

                                 data.frame(cluster = k,
                                            n = nrow(sub),
                                            dominant.shape = names(shapes)[1],
                                            fraction.shape = round(shapes[1] / nrow(sub), 3),
                                            median.amplitude = round(stats::median(sub$amplitude), 3),
                                            median.t.half = round(stats::median(sub$.t.half, na.rm = TRUE), 3),
                                            median.peak.time = round(stats::median(sub$.peak.time), 3),
                                            top.protein = sub$prot.id[which.min(sub$rank.in.cluster)],
                                            row.names = NULL,
                                            stringsAsFactors = FALSE)
                               }))

              ## sorted by median t.half: the clusters appear in the order in which they
              ## respond, which is the reading order for a sequential process
              summary.tb = summary.tb[order(summary.tb$median.t.half),]
              rownames(summary.tb) = NULL

              return(summary.tb)
            })





#' @title DEprot.timecourse plot-method
#'
#' @param x An object of class \code{DEprot.timecourse}.
#' @param y Not used.
#' @param ... Not used.
#'
#' @return The cluster-profile plot stored in the object.
#'
#' @export

setMethod(f = "plot",
          signature = "DEprot.timecourse",
          definition =
            function(x, y, ...) {
              return(x@profile.plot)
            })



#' @title DEprot.timecourse.enrichment show-method
#'
#' @param object An object of class \code{DEprot.timecourse.enrichment}.
#'
#' @return Prints a summary of the enrichment analyses.
#'
#' @export

setMethod(f = "show",
          signature = "DEprot.timecourse.enrichment",
          definition =
            function(object) {

              p = object@parameters

              cat("DEprot.timecourse.enrichment object:\n")
              cat("  Clusters tested:  ", paste(p$clusters, collapse = ", "), "\n", sep = "")
              cat("  Universe size:    ", length(object@universe), " proteins\n", sep = "")
              cat("  Thresholds:       padj < ", p$pvalueCutoff, " (", p$pAdjustMethod,
                  "), qvalue < ", p$qvalueCutoff, "\n", sep = "")

              if (is.null(object@results)) {
                cat("  Enriched genesets: none\n")
              } else {
                signif = object@results[object@results$p.adjust <= p$pvalueCutoff,]
                cat("  Enriched genesets: ", nrow(signif), " (in ",
                    length(unique(signif$cluster)), " cluster(s))\n", sep = "")
              }
            })




#' @title DEprot.timecourse.enrichment plot-method
#'
#' @param x An object of class \code{DEprot.timecourse.enrichment}.
#' @param y Not used.
#' @param ... Not used.
#'
#' @return The dotplot stored in the object.
#'
#' @export

setMethod(f = "plot",
          signature = "DEprot.timecourse.enrichment",
          definition =
            function(x, y, ...) {
              return(x@dotplot)
            })





#' @title DEprot.power show-method
#' @param object Object of class \code{DEprot.power}.
#' @importFrom patchwork wrap_plots
#' @export
setMethod(f = "show",
          signature = "DEprot.power",
          definition =
            function(object) {
              cat(paste0("        Contrast | ", object@params$contrast, ifelse(isTRUE(object@params$paired.test), " (paired)", ""), "\n"))
              cat(paste0("     Counts used | ", object@params$counts.used, " (", object@params$stat.test, ")\n"))
              cat(paste0(" Proteins tested | ", object@params$m, " (m1 = ", object@params$m1, ", pi0 = ", round(object@pi0, 3), ")\n"))
              cat(paste0("     Effect size | ", object@params$effect.size.mode, ", median |d| = ", round(median(object@effect.size, na.rm = TRUE), 3), "\n"))
              cat(paste0("             FDR | ", object@params$fdr, "\n"))
              cat(paste0(" Current n/group | ", object@params$n.current, "\n"))
              cat(paste0("Required n/group | ", ifelse(is.finite(object@n.required),
                                                       paste0(object@n.required, " (average power ", object@params$target.power, ")"),
                                                       paste0("not reached within the range tested")), "\n\n"))
            })




#' @title DEprot.power plot-method
#'
#' @param x An object of class \code{DEprot.power}.
#' @param y Not used.
#' @param ... Passed to \code{patchwork::wrap_plots()}.
#'
#' @importFrom patchwork wrap_plots
#'
#' @return The combination of all the plots available.
#'
#' @export

setMethod(f = "plot",
          signature = "DEprot.power",
          definition =
            function(x, y, ...) {
              print(patchwork::wrap_plots(x@power.plot, x@discoveries.plot, x@effect.size.plot, ...))
            })






## ================================================================================================
## ================================================================================================
## ================================================================================================
## ================================================================================================
##  Instance-aware slot completion for all DEprot S4 classes.
##  Only slots holding a real value (not NULL / not NA / not empty)
##  are offered after `$` and `@`.
## ================================================================================================
## ================================================================================================
## ================================================================================================
## ================================================================================================


# All DEprot S4 classes
.deprot_classes <- c("DEprot", "DEprot.analyses", "DEprot.PCA", "DEprot.PCoA",
                     "DEprot.correlation", "DEprot.upset",
                     "DEprot.contrast.heatmap", "DEprot.counts.heatmap",
                     "DEprot.enrichResult", "DEprot.pvalues",
                     "DEprot.normality", "DEprot.RMSE", "DEprot.SAINTq",
                     "DEprot.missingness", "DEprot.outliers",
                     "DEprot.timecourse", "DEprot.timecourse.enrichment",
                     "DEprot.power")

# A slot counts as "empty" if it is NULL, zero-length, or a single NA.
# Logical flags holding FALSE (e.g. `normalized`) are NOT empty and stay visible.
.deprot_slot_is_empty <- function(v) {
  is.null(v) ||
    length(v) == 0L ||
    (is.atomic(v) && length(v) == 1L && is.na(v))
}

# Completion worker shared by $ and @ (signature matches .DollarNames/.AtNames)
.deprot_complete_slots <- function(x, pattern = "") {
  sn <- methods::slotNames(x)
  full <- sn[!vapply(sn,
                     function(s) .deprot_slot_is_empty(methods::slot(x, s)),
                     logical(1))]
  grep(pattern, full, value = TRUE)
}







# S4 has no `$` by default, so define it (plus `$<-`) for every class,
# otherwise the `$` completions would error the moment they are evaluated.

#' @title Slot access for DEprot classes
#'
#' @description Instance-aware slot access and replacement for all the DEprot S4 classes.
#' These methods allow the slots of a DEprot object to be accessed and modified with the
#' \code{$} operator (in addition to the standard \code{@} operator), enabling slot-name
#' auto-completion in interactive sessions.
#'
#' @param x An object of one of the DEprot classes: \code{DEprot}, \code{DEprot.analyses},
#' \code{DEprot.PCA}, \code{DEprot.PCoA}, \code{DEprot.correlation}, \code{DEprot.counts.heatmap},
#' \code{DEprot.contrast.heatmap}, \code{DEprot.upset}, \code{DEprot.normality}, \code{DEprot.pvalues},
#' \code{DEprot.enrichResult}, \code{DEprot.RMSE}, \code{DEprot.SAINTq}, \code{DEprot.missingness},
#' \code{DEprot.outliers}, \code{DEprot.timecourse}, \code{DEprot.timecourse.enrichment}, \code{DEprot.power}.
#' @param name Name of the slot to access or replace.
#' @param value The value to assign to the slot (only for \code{$<-}).
#'
#' @return For \code{$}, the content of the corresponding slot. For \code{$<-}, the object with the updated slot.
#'
#' @name DEprot-dollar-methods
#' @rdname DEprot-dollar-methods
#'
#' @author Sebastian Gregoricchio
NULL




#' Dollar-style access for DEprot objects
#' @name DEprot-dollar
#' @keywords internal
#' @exportMethod $
#' @exportMethod $<-
NULL


for (.cls in .deprot_classes) {
  setMethod("$", .cls, function(x, name) methods::slot(x, name))
  setMethod("$<-", .cls, function(x, name, value) {
    methods::slot(x, name) <- value
    x
  })
}
rm(.cls)


