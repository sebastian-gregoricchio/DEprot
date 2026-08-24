---
title: "changeLog"
---

# DEprot [<img src="https://sebastian-gregoricchio.github.io/DEprot/DEprot_logo.png" align="right" height = 150/>](https://sebastian-gregoricchio.github.io/DEprot)
![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/DEprot)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/DEprot?style=social)](https://github.com/sebastian-gregoricchio/DEprot/fork)


#### [v2.1.0](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/2.1.0) - August 24<sup>th</sup> 2026
- added `estimate.power` function to compute sample size and power estimation
- updated overview vignette, manual and tests accordingly
- added `heatmap.counts.anno()`: the heatmap of `heatmap.counts()`, with the same data selection (`which.data`, `contrast`, `top.n`, `sample.subset`, `protein.subset`, `group.by.metadata.column`, `scale`), drawn with `ComplexHeatmap` instead of `ggplot2`. It returns the usual `DEprot.counts.heatmap` object, whose `heatmap` slot contains a `Heatmap` object, hence it can be customized, drawn and concatenated with the standard `ComplexHeatmap` syntax
- `ComplexHeatmap` and `circlize` added to the dependencies
- added function `generate.mm` to generate materials and methods paragraph
- updated overview vignette, manual and tests accordingly

<br>

#### [v2.0.0](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/2.0.0) - August 7<sup>th</sup> 2026
- added the protein.info slot to the DEprot and DEprot.analyses objects: an optional annotation table with one row per protein (gene symbol, description, number of peptides, etc.), kept row-by-row aligned with the counts
- `load.counts2()` (and `load.counts()`) accept a protein.info table at loading; the IDs can be given as row names, in a prot.id column or in any column indicated by protein.info.id.column. The table is re-ordered on the counts: unannotated proteins are filled with NA and annotations of proteins absent from the counts are discarded
- added `add.protein.info()` to attach, replace or remove the annotation of an object built previously
- added `get.protein.info()` to extract the annotation table from an object, in the same fashion as get.metadata()
- the annotation is propagated by `filter.proteins()`, `remove.undetected.proteins()`, `harmonize.batches()`, `filter.samples()` and by all the `diff.analyses*()` functions
- `get.results()` gains protein.info.columns and protein.info.prefix: the annotation is appended only on request, using the keywords "none" (default) and "all" or a vector of column names. The columns are bound to the results by protein ID (not by position) and colliding names are made unique
- `export.analyses()` writes the annotated results tables
- objects created with previous versions of the package remain usable: the missing slot is interpreted as an absent annotation
- added `missingness.diagnostic()`: classifies the missing values as MNAR-like or MCAR-like following the same rule of the double-imputation strategy, estimates the limit of detection by a logistic dropout model (`LOD50`), and returns the diagnostic plots (detection *vs* abundance densities, dropout curve, imputation map, Jaccard similarity of the detection patterns with dendrogram, per-sample and per-class summaries, UpSet of the MNAR-like patterns)
- the parameters of `randomize.missing.values` (`group.column`, `percentage.missing`, `tail.percentage`) are retrieved automatically when the randomization has already been run; the counts are selected at the lowest level available (raw > normalized > randomized > imputed)
- contrast-level diagnostics are computed when a `DEprot.analyses` object is provided, including the directional MNAR classification and the `testable` flag for the proteins missing in both groups
- added class `DEprot.missingness` with `show`, `summary` and `plot` methods
- fixed the `ComplexUpset` themes in `plot.upset()`, which were not valid anymore since `ggplot2` v4.0.0
- added time-course analyses: `analyze.timecourse()` returns an object of the new class `DEprot.timecourse` (with `show`, `summary` and `plot` methods). Time is handled as a **numeric** covariate: a natural-spline basis is fitted by `limma` and all the time coefficients are tested jointly by a moderated F-test, resulting in a single test per protein and therefore in no contrast to define
- `analyze.timecourse()`: the spline degrees of freedom are capped automatically at `n.timepoints - 2` (a saturated spline being equivalent to treating the time as a factor), and a linear trend is fitted below 4 timepoints; the `time.transform` option (`log2`, `log10`, `log1p`, `sqrt`) re-spaces log-spaced designs so that the last timepoint does not dominate the fit
- `analyze.timecourse()`: the results table reports the kinetic descriptors of each fitted trajectory (`amplitude`, `initial.slope`, `peak.time`, `trend.shape`: monotone/transient/complex) together with the cluster assignment, the membership and the ranking score
- `analyze.timecourse()`: the trending proteins are soft-clustered (c-means, or PAM) on the Z-scored **fitted** curves rather than on the raw timepoint means, which makes the clustering robust to unequal time spacing and to missing values; the fuzzifier and the number of clusters are estimated from the data when not provided
- added `rank.timecourse()` to re-rank the proteins, globally and within each cluster, without refitting the model, and `get.timecourse.results()` to retrieve the results with cluster/top-N subsetting
- added `plot.timecourse.protein()` (measured points, mean ± SEM and fitted trajectory of individual proteins) and `plot.timecourse.profiles()` (one panel per cluster, lines colored by membership). Both accept `values = "counts"`, `"log2FC"` (relative to `reference.time`) or `"zscore"`
- added `heatmap.timecourse()`: heatmap of the trending proteins with the rows split by cluster, sortable by `rank`, `membership`, `peak.time`, `amplitude` or `hclust`, and displaying either the measured means or the smooth fitted trajectories
- added `timecourse.enrichment()`: over-representation analyses run independently on each cluster, using all the quantified proteins as universe, returning an object of the new class `DEprot.timecourse.enrichment` (with `show` and `plot` methods). The dotplot shows the enrichment as dot size, the significance as color, and the number of proteins of the geneset written inside each dot
- `splines`, `e1071` and `cluster` added to the dependencies
- added `combine.enrichments()`: merges a named list of enrichments in a single dotplot, with the discoveries on the x-axis and the genesets on the y-axis; the dots are sized by fold enrichment or gene ratio, colored by significance, and carry the protein count inside
- added `divergent.enrichment()`: plots two enrichments back-to-back, the second one with the sign inverted, to show together the two sides of the same contrast (*e.g.* ORA of the up- and down-regulated proteins); bar length can be `FoldEnrichment`, `GeneRatio`, `Count`, `NES` or `padj`
- both functions accept a mix of `DEprot.enrichResult`, `DEprot.timecourse.enrichment`, `enrichResult`/`gseaResult` (clusterProfiler) objects and plain result tables; a `DEprot.timecourse.enrichment` is expanded into one discovery per cluster
- added the internal helpers `.parse.ratio()`, `.get.enrichment.table()` and `.collect.enrichments()`, which harmonize the columns of GSEA and ORA results (for a GSEA the leading edge is used as equivalent of the ORA `Count`); `timecourse.enrichment()` now uses the exported `.parse.ratio()` instead of its local copy
- added `export.external()`: converts a `DEprot`/`DEprot.analyses` object into a `SummarizedExperiment`, `QFeatures`, `MSnSet`, `limma::EList` or a plain list of tables, with the shortcuts `as.SummarizedExperiment()`, `as.QFeatures()` and `as.MSnSet()`
- all the count matrices available are exported as assays, the metadata become the column annotation, and `protein.info` together with the differential results become the row annotation
- the parameters that the destination classes cannot store (log base, normalization/imputation methods, contrasts, thresholds) are kept in the metadata of the exported object; `keep.object = TRUE` allows a lossless round-trip
- `SummarizedExperiment`, `QFeatures` and `MSnbase` are optional dependencies: requested only when needed and never installed without confirmation
- added `detect.outliers()`: automatic per-sample QC flagging combining median inter-sample correlation, robust Mahalanobis distance in PC space and missing rate; returns a `DEprot.outliers` object (with `show` and `plot` methods) whose `outliers` slot can be passed directly to `filter.samples()`
- added `diff.analyses.proDA()`: differential analyses on left-censored data using proDA. The missing values are not replaced but modelled as observations below the detection limit, so that the uncertainty on the undetected values propagates into lfcSE instead of being compressed by the imputation. It requires counts that still contain the missing values (which.data = "normalized") and it is the natural continuation of a strongly-MNAR verdict of missingness.diagnostic(), whose output can be passed directly through missingness.object
the results table of `diff.analyses.proDA()` keeps the columns of the other differential functions and adds n.detected.<group> and n.approx; each contrast stores a proDA.fit element with the fitted model and the per-sample dropout curves
- `filter.samples()` recognizes the new engine and re-runs the contrasts with the stored parameters
- `diff.analyses.limma()` now returns the lfdr column when padj.method = "fdrtool", as already stated in the manual proDA added to the dependencies
- updated overview vignette, manual and tests accordingly
- added time-course vignette

<br>

#### [v1.3.1](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/1.3.1) - August 3<sup>rd</sup> 2026
- added `import.external()`: builds a `DEprot` object directly from DIA-NN, Spectronaut, FragPipe, MaxQuant and Proteome Discoverer reports, with the shortcuts `read.diann()`, `read.diann.matrix()`, `read.spectronaut()`, `read.fragpipe()` and `read.maxquant()`
- added `import.msstats()` for the summarized objects of `MSstats` (label-free) and `MSstatsTMT` (isobaric); the metadata are reconstructed from the object annotation when not provided
- `iq` (MaxLFQ summarization) and `nanoparquet` (DIA-NN `.parquet` reports) are optional dependencies: they are requested only when needed and never installed without confirmation
- updated vignette and manual accordingly


<br>

#### [v1.3.0](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/1.3.0) - July 30<sup>th</sup> 2026
- `plot.volcano` and `plot.MA` can automatically plot the top N differential proteins
- `export.report` will plot the top.n proteins in the volcano
- updated vignette and manual accordingly
- added Principal Coordinate Analysis: `perform.PCoA()` with `plot.PCoA.scatter()`, `plot.PCoA.scatter.123()`, `plot.PCoA.cumulative()` and `plot.PCoA.biplot()`, computed from a `DEprot.correlation` object
- updated correlation class, now contains also a slot with the correlation method used
- re-arrangement of internal functions
- update vignette and manual


<br>

#### [v1.2.0](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/1.2.0) - July 12<sup>th</sup> 2026
- `expression.boxplot` function can now show pair-wise comparisons
- `diff.analyses.limma` can use `fdrtool` to adjust the p-values
- added `plot` method for objects of class `DEprot`, `DEprot.analyses`, `DEprot.normality`
- added function `filter.samples`
- added function `export.report` to export an HTML QC-report
- bug-fixed in `check.pvalues` when prolfqua tables are provided
- updated vignette and manual accordingly


<br>

#### [v1.1.3](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/1.1.3) - July 8<sup>th</sup> 2026
- `harmonize.batches()`: added `algorithm` (ComBat/limma), `ComBat.mode`, and `block`; limma helps on sparse designs where ComBat hits a singular matrix
- `harmonize.batches()`: drop uncorrectable proteins with a warning, error on an empty result, and fix the `algorithm`/`block` defaults
- `load.counts2` automatically converts the counts into log2 transformed.
- `diff.analyses`, `diff.analyses.limma` and `diff.analyses.prolfqua` now return, for each group, the standard deviation (`sd.<group>`) and standard error of the mean (`sem.<group>`) of the log2 expression values, plus `lfcSE`, the standard error of the log2(FoldChange). The new columns are appended at the end of the `results` table. `lfcSE` is `NA` when `stat.test = "wilcoxon"`
- updated manual and vignette accordingly to take into account for these exceptions and modifications

<br>

#### [v1.1.2](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/1.1.2) - June 16<sup>th</sup> 2026
- Update labeling of groups for `NES.plot` function
- Bug fixing for heatmap function for showing the dendrogram in the plots

<br>

#### [v1.1.1](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/1.1.1) - June 8<sup>th</sup> 2026
- Bug fixed in the PCA calculation, functions/objects concerned: `perform.PCA`, `diff.analyses*`, `plot.PC.buplot` and the `test.toolbox`, as well as the vignette and manual.

<br>

#### [v1.1.0](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/1.1.0) - May 31<sup>st</sup> 2026
- Added the function `plot.PC.biplot`.
- Added the function `SAINTq` and the `rime.dpo` and `rime.saintq` datasets.
- Updated `plot.PC.scatter` to allow for the separate plotting of x and y zero-lines.
- Updated the vignette to include the new functions.

<br>

#### [v1.0.0](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/1.0.0) - May 23<sup>rd</sup> 2026
- The result from `randomize.missing.values` is now included in a separate slot. Also the a new `boxplot.random` and `randomization.method` slots for the randomized scores have been added.
- In the DEprot.objects the slot `imputation` has been renamed into `imputation.method`. Many functions have been changed accordingly.
- Due to the addition of new slots, multiple functions have been adapted.
- In the differential result tables now degrees of freedom and test statistic columns have been added.
- For power calculation analyses, an estimation of the distribution of the statistics has been added in the results list. For this also the dependencies of the package have been updated.
- For GSEA, gene ranking can be based now also on the test statistic value. Accordingly, the `compare.ranking` function have been updated to compare all three methods.
- Update of the vignette to include the `compare.imp.methods` function and the new updates.
- Fixing a bug in handling single-sample groups for `randomize.missing.values` function.

<br>

#### [v0.1.1](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/0.1.1) - February 16<sup>th</sup> 2026
- Bug fixing on the `check.normality` function, which was inverting the evaluation of the AD's test p-value.
- Update of the vignette
- Update of the CITATION files

<br>

#### [v0.1.0](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/0.1.0) -  January 13<sup>th</sup> 2026
First release.



<br />
<br />

-----------------------------------------------------------------------

##### Contributors
[![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/DEprot?size=50&padding=5&bots=true)](https://sebastian-gregoricchio.github.io/)
