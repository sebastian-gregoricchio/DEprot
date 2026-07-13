---
title: "changeLog"
---

# DEprot [<img src="https://sebastian-gregoricchio.github.io/DEprot/DEprot_logo.png" align="right" height = 150/>](https://sebastian-gregoricchio.github.io/DEprot)
![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/DEprot)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/DEprot?style=social)](https://github.com/sebastian-gregoricchio/DEprot/fork)


#### [v1.2.1](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/1.2.1) - July 13<sup>th</sup> 2026
- `plot.volcano` and `plot.MA` can automatically plot the top N differential proteins
- `export.report` will plot the top.n proteins in the volcano
- updated vignette and manual accordingly


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
