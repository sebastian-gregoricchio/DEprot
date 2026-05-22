---
title: "changeLog"
---

# DEprot [<img src="https://sebastian-gregoricchio.github.io/DEprot/DEprot_logo.png" align="right" height = 150/>](https://sebastian-gregoricchio.github.io/DEprot)
![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/DEprot)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/DEprot?style=social)](https://github.com/sebastian-gregoricchio/DEprot/fork)


#### [v1.0.0](https://github.com/sebastian-gregoricchio/DEprot/releases/tag/1.0.0) - May 22<sup>nd</sup> 2026
- The result from `randomize.missing.values` is now included in a separate slot. Also the boxplot is in a separate slot.
- In the DEprot.objects the slot `imputation` has been changed into `imputation.method`. Many functions have been changed accordingly.
- Due to the addition of new slots, multiple functions have been adapted.
- In the differential result tables now degrees of freedom and test statistic columns have been added.
- For power calculation analyses, an estimation of the distribution of the statistics has been added in the results list. For this also the dependencies of the package have been updated.
- In the GSEA gene ranking can be based now also on the test statistic value. Accordingly, the `compare.ranking` function have been updated to compare all three methods.
- Update of the vignette to include the `compare.imp.methods` function

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
