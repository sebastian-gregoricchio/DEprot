![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/DEprot)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://sebastian-gregoricchio.github.io/DEprot/LICENSE.md/LICENSE)
[![R-CMD-check](https://github.com/sebastian-gregoricchio/DEprot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sebastian-gregoricchio/DEprot/actions/workflows/R-CMD-check.yaml)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/DEprot?style=social)](https://github.com/sebastian-gregoricchio/DEprot/fork)
<!-- ![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/DEprot)
![visits](https://badges.pufler.dev/visits/sebastian-gregoricchio/DEprot) -->
<!---![downloads](https://img.shields.io/github/downloads/sebastian-gregoricchio/DEprot/total.svg)--->

# DEprot [<img src="https://sebastian-gregoricchio.github.io/DEprot/DEprot_logo.png" align="right" height = 150/>](https://sebastian-gregoricchio.github.io/DEprot)

R package to impute, normalize and perform differential analyses on proteomics data (LFQ-based)


## Introduction
The concept behind `DEprot` (Differential Expression proteomics) is to provide a toolkit that allows for the normalization, imputation and analyses of the differential protein expression in proteomics data. The data are assumed to be LFQ (label-free quantitation) values.

[<img src="https://github.com/sebastian-gregoricchio/DEprot/blob/main/resources/DEprot_workflow.png" align="center" height=400 class="center"/>](https://sebastian-gregoricchio.github.io/DEprot)


### Citation
If you use this package, please cite:

<div class="warning" style='padding:2.5%; background-color:#ffffee; color:#787878; margin-left:5%; margin-right:5%; border-radius:15px;'>
<span>
<font size="-0.5">

<div style="margin-left:2%; margin-right:2%; text-align: justify">
*No publication associated yet.*
</div>
</font>

</span>
</div>

<br>


## Installation
```r
## Install devtools from CRAN (if not installed)
# install.packages("devtools")

## Or the development version from GitHub:
# install.packages("devtools")
# devtools::install_github("r-lib/devtools")

# Install the DEprot package
devtools::install_github("sebastian-gregoricchio/DEprot",
                         build_manual = TRUE,
                         build_vignettes = TRUE)
```

### Possibile installation issues
The package [`aPEAR`](https://github.com/kerseviciute/aPEAR) ([Kerseviciute & Gordevicius, Bioinformatics 2023](https://doi.org/10.1093/bioinformatics/btad672)) is required for the ORA/GSEA analyses. However, recently, it is has been removed from the CRAN, but it can be installed using:

```r
devtools::install_github("kerseviciute/aPEAR",
                         build_manual = FALSE,
                         build_vignettes = FALSE)
```


<br />

## Documentation
With the package a [web-manual](https://sebastian-gregoricchio.github.io/DEprot/reference/index.html) and a [vignette](https://sebastian-gregoricchio.github.io/DEprot/doc/DEprot.overview.vignette.html) are available.
The vignette can be inspected on R as well by typing `browseVignettes("DEprot")`.

Other examples of `DEprot` usage can be found on [Zenodo](https://doi.org/10.5281/zenodo.14823763).


<br />

## Package history and releases
A list of all releases and respective description of changes applied could be found [here](https://sebastian-gregoricchio.github.io/DEprot/NEWS).

<br />

-----------------
## Contact
For any suggestion, bug fixing, commentary please report it in the [issues](https://github.com/sebastian-gregoricchio/DEprot/issues)/[request](https://github.com/sebastian-gregoricchio/DEprot/pulls) tab of this repository.

## License
This package is under a [GNU General Public License (version 3)](https://sebastian-gregoricchio.github.io/DEprot/LICENSE.md/LICENSE).

<br />

#### Contributors
![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/DEprot?size=50&padding=5&bots=true)
