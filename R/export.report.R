#' @title export.report
#'
#' @description
#' Builds a self-contained HTML report from a \code{DEprot} or
#' \code{DEprot.analyses} object. The report always includes a quality-control
#' (QC) section (analysis parameters, value distributions, PCA of PC1/PC2/PC3
#' and sample correlation heatmaps) and, when differential analyses are present
#' in the object (\code{DEprot.analyses}), a results section reporting the
#' parameters used, a cross-contrast summary and, for each contrast, the volcano
#' plot, MA-plot and the top differential proteins. A final section lists the R
#' session and the versions of the packages used to generate the report.
#'
#' The function writes a parameterised R Markdown document to a temporary file
#' and renders it with \code{rmarkdown::render()}; all plots are produced with
#' the native DEprot plotting functions (\code{\link{perform.PCA}},
#' \code{\link{plot.PC.scatter.123}}, \code{\link{plot.PC.cumulative}},
#' \code{\link{plot.correlation.heatmap}}, \code{\link{plot.volcano}} and
#' \code{\link{plot.MA}}) or the plots already stored inside the object.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param output.file String with the path (or file name) of the HTML file to create. A \code{.html} extension is added when missing. Default: \code{"DEprot.report.html"}.
#' @param report.title String used as the report title. Default: \code{"DEprot report"}.
#' @param author.name String used as the report author. Default: \code{"DEprot"}.
#' @param which.data String indicating which counts to use for the QC section.
#'   One among 'raw', 'normalized'/'norm', 'randomized'/'random',
#'   'imputed'/'imp'. Default: \code{NULL}, in which case the best available data
#'   are chosen automatically (imputed > normalized > randomized > raw).
#' @param correlation.method Character vector with the correlation method(s) to display in the QC section. Any of 'pearson', 'spearman', 'kendall'. Default: \code{c("pearson", "spearman")}.
#' @param PCA.color.column String with the metadata column used to color the PCA points. Default: \code{"column.id"} (one color per sample).
#' @param PCA.shape.column String with the metadata column used for the PCA point shapes. Default: \code{NULL}.
#' @param PCA.label.column String with the metadata column used to label the PCA points. Default: \code{NULL}.
#' @param protein.summary.group.column String with the metadata column used to group samples in the protein-summary barplot (absolute protein counts). Default: \code{"column.id"} (one bar per sample).
#' @param volcano.use.uncorrected.pvalue Logical; if \code{TRUE} the volcano plots are regenerated using the uncorrected p-value. Default: \code{FALSE}.
#' @param show.MA.plot Logical; whether to include the MA-plot for each contrast. Default: \code{TRUE}.
#' @param include.contrast.qc Logical indicating whether to additionally include the per-contrast PCA and correlation plots stored in the object. Default: \code{TRUE}.
#' @param top.n.proteins Integer; number of top proteins (ranked by adjusted p-value) shown in the per-contrast results table. Default: \code{25}.
#' @param plot.width,plot.height Numeric; default figure width and height (in inches) used in the report. Default: \code{8} and \code{5.5}.
#' @param self.contained Logical indicating whether the HTML embeds all resources into a single portable file. Default: \code{TRUE}.
#' @param keep.Rmd Logical; if \code{TRUE} the intermediate \code{.Rmd} source is copied next to the output file. Default: \code{FALSE}.
#' @param quiet Logical indicating whether to suppress the rendering log. Default: \code{TRUE}.
#'
#' @return Invisibly, the path to the generated HTML file.
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' \dontrun{
#' # QC-only report (DEprot object)
#' export.report(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'               output.file = "QC.report.html",
#'               report.title = "Proteomics QC",
#'               author.name = "Jane Doe")
#'
#' # Full report with differential analyses (DEprot.analyses object)
#' export.report(DEprot.object = DEprot::test.toolbox$diff.exp.limma,
#'               output.file = "DE.report.html",
#'               report.title = "FBS vs DCC differential analysis",
#'               author.name = "Jane Doe",
#'               PCA.color.column = "condition",
#'               PCA.shape.column = "replicate")
#' }
#'
#'
#' @importFrom knitr kable opts_chunk
#' @importFrom rmarkdown render pandoc_available find_pandoc
#' @importFrom methods is slot
#' @importFrom stats setNames
#' @importFrom utils capture.output packageVersion sessionInfo str
#'
#' @export


export.report <- function(DEprot.object,
                          output.file = "DEprot.report.html",
                          report.title = "DEprot report",
                          author.name = "DEprot",
                          which.data = NULL,
                          correlation.method = c("pearson", "spearman"),
                          PCA.color.column = "column.id",
                          PCA.shape.column = NULL,
                          PCA.label.column = NULL,
                          protein.summary.group.column = "column.id",
                          volcano.use.uncorrected.pvalue = FALSE,
                          show.MA.plot = TRUE,
                          include.contrast.qc = TRUE,
                          top.n.proteins = 25,
                          plot.width = 8,
                          plot.height = 5.5,
                          self.contained = TRUE,
                          keep.Rmd = FALSE,
                          quiet = TRUE) {

  ## checks ------------------------------------------------------------------------------
  for (pkg in c("rmarkdown", "knitr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("The package '", pkg, "' is required to build the report. ",
           "Please install it with install.packages('", pkg, "').", call. = FALSE)
    }
  }
  if (!rmarkdown::pandoc_available()) {
    stop("Pandoc is required to render the HTML report but was not found. ",
         "Install RStudio or a standalone Pandoc (see rmarkdown::find_pandoc()).",
         call. = FALSE)
  }
  if (!methods::is(DEprot.object, "DEprot") &&
      !methods::is(DEprot.object, "DEprot.analyses")) {
    stop("'DEprot.object' must be an object of class 'DEprot' or 'DEprot.analyses'.",
         call. = FALSE)
  }

  ## are differential analyses available?
  is.analyses <- methods::is(DEprot.object, "DEprot.analyses") &&
    !is.null(DEprot.object@analyses.result.list) &&
    length(DEprot.object@analyses.result.list) > 0

  ## resolve data set --------------------------------------------------------------------
  get_mat <- function(slot.name) {
    m <- tryCatch(methods::slot(DEprot.object, slot.name), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    m <- suppressWarnings(try(as.matrix(m), silent = TRUE))
    if (inherits(m, "try-error") || is.null(m) || nrow(m) == 0 || ncol(m) == 0) return(NULL)
    m
  }
  available <- c(imputed    = !is.null(get_mat("imputed.counts")),
                 normalized = !is.null(get_mat("norm.counts")),
                 randomized = !is.null(get_mat("random.counts")),
                 raw        = !is.null(get_mat("raw.counts")))

  canon <- c(raw = "raw", normalized = "normalized", norm = "normalized",
             randomized = "randomized", random = "randomized",
             imputed = "imputed", imp = "imputed")
  if (!is.null(which.data)) {
    which.data <- match.arg(tolower(which.data), names(canon))
    which.data <- unname(canon[which.data])
    if (!isTRUE(available[[which.data]])) {
      warning("Requested 'which.data' = '", which.data,
              "' is not available in the object; using the best available data instead.",
              call. = FALSE)
      which.data <- NULL
    }
  }
  if (is.null(which.data)) {
    priority   <- c("imputed", "normalized", "randomized", "raw")
    which.data <- priority[which(available[priority])][1]
  }
  if (length(which.data) == 0 || is.na(which.data)) {
    stop("The object does not contain any usable counts matrix (raw/normalized/imputed).",
         call. = FALSE)
  }
  primary.slot <- c(imputed = "imputed.counts", normalized = "norm.counts",
                    randomized = "random.counts", raw = "raw.counts")[[which.data]]

  ## config ------------------------------------------------------------------------------
  rt <- if (length(report.title) >= 1 && nzchar(as.character(report.title)[1]))
    as.character(report.title)[1] else "DEprot report"
  an <- if (length(author.name) >= 1 && nzchar(as.character(author.name)[1]))
    as.character(author.name)[1] else "DEprot"

  cfg <- list(
    report.title   = rt,
    author.name    = an,
    generated.time = Sys.time(),
    which.data     = which.data,
    primary.slot   = primary.slot,
    correlation.method = tolower(as.character(correlation.method)),
    PCA.color.column   = PCA.color.column,
    PCA.shape.column   = PCA.shape.column,
    PCA.label.column   = PCA.label.column,
    protein.summary.group.column = protein.summary.group.column,
    volcano.use.uncorrected.pvalue = isTRUE(volcano.use.uncorrected.pvalue),
    show.MA.plot        = isTRUE(show.MA.plot),
    include.contrast.qc = isTRUE(include.contrast.qc),
    top.n.proteins      = as.integer(top.n.proteins),
    is.analyses         = isTRUE(is.analyses),
    plot.width          = plot.width,
    plot.height         = plot.height,
    logo.src            = "",
    self.contained      = isTRUE(self.contained),
    DEprot.version      = tryCatch(as.character(utils::packageVersion("DEprot")),
                                   error = function(e) NA_character_)
  )

  ## R Markdown template -----------------------------------------------------------------
  rmd.template <- r"---(---
title: "@@REPORT_TITLE@@"
author: "@@AUTHOR_NAME@@"
date: "@@GENERATED_DATE@@"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: true
    toc_depth: 3
    number_sections: true
    df_print: paged
    theme: flatly
    highlight: tango
    self_contained: @@SELF_CONTAINED@@
---

<style>
body { max-width: 1150px; }
.table-wrap { overflow-x: auto; margin-bottom: 1rem; }
table.table { font-size: 0.88em; }
.badge-row { margin: 0.4rem 0 1.4rem 0; }
.badge-box { display: inline-block; padding: 0.5rem 0.9rem; margin: 0.2rem 0.3rem 0.2rem 0;
  border-radius: 0.4rem; background: #eef3f7; border: 1px solid #d6e0ea; text-align: center; }
.badge-box .n { font-size: 1.35em; font-weight: 700; color: #2c3e50; display: block; }
.badge-box .l { font-size: 0.72em; color: #667; text-transform: uppercase; letter-spacing: 0.03em; }
img { max-width: 100%; height: auto; }
.note { color: #778; font-size: 0.9em; }
h1 { margin-top: 1.6rem; }
h1.title { color: #000080; font-weight: 700; }
.citation-box { background:#f4f7fa; border-left:4px solid #000080; padding:0.8rem 1rem;
  margin:0.6rem 0 1.2rem 0; border-radius:0.3rem; font-size:0.92em; line-height:1.5; }
.citation-box .cite-label { font-weight:700; color:#000080; display:block; margin-bottom:0.4rem; }
.side-panel { position: fixed; top: 70px; right: 18px; width: 220px; z-index: 20; }
.side-panel .side-logo { display:block; width: 82%; max-width: 165px; height:auto;
  margin: 0 auto 0.7rem; }
.side-panel .citation-box { font-size: 0.78em; margin: 0; }
/* On viewports too narrow for a right margin, drop the panel into the normal
   flow at the top of the report instead of overlapping the content. */
@media (max-width: 1000px) {
  .side-panel { position: static; width: auto; margin: 0.2rem 0 1.2rem 0; }
  .side-panel .side-logo { margin: 0 0 0.6rem 0; max-width: 150px; }
  .side-panel .citation-box { font-size: 0.92em; }
}
</style>

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,
                      fig.align = "center", dev = "png", dpi = 150,
                      fig.width = .cfg$plot.width, fig.height = .cfg$plot.height)
options(knitr.kable.NA = "")
suppressPackageStartupMessages(library(DEprot))

dpo <- get(".dpo")
cfg <- get(".cfg")

have <- function(p) requireNamespace(p, quietly = TRUE)
if (have("ggplot2")) suppressPackageStartupMessages(library(ggplot2))
if (have("legendry")) suppressWarnings(try(suppressPackageStartupMessages(require(legendry)), silent = TRUE))

## --- helpers -------------------------------------------------------------
esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;",  s, fixed = TRUE)
  gsub(">", "&gt;",  s, fixed = TRUE)
}

fmt_val <- function(v) {
  if (is.null(v)) return("&mdash;")
  if (is.logical(v) && length(v) == 1) return(if (isTRUE(v)) "yes" else "no")
  if (is.atomic(v)) return(paste(format(v, trim = TRUE), collapse = ", "))
  paste(trimws(utils::capture.output(utils::str(v, give.attr = FALSE, no.list = TRUE))),
        collapse = " | ")
}

show_table <- function(df, caption = NULL, max_rows = NULL, striped = TRUE) {
  if (is.null(df) || (is.data.frame(df) && nrow(df) == 0)) {
    cat("<p class='note'>No data available.</p>\n"); return(invisible())
  }
  df <- as.data.frame(df, check.names = FALSE)
  note <- NULL
  if (!is.null(max_rows) && nrow(df) > max_rows) {
    note <- sprintf("<p class='note'>Showing top %d of %s rows.</p>",
                    max_rows, format(nrow(df), big.mark = ","))
    df <- df[seq_len(max_rows), , drop = FALSE]
  }
  cls <- paste0("table table-sm", if (striped) " table-striped" else "")
  cat("<div class='table-wrap'>\n")
  cat(knitr::kable(df, format = "html", escape = TRUE, caption = caption,
                   row.names = FALSE,
                   table.attr = sprintf('class="%s"', cls)), sep = "\n")
  cat("\n</div>\n")
  if (!is.null(note)) cat(note, "\n")
}

kv_table <- function(x, caption = NULL) {
  if (is.null(x) || length(x) == 0) {
    cat("<p class='note'>No parameters recorded.</p>\n"); return(invisible())
  }
  nm <- names(x)
  if (is.null(nm)) nm <- rep("", length(x))
  empty <- is.na(nm) | nm == ""
  if (any(empty)) nm[empty] <- paste0("param_", which(empty))
  df <- data.frame(Parameter = nm,
                   Value = vapply(unname(x), fmt_val, character(1)),
                   check.names = FALSE, stringsAsFactors = FALSE)
  show_table(df, caption = caption)
}

round_df <- function(df, digits = 3) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(z) signif(z, digits))
  df
}

is_plot <- function(x) inherits(x, c("gg", "ggplot", "patchwork"))

safe_print <- function(expr, what = "plot") {
  p <- tryCatch(expr, error = function(e) e)
  if (inherits(p, "error")) {
    cat(sprintf("<p class='note'>Could not generate the %s: %s</p>\n",
                what, esc(conditionMessage(p))))
    return(invisible())
  }
  if (methods::is(p, "DEprot.correlation")) { print(p@heatmap); return(invisible()) }
  print(p)
}

## The @contrasts slot stores, per contrast, either a named list
## (metadata.column / var.1 / var.2) or a 3-element c(column, groupA, groupB)
## vector. This normalises both to $column / $A / $B (all character scalars).
contrast_parts <- function(i) {
  ct <- tryCatch(dpo@contrasts[[i]], error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  col <- A <- B <- NULL
  if (is.list(ct)) {
    A <- ct[["var.1"]]; B <- ct[["var.2"]]; col <- ct[["metadata.column"]]
    if (is.null(A) || is.null(B)) {           # unnamed list fallback
      if (length(ct) >= 3) { col <- ct[[1]]; A <- ct[[2]]; B <- ct[[3]] }
    }
  } else if (length(ct) >= 3) {
    col <- ct[1]; A <- ct[2]; B <- ct[3]
  }
  if (is.null(A) || is.null(B)) return(NULL)
  list(column = if (is.null(col)) NA_character_ else as.character(col)[1],
       A = as.character(A)[1], B = as.character(B)[1])
}

contrast_label <- function(i) {
  nms <- names(dpo@analyses.result.list)
  if (!is.null(nms) && length(nms) >= i && !is.na(nms[i]) && nzchar(nms[i])) return(nms[i])
  cp <- contrast_parts(i)
  if (!is.null(cp)) return(sprintf("%s vs %s (%s)", cp$A, cp$B, cp$column))
  paste("Contrast", i)
}
```

```{r side-panel, results='asis'}
## Fixed right-side panel (opposite the floating TOC) holding the DEprot logo and
## the citation. Emitted as one raw-HTML block so the DOI/PMID render as clean links.
logo.html <- if (nzchar(cfg$logo.src)) {
  sprintf("<img src='%s' alt='DEprot logo' class='side-logo'>", cfg$logo.src)
} else ""
cat(
  "<div class='side-panel'>",
  logo.html,
  "<div class='citation-box'>",
  "<div class='cite-label'>If you use DEprot, please cite:</div>",
  "Eickhoff N, Hoekman L, Bleijerveld O, Bergman AM, Zwart W, Gregoricchio S (2026). ",
  "&ldquo;DEprot: a comprehensive R-package for the analyses of label-free quantitation ",
  "mass-spectrometry data.&rdquo; <em>NAR Genomics and Bioinformatics</em>, 8(1), 1&ndash;12.<br>",
  "PMID: <a href='https://pubmed.ncbi.nlm.nih.gov/41685351/'>41685351</a>",
  " &nbsp;&middot;&nbsp; ",
  "doi: <a href='https://doi.org/10.1093/nargab/lqag015'>10.1093/nargab/lqag015</a>",
  "</div>",
  "</div>",
  sep = ""
)
```

<br>

<hr style="border:2px solid #C9C9C9">

# **Overview**

```{r overview, results='asis'}
mat     <- tryCatch(as.matrix(methods::slot(dpo, cfg$primary.slot)), error = function(e) NULL)
n.prot  <- if (!is.null(mat)) nrow(mat) else NA
n.samp  <- if (!is.null(mat)) ncol(mat) else NA
n.contr <- if (isTRUE(cfg$is.analyses)) length(dpo@analyses.result.list) else 0L

badge <- function(n, l) sprintf(
  "<span class='badge-box'><span class='n'>%s</span><span class='l'>%s</span></span>", n, l)

cat("<div class='badge-row'>")
cat(badge(if (is.na(n.prot)) "&mdash;" else format(n.prot, big.mark = ","), "proteins"))
cat(badge(if (is.na(n.samp)) "&mdash;" else n.samp, "samples"))
if (isTRUE(cfg$is.analyses)) cat(badge(n.contr, "contrasts"))
cat(badge(cfg$which.data, "data used"))
cat("</div>\n\n")

cat("This report was generated with **DEprot**",
    if (is.na(cfg$DEprot.version)) "" else paste0(" v", cfg$DEprot.version),
    ". It summarises the quality-control metrics of the dataset",
    if (isTRUE(cfg$is.analyses))
      " together with the differential-expression analyses stored in the object." else ".",
    "\n\n", sep = "")
```

<br>

## Analysis parameters

```{r params-qc, results='asis'}
## Blocks are shown in processing order: raw -> normalization -> randomization -> imputation.
## For the normalization/imputation methods only the parameters actually passed to the
## corresponding DEprot function are kept (diagnostic/timing outputs are dropped).
diagnostic.fields <- c("OOBerror", "variable.wise.OOBerror", "processing.time")

## --- Raw data ---
cat("**Raw data**\n\n")
kv_table(list("Log-transformed"    = dpo@log.transformed,
              "Log base"           = dpo@log.base,
              "Counts used for QC" = cfg$which.data))

## --- Normalization ---
cat("\n**Normalization**\n\n")
if (isTRUE(dpo@normalized) && !is.null(dpo@normalization.method)) {
  nm <- dpo@normalization.method
  if (is.data.frame(nm) && ncol(nm) >= 2) {
    df <- data.frame(Parameter = as.character(nm[[1]]),
                     Value     = vapply(nm[[2]], fmt_val, character(1), USE.NAMES = FALSE),
                     check.names = FALSE, stringsAsFactors = FALSE)
    show_table(df)
  } else {
    kv_table(nm)
  }
} else {
  cat("<p class='note'>Data were not normalized within DEprot.</p>\n")
}

## --- Randomization ---
cat("\n**Randomization**\n\n")
if (isTRUE(dpo@randomized) && !is.null(dpo@randomization.method) && length(dpo@randomization.method) > 0) {
  kv_table(dpo@randomization.method)
} else {
  cat("<p class='note'>No randomization of missing values was applied.</p>\n")
}

## --- Imputation ---
cat("\n**Imputation**\n\n")
if (isTRUE(dpo@imputed) && !is.null(dpo@imputation.method) && length(dpo@imputation.method) > 0) {
  im <- dpo@imputation.method
  if (is.list(im)) im <- im[setdiff(names(im), diagnostic.fields)]
  kv_table(im)
} else {
  cat("<p class='note'>Data were not imputed.</p>\n")
}
```

## Sample metadata

```{r metadata, results='asis'}
md <- tryCatch(as.data.frame(dpo@metadata, check.names = FALSE), error = function(e) NULL)
show_table(md)
```

## Value distributions {.tabset .tabset-pills}

```{r distributions, results='asis', fig.width=9, fig.height=5}
box.slots <- list("Raw" = "boxplot.raw", "Normalized" = "boxplot.norm",
                  "Randomized" = "boxplot.random", "Imputed" = "boxplot.imputed")
any.box <- FALSE
for (lab in names(box.slots)) {
  bp <- tryCatch(methods::slot(dpo, box.slots[[lab]]), error = function(e) NULL)
  if (!is.null(bp) && is_plot(bp)) {
    any.box <- TRUE
    cat(sprintf("\n\n### %s\n\n", lab))
    safe_print(bp, what = paste(tolower(lab), "distribution plot"))
    cat("\n\n")
  }
}
if (!any.box) cat("<p class='note'>No distribution plots are stored in this object.</p>\n")
```

## Protein counts per sample

```{r protein-summary, results='asis', fig.width=9, fig.height=5}
## Grouped by the chosen metadata column and shown as absolute counts (show.frequency = FALSE).
ps <- tryCatch(
  protein.summary(DEprot.object = dpo,
                  group.column   = cfg$protein.summary.group.column,
                  show.frequency = FALSE),
  error = function(e) e)
if (inherits(ps, "error")) {
  cat(sprintf("<p class='note'>The protein summary could not be generated: %s</p>\n",
              esc(conditionMessage(ps))))
} else {
  safe_print(ps, what = "protein-summary barplot")
}
```

<br>

--------------------------------------------------------

# **Quality control**
## Principal component analysis (PC1&ndash;PC2&ndash;PC3)

```{r pca-compute}
pca.obj <- tryCatch(
  perform.PCA(DEprot.object = dpo, which.data = cfg$which.data),
  error = function(e) e)
```

```{r pca-scatter, results='asis', fig.width=9.45, fig.height=3.7}
if (inherits(pca.obj, "error")) {
  cat(sprintf("<p class='note'>PCA could not be computed: %s</p>\n", esc(conditionMessage(pca.obj))))
} else {
  safe_print(
    plot.PC.scatter.123(DEprot.PCA.object = pca.obj,
                        color.column = cfg$PCA.color.column,
                        shape.column = cfg$PCA.shape.column,
                        label.column = cfg$PCA.label.column),
    what = "PCA scatter (PC1-PC2 / PC2-PC3)")
  cat("\n\n")
}
```

```{r pca-cumulative, results='asis', fig.width=8, fig.height=4.5}
if (!inherits(pca.obj, "error")) {
  safe_print(plot.PC.cumulative(pca.obj), what = "PCA cumulative-variance plot")
}
```

## Sample correlation {.tabset .tabset-pills}

```{r correlation, results='asis', fig.width=7.5, fig.height=6.5}
methods.vec <- cfg$correlation.method
if (is.null(methods.vec) || length(methods.vec) == 0) methods.vec <- "pearson"
for (method in methods.vec) {
  cat(sprintf("\n\n### %s\n\n", tools::toTitleCase(method)))
  corr <- tryCatch(
    plot.correlation.heatmap(DEprot.object = dpo,
                             correlation.method = method,
                             which.data = cfg$which.data,
                             correlation.scale.limits = c(NA, 1)),
    error = function(e) e)
  if (inherits(corr, "error")) {
    cat(sprintf("<p class='note'>Could not compute the %s correlation: %s</p>\n",
                method, esc(conditionMessage(corr))))
  } else {
    safe_print(corr, what = paste(method, "correlation heatmap"))
  }
  cat("\n\n")
}
```
<br>

--------------------------------------------------------

```{r results-section, eval=cfg$is.analyses, results='asis', fig.width=8.5, fig.height=6.2}
cat("# **Differential expression results**\n\n")

## ---- parameters used for the differential analyses ----
cat("## Differential analysis parameters\n\n")
kv_table(dpo@differential.analyses.params)

## ---- cross-contrast summary ----
cat("\n## Summary of differential analyses\n\n")
res.list <- dpo@analyses.result.list

## The status column is 'diff.status'; its levels are the two contrast group names
## (proteins enriched in each group) plus 'unresponsive' and 'null'. Prefer the
## precomputed per-contrast counts stored in '$n.diff' when present.
count_status <- function(rl) {
  nd <- rl$n.diff
  if (!is.null(nd) && is.data.frame(nd) && all(c("diff.status", "n") %in% names(nd))) {
    return(stats::setNames(as.integer(nd$n), as.character(nd$diff.status)))
  }
  r <- rl$results
  if (is.null(r) || is.null(r$diff.status)) return(NULL)
  tb <- table(as.character(r$diff.status))
  stats::setNames(as.integer(tb), names(tb))
}

## Map each contrast's status counts onto generic columns so the table stays at
## four columns regardless of how many contrasts (or group names) are present:
## group.A / group.B are the first and second groups of the contrast (see 'Contrast').
pick <- function(z, key) {
  if (is.null(z) || length(key) == 0) return(0L)
  key <- as.character(key)[1]
  if (is.na(key) || !key %in% names(z)) return(0L)
  v <- z[[key]]; if (is.null(v) || is.na(v)) 0L else as.integer(v)
}
summ.row <- function(i, rl) {
  z  <- count_status(rl)
  cp <- contrast_parts(i)
  if (!is.null(cp)) {
    a <- pick(z, cp$A); b <- pick(z, cp$B)
  } else {                     # fallback: the two non-tail statuses, in stored order
    grp <- setdiff(names(z), c("unresponsive", "null"))
    a <- if (length(grp) >= 1) pick(z, grp[1]) else 0L
    b <- if (length(grp) >= 2) pick(z, grp[2]) else 0L
  }
  data.frame("group.A" = a, "group.B" = b,
             "unresponsive" = pick(z, "unresponsive"),
             "null" = pick(z, "null"),
             check.names = FALSE)
}

overview <- data.frame(Contrast = vapply(seq_along(res.list), contrast_label, character(1)),
                       check.names = FALSE, stringsAsFactors = FALSE)
counts.df <- do.call(rbind, lapply(seq_along(res.list), function(i) summ.row(i, res.list[[i]])))
overview  <- cbind(overview, counts.df)
overview[["Total tested"]] <- rowSums(counts.df)
show_table(overview, caption = paste(
  "Number of proteins per differential-expression status and contrast.",
  "'group.A' and 'group.B' are, respectively, the first and second group of each",
  "contrast (as shown in the 'Contrast' column) and hold the proteins significantly",
  "enriched in that group; 'unresponsive' and 'null' are non-significant proteins."))

## ---- per-contrast detail (one tab per contrast) ----
cat("\n## Results per contrast {.tabset .tabset-pills}\n\n")
for (i in seq_along(res.list)) {
  rl <- res.list[[i]]
  cat(sprintf("\n\n### %s\n\n", contrast_label(i)))

  ## Volcano plot
  cat("\n#### Volcano plot\n\n")
  volc <- rl$volcano
  if (!isTRUE(cfg$volcano.use.uncorrected.pvalue) && !is.null(volc) && is_plot(volc)) {
    safe_print(volc, what = "volcano plot")
  } else {
    safe_print(plot.volcano(DEprot.analyses.object = dpo, contrast = i,
                            use.uncorrected.pvalue = isTRUE(cfg$volcano.use.uncorrected.pvalue)),
               what = "volcano plot")
  }
  cat("\n\n")

  ## MA-plot
  if (isTRUE(cfg$show.MA.plot)) {
    cat("\n#### MA-plot\n\n")
    ma <- rl$MA.plot
    if (!is.null(ma) && is_plot(ma)) safe_print(ma, what = "MA-plot")
    else safe_print(plot.MA(DEprot.analyses.object = dpo, contrast = i), what = "MA-plot")
    cat("\n\n")
  }

  ## Top differential proteins
  cat("\n#### Top proteins\n\n")
  r <- rl$results
  if (!is.null(r) && nrow(r) > 0) {
    r <- as.data.frame(r, check.names = FALSE)
    ord.col <- intersect(c("padj", "p.adj", "p_adj", "padjusted", "adj.pvalue",
                           "pvalue", "p.value", "p_value"), names(r))
    if (length(ord.col) > 0) r <- r[order(r[[ord.col[1]]]), , drop = FALSE]
    show_table(round_df(r, 3), max_rows = cfg$top.n.proteins)
    cat("<p class='note'>The complete results table can be exported with ",
        "<code>DEprot::export.analyses()</code> or <code>DEprot::get.results()</code>.</p>\n")
  } else {
    cat("<p class='note'>No results table available for this contrast.</p>\n")
  }

  ## Optional per-contrast QC (stored in the object)
  if (isTRUE(cfg$include.contrast.qc)) {
    if (!is.null(rl$PCA.plots) && is_plot(rl$PCA.plots)) {
      cat("\n#### PCA (contrast subset)\n\n")
      safe_print(rl$PCA.plots, what = "contrast PCA"); cat("\n\n")
    }
    if (!is.null(rl$correlations) && is_plot(rl$correlations)) {
      cat("\n#### Correlation (contrast subset)\n\n")
      safe_print(rl$correlations, what = "contrast correlation"); cat("\n\n")
    }
  }
}
```

<br>

--------------------------------------------------------

# **Packages and session information**

```{r session, results='asis'}
si <- utils::sessionInfo()
cat(sprintf(
  "<p><strong>R version:</strong> %s &nbsp;|&nbsp; <strong>Platform:</strong> %s &nbsp;|&nbsp; <strong>DEprot:</strong> %s</p>\n",
  esc(si$R.version$version.string), esc(si$platform),
  if (is.na(cfg$DEprot.version)) "NA" else cfg$DEprot.version))

pkg_df <- function(lst, status) {
  if (is.null(lst) || length(lst) == 0) return(NULL)
  data.frame(Package = names(lst),
             Version = vapply(lst, function(p) as.character(p$Version), character(1)),
             Loaded  = status, check.names = FALSE, stringsAsFactors = FALSE)
}
tab <- rbind(pkg_df(si$otherPkgs, "attached"), pkg_df(si$loadedOnly, "namespace"))
if (!is.null(tab) && nrow(tab) > 0) {
  tab <- tab[order(tolower(tab$Package)), , drop = FALSE]
  rownames(tab) <- NULL
  show_table(tab, caption = "Packages (with versions) used to build this report.")
}

cat("\n<details>\n<summary>Full session information</summary>\n<pre>\n")
cat(esc(paste(utils::capture.output(print(si)), collapse = "\n")))
cat("\n</pre>\n</details>\n")
```

---
<p class='note'>Report generated by <code>DEprot::export.report()</code> on `r format(.cfg$generated.time, '%Y-%m-%d %H:%M %Z')`.</p>
)---"

## write & render ------------------------------------------------------------------------
rmd.path <- tempfile(pattern = "DEprot.report.", fileext = ".Rmd")
## Bake YAML front-matter values in directly (avoids fragile inline R in
## YAML; an unquoted `r ...` there is invalid YAML). Escape for a
## double-quoted YAML scalar: backslash and double-quote.
yaml_dq <- function(x) {
  x <- gsub("\\", "\\\\", as.character(x)[1], fixed = TRUE)
  gsub("\"", "\\\"", x, fixed = TRUE)
}

ord <- function(d) {
  d <- as.integer(d)
  if (d %% 100 %in% 11:13) "th" else c("th","st","nd","rd","th","th","th","th","th","th")[d %% 10 + 1]
}
gt <- cfg$generated.time
date.str <- paste0(format(gt, "%d"), "<sup>", ord(format(gt, "%d")),
                   "</sup> ", format(gt, "%B %Y, %H:%M %Z"))

rmd.template <- gsub("@@REPORT_TITLE@@",  yaml_dq(cfg$report.title), rmd.template, fixed = TRUE)
rmd.template <- gsub("@@AUTHOR_NAME@@",   yaml_dq(cfg$author.name),  rmd.template, fixed = TRUE)
rmd.template <- gsub("@@GENERATED_DATE@@", yaml_dq(date.str), rmd.template, fixed = TRUE)
rmd.template <- gsub("@@SELF_CONTAINED@@",
                     if (isTRUE(self.contained)) "true" else "false",
                     rmd.template, fixed = TRUE)

cat(rmd.template, file = rmd.path)

output.file <- path.expand(output.file)
if (!grepl("\\.html?$", output.file, ignore.case = TRUE)) {
  output.file <- paste0(output.file, ".html")
}
out.dir  <- dirname(output.file)
out.base <- basename(output.file)
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

## Download the DEprot logo into the output folder as a temporary file, embed it
## as a data URI (so the report stays self-contained), then remove the temp file.
logo.url <- "https://raw.githubusercontent.com/sebastian-gregoricchio/DEprot/main/DEprot_logo.png"
logo.tmp <- tempfile(pattern = "DEprot_logo_", tmpdir = out.dir, fileext = ".png")
cfg$logo.src <- tryCatch({
  suppressWarnings(utils::download.file(logo.url, destfile = logo.tmp,
                                        mode = "wb", quiet = TRUE))
  if (file.exists(logo.tmp) && file.info(logo.tmp)$size > 0) {
    knitr::image_uri(logo.tmp)
  } else ""
}, error = function(e) "")
if (file.exists(logo.tmp)) unlink(logo.tmp)

report.env <- new.env(parent = globalenv())
assign(".dpo", DEprot.object, envir = report.env)
assign(".cfg", cfg,           envir = report.env)

out.path <-
  rmarkdown::render(input = rmd.path,
                    output_file = out.base,
                    output_dir = out.dir,
                    intermediates_dir = tempdir(),
                    envir = report.env,
                    quiet = quiet)

if (isTRUE(keep.Rmd)) {
  kept <- file.path(out.dir, sub("\\.html?$", ".Rmd", out.base, ignore.case = TRUE))
  file.copy(rmd.path, kept, overwrite = TRUE)
  if (!quiet) message("R Markdown source kept at: ", kept)
}
if (!quiet) message("Report written to: ", normalizePath(out.path))

invisible(out.path)
}
