% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.MA.R
\name{plot.MA}
\alias{plot.MA}
\title{plot.MA}
\usage{
\method{plot}{MA}(
  DEprot.analyses.object,
  contrast = 1,
  up.color = "indianred",
  down.color = "steelblue",
  density.colors = colorRampPalette(colors = RColorBrewer::brewer.pal(9, "Blues"))(101),
  point.size = 2,
  point.alpha = 0.5,
  title = NULL,
  use.uncorrected.pvalue = FALSE,
  symmetric.y = TRUE,
  dot.labels = NULL,
  labels.in.boxes = FALSE,
  label.font.size = 2,
  label.max.overlaps = 100,
  min.segment.length.labels = 0
)
}
\arguments{
\item{DEprot.analyses.object}{An object of class \code{DEprot.analyses}.}

\item{contrast}{Number indicating the position of the contrast to use for the plotting.}

\item{up.color}{String indicating the color to use for up-regulated proteins in the plots. Default: \code{"indianred"}.}

\item{down.color}{String indicating the color to use for up-regulated proteins in the plots. Default: \code{"steelblue"}.}

\item{density.colors}{List of colors, passed to \code{scale_fill_gradientn}, to use for the density gradient. Default: \code{"colorRampPalette(colors = RColorBrewer::brewer.pal(9, "Blues"))(101)"}.}

\item{point.size}{Numeric value indicating the size of the dots. Default: \code{2}.}

\item{point.alpha}{Numeric value between 0 and 1 to indicate the transparency (alpha) of the dots. Default: \code{0.5}.}

\item{title}{String indicating the title to use. Default: \code{NULL} (automatic title).}

\item{use.uncorrected.pvalue}{Logical value indicating whether it should be used the normal p-value instead of the adjusted one (differential proteins numbers are recomputed). Default: \code{FALSE}, padj is used.}

\item{dot.labels}{String vector indicating labels to show on the plot that should correspond to \code{prot.id} column values. Default: \code{NULL} (no labels shown).}

\item{labels.in.boxes}{Logical value indicating whether the labels should be visualized as boxes. Default: \code{FALSE}.}

\item{label.font.size}{Numeric value indicating the size to use for the dot labels. Default: \code{2}.}

\item{label.max.overlaps}{Numeric value indicating the maximum number of overlaps allowed between labels. Default: \code{100}.}

\item{min.segment.length.labels}{Numeric value indicating the minimal length of the segments that connect the labels to the points. Default: \code{0} (segment always shown).}

\item{symmetric.x}{Logical values indicating whether the x-axis scale should be symmetric or not. Default: \code{TRUE}.}
}
\value{
A ggplot object.
}
\description{
Plots a MA plot log2(basemean) x log2(FoldChange) of differential expression results
}
