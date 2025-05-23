% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/protein.summary.R
\name{protein.summary}
\alias{protein.summary}
\title{protein.summary}
\usage{
protein.summary(
  DEprot.object,
  group.column = "column.id",
  sample.subset = NULL,
  n.labels = NULL,
  label.color = "white",
  x.label.angle = 30,
  show.frequency = FALSE,
  colors = NULL,
  title = NULL,
  subtitle = NULL
)
}
\arguments{
\item{DEprot.object}{An object of class \code{DEprot} or \code{DEprot.analyses}.}

\item{group.column}{String indicating a column among the ones in the metadata table. Default: \code{"column.id"} (no groups).}

\item{sample.subset}{String vector indicating the column names (samples) to keep in the counts table (the 'column.id' in the metadata table). If subset is applied, the proteins with all NAs will be removed (total number of protains will be alterated). Default: \code{NULL} (no subsetting).}

\item{n.labels}{String indicating the type of values to display on the barplot. One among: \code{NULL} (no labels), "frequency", "percentage", "counts". Default: \code{NULL} (no labels).}

\item{label.color}{String indicating the font-color to use for the barplot values labels. Default: \code{"white"}.}

\item{x.label.angle}{Numeric value indicating the rotation angle to use for the x-axis labels. Default: \code{30} (degrees).}

\item{show.frequency}{Logical value indicating whether the barplot should show the y-axis as frequency or as absolute counts. Default: \code{FALSE} (counts).}

\item{colors}{String-vector indicating the colors to use for the degrees of protein-presence. Default: \code{NULL} (colors will be assigned automatically).}

\item{title}{String indicating the title to use (markdown annotation supported). Default: \code{NULL}.}

\item{subtitle}{String indicating the subtitle to use (markdown annotation supported). Default: \code{NULL}.}
}
\value{
A ggplot object.
}
\description{
Plots stacked-barplot recapitulating
}
