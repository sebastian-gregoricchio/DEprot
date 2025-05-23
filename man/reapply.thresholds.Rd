% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reapply.thresholds.R
\name{reapply.thresholds}
\alias{reapply.thresholds}
\title{reapply.thresholds}
\usage{
reapply.thresholds(
  DEprot.analyses.object,
  linear.FC = 2,
  p.adjusted = 0.05,
  linear.FC.unresp.range = c(1/1.1, 1.1),
  up.color = "indianred",
  down.color = "steelblue",
  unresponsive.color = "purple",
  null.color = "gray"
)
}
\arguments{
\item{DEprot.analyses.object}{An object of class \code{DEprot.analyses}.}

\item{linear.FC}{Number indicating the (absolute) fold change threshold (linear scale) to use to define differential proteins. Default: \code{2}.}

\item{p.adjusted}{Numeric value indicating the p.adjusted threshold to apply to the differential analyses. Default: \code{0.05}.}

\item{linear.FC.unresp.range}{A numeric 2-elements vector indicating the range (linear scale) used to define the unresponsive fold changes. Default: \code{c(1/1.1, 1.1)}.}

\item{up.color}{String indicating the color to use for up-regulated proteins in the plots. Default: \code{"indianred"}.}

\item{down.color}{String indicating the color to use for up-regulated proteins in the plots. Default: \code{"steelblue"}.}

\item{unresponsive.color}{String indicating the color to use for unresponsive proteins in the plots. Default: \code{"purple"}.}

\item{null.color}{String indicating the color to use for null proteins in the plots. Default: \code{"gray"}.}
}
\description{
Allows for the re-computation of the differential status and re-plotting volcano and MA plot for each contrast.
}
