% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.and.methods.R
\docType{class}
\name{DEprot.analyses}
\alias{DEprot.analyses}
\title{DEprot.analyses class}
\description{
DEprot.analyses class
}
\section{Slots}{

 \describe{
   \item{\code{metadata}:}{The data.frame corresponding to the metadata table describing the samples. Class: \code{"ANY"}.}
   \item{\code{raw.counts}:}{Numeric matrix (rows: proteins, columns: samples) of the raw counts. Class: \code{"ANY"}.}
   \item{\code{norm.counts}:}{Numeric matrix (rows: proteins, columns: samples) of the normalized counts. Class: \code{"ANY"}.}
   \item{\code{imputed.counts}:}{Numeric matrix (rows: proteins, columns: samples) of the imputed counts. Class: \code{"ANY"}.}
   \item{\code{log.base}:}{Numeric value indicating the base of the logarithm expressing the values in the loaded data. Class: \code{"ANY"}.}
   \item{\code{log.transformed}:}{Logical value indicating whether the data are log-transformed or not. Class: \code{"logical"}.}
   \item{\code{imputed}:}{Logical value indicating whether the data are imputed. Class: \code{"logical"}.}
   \item{\code{imputation}:}{String (or any other class) value indicating the imputation method. Class: \code{"ANY"}.}
   \item{\code{normalized}:}{Logical value indicating whether the data are normalized. Class: \code{"logical"}.}
   \item{\code{normalization.method}:}{String (or any other class) value indicating the normalization method. Class: \code{"ANY"}. Class: \code{"ANY"}.}
   \item{\code{boxplot.raw}:}{Ggplot object showing the distribution of the raw values per sample. Class: \code{"ANY"}.}
   \item{\code{boxplot.norm}:}{Ggplot object showing the distribution of the normalized values per sample. Class: \code{"ANY"}.}
   \item{\code{boxplot.imputed}:}{Ggplot object showing the distribution of the imputed values per sample. Class: \code{"ANY"}.}
   \item{\code{analyses.result.list}:}{List containing the differential results for each contrast. Class: \code{"ANY"}. The list contains the following elements:
     \itemize{
       \item{\code{results}: }{a data.frame containing the results of the analyses; includes average expression of each group, basemean, foldchange, pvalue and p.adj, differential.status}
       \item{\code{n.diff}: }{a summary table showing the number of proteins in each differential expression status (up/down/unresponsive, null)}
       \item{\code{PCA.data}: }{output of \link{perform.PCA} for the subset of samples analyzed in a specific contrast}
       \item{\code{PCA.plots}: }{combination of 3 plots: scatter PC1-vs-PC2, scatter PC2-vs-PC3, and cumulative bar plot}
       \item{\code{correlations}: }{combination of Pearson and Spearman correlation heatmaps (obtained by \link{plot.correlation.heatmap}) for the subset of samples analyzed in a specific contrast}
       \item{\code{volcano}: }{volcano plot showing the log2(FoldChange) x -log10(p.adjusted) of differential expression results; it can be regenerated using \link{plot.volcano}}
       \item{\code{MA.plot}: }{MA-plot showing the log2(basemean) x log2(FoldChange) of differential expression results; it can be regenerated using \link{plot.MA}}}}
   \item{\code{contrasts}:}{List of contrasts. each contrast is a vector indicating, in the order: metadata.table.column - groupA - groupB; (groupA / group B). Class: \code{"ANY"}.}
   \item{\code{differential.analyses.params}:}{List of parameters used to run the differential analyses (fold change thresholds, p-value threshold, p-adjustement method, etc.). Class: \code{"ANY"}.}
 }
}

