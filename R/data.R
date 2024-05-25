#' @title Unimputed proteomics (LFQ) data
#' @description Dummy example of full proteomics data (LFQ values in log2). Not imputed.
#' @format A data frame with 13,239 rows and 12 columns:
#' \describe{
#'   \item{\code{rows}}{proteins}
#'   \item{\code{columns}}{samples}
#'}
#' @source Simulated data
"unimputed.counts"


#' @title Example of metadata table
#' @description Dummy table for metadata corresponding to the sample configuratiuon of \link{unimputed.counts}.
#' @format A data frame with 12 rows and 6 columns:
#' \describe{
#'   \item{\code{column.id}}{IDs of columns in the 'unimputed.counts'}
#'   \item{\code{sample.id}}{Actual IDs of the samples}
#'   \item{\code{cell}}{cell line ID}
#'   \item{\code{condition}}{Culture and treatment conditions}
#'   \item{\code{combined.id}}{ID combining cell and condition columns}
#'   \item{\code{replicate}}{biological replicate ID}
#'}
#' @source Simulated data
"sample.config"
