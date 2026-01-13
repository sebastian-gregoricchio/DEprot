.get_option <- function(name, default) {
  getOption(paste0("DEprot.", name), default)
}


#' @title Package options
#'
#' @description The \code{DEprot} package uses the following global options
#' (set via \code{\link[base:options]{options()}}):
#'
#' \describe{
#'   \item{\code{DEprot.update_check}}{
#'     Logical. Whether DEprot checks for package updates.
#'     Default is \code{TRUE}.
#'   }
#' }
#'
#' @name DEprot-options
#' @docType data
#' @keywords internal
NULL
