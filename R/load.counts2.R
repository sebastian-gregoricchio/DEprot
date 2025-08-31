#' @title load.counts2
#'
#' @description Function used to generate a \code{DEprot} object starting from counts and metadata.
#'
#' @param counts A data.frame or a matrix in which the rownames are the proteins and the columns the samples.
#' @param metadata A data.frame containing at least one column called \code{column.id} which corresponds to the colnames of \code{counts}. Any other column can be added and will correspond to a "feature"of each sample.
#' @param data.type String indicating the type of data that are loaded. One among: 'raw', 'normalized', 'imputed'.
#' @param log.base Number indicating the base of the log used to transform the counts. If none transformation is applied, indicate the base \code{1}.
#' @param normalization.method String or list indicating the normalization method used. If none, use the default value \code{NA}.
#' @param imputation.method A string indicating the imputation method used. If none, use the default value \code{NA}.
#' @param column.id String indicating the name of the column to use as "column.id" from the metadata data.frame. This column must contain all the colnames of \code{counts}.
#'
#' @return A \code{DEprot} object (S4 vector).
#'
#' @import dplyr
#' @import ggplot2
#' @import methods
#' @importFrom data.table fread
#' @importFrom reshape2 melt
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' dpo <- load.counts2(counts = DEprot::unimputed.counts,
#'                     metadata = DEprot::sample.config,
#'                     log.base = 2,
#'                     data.type = "raw")
#'
#'
#' @export load.counts2

load.counts2 =
  function(counts,
           metadata,
           data.type,
           log.base,
           normalization.method = NA,
           imputation.method = NA,
           column.id = "column.id") {

    # ### Libraries
    # require(dplyr)
    # require(ggplot2)
    # require(methods)


    ### Check counts
    check.rownames =
      function(x) {
        ### Check if "" (empty names) are present
        if ("" %in% rownames(x) | NA %in% rownames(x)) {
          stop("The rownames of the counts contain missing values ('') or NAs. Replace it with actual values.")
          #return(return())
          ### Check if rownames are duplicated
        } else if (TRUE %in% duplicated(rownames(counts))) {
          message("One or more rownames in the counts tables are duplicated. Only unique values are allowed. The `make.unique` function is applied on counts rownames unisng a '.' as separator.")
          rownames(x) = make.unique(names = rownames(x), sep = ".")
          return(x)
        } else {
          return(x)
        }
      }


    # Check counts type
    if (!(tolower(data.type) %in% c("raw", "r", "n", "nor", "norm", "normalized", "i", "imp", "im", "imputed"))) {
      stop("The `data.type` must be one among: 'raw', 'normalized', 'imputed'")
      #return()
    }


    ### Load intensities
    if ("data.frame" %in% class(counts) | "matrix" %in% class(counts)) {
      cnt = as.matrix(check.rownames(counts))
    } else{
      stop("The 'counts' table must be either a matrix or a data.frame. Rows are the protein.IDs and columns the samples.")
      #return()
    }


    # Remove rows with all NA values in the intensity matrix
    cnt[is.nan(cnt)] = NA

    pre.clean.nrow = nrow(cnt)
    cnt = cnt[rowSums(abs(cnt), na.rm = TRUE) > 0,]

    if (nrow(cnt) != pre.clean.nrow) {
      n.del.rows = pre.clean.nrow - nrow(cnt)
      message(paste("The counts matrix contained", n.del.rows, ifelse(n.del.rows == 1, yes = "row", no = "rows"), "with only NA values.\nThe latter", ifelse(n.del.rows == 1, yes = "has", no = "have"), "been removed from the matrix."))
    }



    ### Check metadata
    if ("data.frame" %in% class(metadata) | "matrix" %in% class(metadata)) {
      meta = as.data.frame(metadata)
    } else if ("character" %in% class(metadata)) {
      meta = data.table::fread(metadata, data.table = FALSE)
    } else {
      stop("The 'metadata' table must be either a matrix/data.frame or a string path. At least one column ('column.id') should contain the column names of the counts table.")
      #return()
    }

    # check that column.id is present in metadata
    if (!(column.id %in% colnames(meta))) {
      stop("The 'metadata' table must be either a matrix/data.frame or a string path. At least one column ('column.id') should contain the column names of the counts table.")
      #return()
    } else if (!all(sort(colnames(cnt)) == sort(meta[,column.id]))) {
      stop(paste0("Not all column names of the counts table correspond to the IDs indicated in the column '",
                  column.id,"' of the metadata table:\n",
                  "- counts IDs:\n", paste0(colnames(cnt), collapse = ", "), "\n\n",
                  "- metadata IDs ('",column.id,"'):\n", paste0(meta[,column.id], collapse = ", ")))
      #return()
    }


    # check whether the column.id is column.id, otherwise add it
    if (column.id != "column.id") {
      meta$column.id = meta[,column.id]
      meta = dplyr::relocate(.data = meta, column.id, .before = colnames(meta)[1])
    }




    ### Generate boxplot of counts
    # melt counts table
    melt.cnt =
      suppressMessages(reshape2::melt(as.data.frame(cnt))) %>%
      dplyr::mutate(variable = factor(variable, levels = colnames(cnt)))

    # compute stats
    cnt.stats =
      melt.cnt %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(min = min(value, na.rm = TRUE),
                       max = max(value, na.rm = TRUE))

    boxplot =
      ggplot() +
      geom_violin(data = melt.cnt,
                  mapping = aes(x = variable,
                                y = value,
                                group = variable),
                  width = 0.75,
                  alpha = 0.75,
                  fill = "darkorange",
                  color = NA) +
      geom_boxplot(data = melt.cnt,
                   mapping = aes(x = variable,
                                 y = value,
                                 group = variable),
                   fill = "white",
                   color = "darkorange3",
                   width = 0.15,
                   outlier.color = "black",
                   outlier.stroke = NA,
                   outlier.size = 2,
                   outlier.alpha = 0.25) +
      geom_line(data = data.frame(cnt.stats),
                mapping = aes(x = variable,
                              y = max,
                              group = 1),
                color = "indianred",
                linetype = 2,
                inherit.aes = FALSE) +
      geom_line(data = data.frame(cnt.stats),
                mapping = aes(x = variable,
                              y = min,
                              group = 1),
                color = "steelblue",
                linetype = 2,
                inherit.aes = FALSE) +
      ylab(ifelse(is.na(log.base),
                  yes = "Intensity",
                  no = paste0(ifelse(log.base == exp(1),
                                     yes = "ln", no = paste0("log<sub>",log.base,"</sub>")),
                              "(Intensity)"))) +
      ggtitle(ifelse(is.null(normalization.method),
                     yes = "**Unnormalized**",
                     no = paste0("**Normalized**<br>(",normalization.method,")"))) +
      xlab("Sample") +
      theme_classic() +
      theme(axis.text.y = element_text(color = "black"),
            axis.text.x = element_text(color = "black", hjust = 1, angle = 30),
            axis.title = ggtext::element_markdown(color = "black"),
            axis.ticks.y = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            plot.title = ggtext::element_markdown(hjust = 0.5),
            aspect.ratio = 10/ncol(cnt))


    ### Define DEprot object values depending on counts input type
    ## base values
    boxplot.raw = NA
    boxplot.norm = NA
    boxplot.imputed = NA
    raw.counts = NULL
    norm.counts = NULL
    imputed.counts = NULL
    imputed = FALSE
    normalized = FALSE


    # ------------------------------------------------------------------------------------
    ### Build final object
    if (tolower(data.type) %in% c("raw", "r")) {
      boxplot.raw = boxplot
      #boxplot.norm = NA
      #boxplot.imputed = NA
      raw.counts = cnt
      #norm.counts = NULL
      #imputed.counts = NULL
      #imputed = F
      #normalized = F
      #normalization.method = "none"
      #imputation.method = "none"

    } else if (tolower(data.type) %in% c("n", "nor", "norm", "normalized")) {
      #boxplot.raw = NA
      boxplot.norm = boxplot
      #boxplot.imputed = NA
      #raw.counts = NULL
      #imputed.counts = NULL
      #imputed = F
      norm.counts = cnt
      normalized = TRUE
      #imputation.method = "none"

    } else {
      #boxplot.raw = NA
      #boxplot.norm = NA
      boxplot.imputed = boxplot
      #raw.counts = NULL
      #norm.counts = NULL
      imputed.counts = cnt
      imputed = TRUE
      normalized = TRUE
    }

    # ------------------------------------------------------------------------------------



    ##################
    ### Building S4vector (DEprot object)
    DEprot.object =
      new(Class = "DEprot",
          metadata = meta,
          raw.counts = raw.counts,
          norm.counts = norm.counts,
          imputed.counts = imputed.counts,
          log.base = log.base,
          log.transformed = ifelse(is.null(log.base), yes = FALSE, no = TRUE),
          imputed = imputed,
          imputation = imputation.method,
          normalized = normalized,
          normalization.method = normalization.method,
          boxplot.raw = boxplot.raw,
          boxplot.norm = boxplot.norm,
          boxplot.imputed = boxplot.imputed,
          analyses.result.list = NULL,
          contrasts = NA,
          differential.analyses.params = NULL)


    return(DEprot.object)
  } # END function
