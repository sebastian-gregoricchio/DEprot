#' @title Sample reorderer computeMatrix file
#'
#' @description A tool to reorder the samples in a computeMatrix file avoiding the re-computation of the latter.
#'
#'
#' @param matrix.file String with the full path to a deeptools computeMatrix \code{.gz} file.
#' @param new.sample.order String vector with the sample labels in the order in which they should appear in the matrix. By default {NULL}, which returns a message with the sample labels in the original order. Further, it returns the vector with the original sample order.
#' @param reordered.matrix.path String with full path to for the file of the reorder matrix (.gz). By default the output name will be \code{<original.matrix.name>_reordered.gz}.
#' @param ignore.header.error Logical value to indicate whether the error of sample_label reassignment in the header should be ignored. The plotted labels can be changed during the plotting.
#' @param verbose Logical value to indicate whether the final message should be printed. By default \code{TRUE}.
#'
#' @return The output is a computeMatrix file (.gz format) with the samples chucks re-shuffled to be in the order provided by the user.
#' 
#' @export reorder.samples.computeMatrix



reorder.samples.computeMatrix =
  function(matrix.file,
           new.sample.order = NULL,
           reordered.matrix.path = gsub(".gz", "_reordered.gz", matrix.file),
           ignore.header.error = FALSE,
           verbose = TRUE) {
    
    ######################################################################################
    ### Required libraries
    # require(data.table)
    # require(stringr)
    require(dplyr)
    
    # Check if Rseb is up-to-date #
    Rseb::actualize(update = F, verbose = F)
    ######################################################################################
    
    
    ######## Check whether the matrix file exits and whether I can write in the output folder
    
    if (!file.exists(matrix.file)) {
      return(warning("The matrix file provided does not exists."))
    } else if (file.access(dirname(reordered.matrix.path), 2)[[1]] == -1) {
      return(warning("The user does not own writing permissions for the output folder. The output matrix can't be exported in the location provided."))
    }
    
  
    
    ######## Read the matrix and get the sample/data_boundaries info
    matrix = Rseb::read.computeMatrix.file(matrix.file)
    
    sample_start_column = as.numeric(as.vector(stringr::str_split(string = dplyr::filter(.data = matrix$metadata, parameters == "sample_boundaries")$values, pattern = ",")[[1]]))
    sample_names = as.vector(stringr::str_split(string = dplyr::filter(.data = matrix$metadata, parameters == "sample_labels")$values, pattern = ",")[[1]])
    
    # group_start_row = as.numeric(as.vector(stringr::str_split(string = dplyr::filter(.data = matrix$metadata, parameters == "group_boundaries")$values, pattern = ",")[[1]]))
    # group_names = as.vector(stringr::str_split(string = dplyr::filter(.data = matrix$metadata, parameters == "group_labels")$values, pattern = ",")[[1]])
    
    number.info.columns = ncol(matrix$matrix.data) - sample_start_column[length(sample_start_column)]
    
    
    ######## Check that the new names are in the previous one
    if (!is.null(new.sample.order)) {
      if (unique(!(new.sample.order %in% sample_names)) > 1 |
          unique(!(sample_names %in% new.sample.order)) > 1 |
          sort(unique(sample_names %in% new.sample.order))[1] == F |
          sort(unique(new.sample.order %in% sample_names))[1] == F) {
        
        original_order = paste0(paste(paste0("   ",1:length(sample_names), "."), sample_names), collapse = "\n")
        warning("The labels in the new.sample.order are not corresponding to the sample labels of the current matrix:\n",
                original_order)
        return(sample_names)
      }
    } else {
      warning("Now new.sample.order provided. The current order is:\n",
              original_order)
      return(sample_names)
    }
    
    
    
    ######## Generate a list with a table per sample
    samples.table.list = list()
    
    for (i in 1:length(sample_names)) {
      # define the limits of the sample keeping the first lines
      start.col = sample_start_column[i] + number.info.columns + 1
      end.col = sample_start_column[i+1] + number.info.columns
      
      sample.table =
        matrix$matrix.data %>%
        dplyr::select(c(1:number.info.columns, start.col:end.col)) %>%
        # adding a column with the groups names
        dplyr::mutate(group = rep(group_names,
                                  times = c(group_start_row[-1] - group_start_row[-length(group_start_row)])))
      
      # Add the single table to a list
      samples.table.list[[i]] = sample.table
    }
    
    names(samples.table.list) = sample_names
    
    
    
    
    ######## Reorder samples in matrix
    reordered.matrix = matrix$matrix.data[,1:number.info.columns]
    
    for (i in new.sample.order) {
      reordered.matrix = cbind(reordered.matrix,
                               samples.table.list[[i]][,-c(1:number.info.columns,ncol(samples.table.list[[i]]))])
    }
    
    if (ncol(reordered.matrix) != ncol(matrix$matrix.data)){
      return(warning("Something went wrong and the number of original columns in the matrix does not correspond with the number of columns of the reordered one!!!"))
    }
    
    
    
    ######## Change sample_labels in the metadata
    new.header = gsub(pattern = paste0('[\"', paste0(names(samples.table.list), collapse = '\",\"'), '\"]'),
                      replacement = paste0('[\"', paste0(new.sample.order, collapse = '\",\"'), '\"]'),
                      x = paste0("@",data.table::fread(matrix.file, nrows = 1, stringsAsFactors = F, sep = "@", h = F)$V2),
                      fixed = T)
    
    if (nchar(original.header) != nchar(new.header) & isFALSE(ignore.header.error)){
      return(warning("Something went wrong with the replacement of the names in the header replacement!\nYou can ignore this by setting 'ignore.header.error = TRUE'."))
    }
    
    
    
    
    
    ######## Export the table and gzip it
    write(file = gsub(".gz", ".txt", reordered.matrix.path), x = new.header)
    
    data.table::fwrite(x = reordered.matrix,
                       file = gsub(".gz", ".txt", reordered.matrix.path),
                       quote = F, sep = "\t", row.names = F, col.names = F,
                       append = T)
    
    
    R.utils::gzip(filename = gsub(".gz", ".txt", reordered.matrix.path),
                  destname = reordered.matrix.path,
                  overwrite = T)
    
    
    ######## Print new order as message
    if (isTRUE(verbose)) {
      new_order = paste0(paste(paste0("   ",1:length(new.sample.order), "."), new.sample.order), collapse = "\n")
      message("The reorded matrix has been successfully exported with the following sample order:\n",
              new_order)
    }
  } # END function

