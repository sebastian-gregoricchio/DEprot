# data.table


bed.overlap = function(reference.regions,
                       reference.regions.table.name = "referenceRegions",
                       target.regions,
                       export.table.file = NULL,
                       return.table = TRUE,
                       collapse.regions = FALSE,
                       verbose = TRUE) {
  
  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)
  #-----------------------------#
  
  # Check execution
  if (is.null(export.table.file) & isFALSE(return.table)) {
    return(warning("The parameters 'export.table.file' and 'return.table' are NULL and FALSE respectively.\nIn this way the output won't be saved. The execution is therefore interrupted."))
  }
  
  
  
  # Load reference.regions
  if ("character" %in% class(reference.regions)) {
    reference.regions.sorted = data.frame(data.table::fread(reference.regions))
    reference.regions.sorted = Rseb::sort.bed(reference.regions.sorted, return.bed = T)
    if (collapse.regions == T) {reference.regions.sorted = Rseb::collapse.bed(reference.regions.sorted, return.bed = T)}
  } else if ("data.frame" %in% class(reference.regions)) {
    reference.regions.sorted = Rseb::sort.bed(reference.regions, return.bed = T)
    if (collapse.regions == T) {reference.regions.sorted = Rseb::collapse.bed(reference.regions.sorted, return.bed = T)}
  } else {
    return(warning("The 'reference.regions' option must be either a character vector with the full path to the a bed file to load or a data.frames in at least BED3 format."))
  }
  
  colnames(reference.regions.sorted)[1:3] = c("chr", "start", "end")
  
  
  
  # Loading of target.regions
  if ("character" %in% class(target.regions)) {
    target.regions.sorted = list()
    for (i in 1:length(target.regions)){
      target.regions.sorted[[i]] = data.frame(data.table::fread(target.regions[i]))
      if (collapse.regions == T) {target.regions.sorted[[i]] = Rseb::collapse.bed(target.regions.sorted[[i]], return.bed = T)}
    }
  } else if ("list" %in% class(target.regions)) {
    target.regions.sorted = list()
    for (i in 1:length(target.regions)) {
      if (collapse.regions == T) {target.regions.sorted[[i]] = Rseb::collapse.bed(target.regions[[i]], return.bed = T)}
    }
  } else {
    return(warning("The 'target.regions' option must be either a character vector with the full path to the bed files to load or a list of data.frames."))
  }
  
  ## Assign names to target.regions
  if (is.null(names(target.regions))) {
    names(target.regions.sorted) = paste0("targetRegions.", LETTERS[1:length(target.regions.sorted)])
  } else {
    names(target.regions.sorted) = names(target.regions)
  }
  
  ## Merge target.regions tables
  target.regions.sorted = dplyr::bind_rows(target.regions.sorted, .id = "group.ID")
  target.regions.sorted = Rseb::move.df.col(target.regions.sorted, "group.ID last")
  target.regions.sorted = Rseb::sort.bed(target.regions.sorted, return.bed = T)
  
  colnames(target.regions.sorted)[1:3] = c("chr", "start", "end")
  
  
  # Merge and sort targets with reference
  all.regions.sorted = rbind(target.regions.sorted[,c(1:3,ncol(target.regions.sorted))],
                             dplyr::mutate(.data = reference.regions.sorted[,1:3],
                                           group.ID = "referenceRegions"))
  
  all.regions.sorted =
    dplyr::arrange(.data = all.regions.sorted, chr, start, end) %>%
    dplyr::mutate(rowNum = 1:nrow(all.regions.sorted))
  
  
  
  
  
} # END function