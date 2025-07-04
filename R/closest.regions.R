# data.table


closest.regions = function(reference.regions,
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
  if ("data.frame" %in% class(target.regions)) {
    target.regions = list(target.regions)
  }
    
  if ("character" %in% class(target.regions)) {
    target.regions.sorted = list()
    for (i in 1:length(target.regions)){
      target.regions.sorted[[i]] = data.frame(data.table::fread(target.regions[i]))
      if (collapse.regions == T) {target.regions.sorted[[i]] = Rseb::collapse.bed(target.regions.sorted[[i]], return.bed = T)}
    }
  } else if ("list" %in% class(target.regions)) {
    target.regions.sorted = list()
    for (i in 1:length(target.regions)) {
      if (collapse.regions == T) {target.regions.sorted[[i]] = Rseb::collapse.bed(target.regions[[i]], return.bed = T)
      } else {
        target.regions.sorted[[i]] = data.frame(target.regions[i])
      }
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
                             dplyr::mutate(.data = reference.regions.sorted[,1:3], group.ID = "referenceRegions"))
  
  all.regions.sorted =
    dplyr::arrange(.data = all.regions.sorted, chr, start, end) %>%
    dplyr::mutate(rowNum = 1:nrow(all.regions.sorted))
  
  
  
  # Detect the regions before and after each referenceRegion
  final.table = data.frame()
  cubes.shown = -1
  
  for (i in 1:nrow(reference.regions.sorted)) {
    # Print progression bar if required
    if ((verbose == T)) {
      progression = i/nrow(reference.regions.sorted)
      perc.progression = Rseb::floating.floor(progression*100, 1)
      n.cubes = 75
      cubes.to.show = floor(progression*n.cubes)
      
      if (cubes.to.show > cubes.shown) {
        message((paste0("|", paste0(rep("\u25A0", cubes.to.show), collapse = ""),
                        paste0(rep("-", (n.cubes - cubes.to.show)), collapse = ""),
                        "| ", perc.progression, "%")))
        cubes.shown = cubes.to.show
      }
    }
    
    
    # identification of the position of a reference.region in the all.regions.table
    ref.position = (all.regions.sorted %>% dplyr::filter(chr == reference.regions.sorted$chr[i],
                                                         start == reference.regions.sorted$start[i],
                                                         end == reference.regions.sorted$end[i],
                                                         group.ID == "referenceRegions"))$rowNum[1]
    
    # Get upstream region position
    if (ref.position > 1) {
      continue = TRUE
      k = ref.position
      while (continue == TRUE & k>1) {
        k = k-1
        if ((all.regions.sorted$group.ID[k] != "referenceRegions") & (all.regions.sorted$chr[k] == reference.regions.sorted$chr[i])) {
          continue = FALSE
          upstream.region.pos = k
        } else if ((k-1) > 1) {
          if (all.regions.sorted$chr[k-1] == reference.regions.sorted$chr[i]) {continue = TRUE}
        } else {
          upstream.region.pos = NA
          continue = FALSE
        }
      }#w
    } else {
      upstream.region.pos = NA
      continue = FALSE
    } # end upstream
    
    
    # Get downstream region position
    if (ref.position < nrow(all.regions.sorted)) {
      continue = TRUE
      k = ref.position
      while ((continue == TRUE) & (k < nrow(all.regions.sorted))) {
        k = k+1
        if ((all.regions.sorted$group.ID[k] != "referenceRegions") & (all.regions.sorted$chr[k] == reference.regions.sorted$chr[i])) {
          continue = FALSE
          downstream.region.pos = k
        } else if ((k+1) < nrow(all.regions.sorted)) {
          if (all.regions.sorted$chr[k+1] == reference.regions.sorted$chr[i]) {continue = TRUE}
        } else {
          downstream.region.pos = NA
          continue = FALSE
        }
      }#w
    } else {
      downstream.region.pos = NA
      continue = FALSE
    } # end downstream
    
    
    
    # set NA position to nrow+1 in order to get NA values for the final table (next step)
    if (is.na(upstream.region.pos)) {upstream.region.pos = nrow(all.regions.sorted)+1}
    if (is.na(downstream.region.pos)) {downstream.region.pos = nrow(all.regions.sorted)+1}
    
    
    # Get refRegion + upstream + downstream
    current.reference = data.frame(chr.ref = all.regions.sorted$chr[ref.position],
                                   start.ref = all.regions.sorted$start[ref.position],
                                   end.ref = all.regions.sorted$end[ref.position],
                                   reference.group = "referenceRegions",
                                   
                                   chr.upstream = all.regions.sorted$chr[upstream.region.pos],
                                   start.upstream = all.regions.sorted$start[upstream.region.pos],
                                   end.upstream = all.regions.sorted$end[upstream.region.pos],
                                   upstream.group = all.regions.sorted$group.ID[upstream.region.pos],
                                   
                                   chr.downstream = all.regions.sorted$chr[downstream.region.pos],
                                   start.downstream = all.regions.sorted$start[downstream.region.pos],
                                   end.downstream = all.regions.sorted$end[downstream.region.pos],
                                   downstream.group = all.regions.sorted$group.ID[downstream.region.pos])
    
    
    # Append current reference region to general table
    final.table = rbind(final.table, current.reference)
  } # FOR end
  
  if (nrow(final.table) != nrow(reference.regions.sorted)) {
    return(warning("Something went wrong: the number of regions in the resulting table is different from the original reference one."))
  }
  
  
  
  # Export table if required
  if (!is.null(export.table.file)) {
    write.table(x = final.table, file = export.table.file,
                sep = "\t", quote = F, col.names = T, row.names = F)
  }
  
  if (return.table == T) {
    return(final.table)
  } else {return(invisible(NULL))}
  
} # END function