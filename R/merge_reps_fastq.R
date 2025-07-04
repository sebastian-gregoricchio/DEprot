organoids.merge.fastq.reps(samples = c("ITO68T", "ITO124T"),
                           input.folder = "/home/s.gregoricchio/DATA_seb/Organoids/ATAC/00_runs/",
                           output.folder = "/home/s.gregoricchio/DATA_seb/Organoids/ATAC_merged_reps/00_runs/")

organoids.merge.fastq.reps =
  function(samples, input.folder, output.folder) {
    
    input.folder = gsub("/$", "", input.folder)
    output.folder = gsub("/$", "", output.folder)
    dir.create(path = paste0(output.folder, "/"), showWarnings = F, recursive = T)
    
    for (i in 1:length(samples)) {
      message(paste0("[", i, "/", length(samples), "] Merging sample: ", samples[i]))
      
      R1.list = paste(list.files(path = input.folder, pattern = paste0(samples[i], ".*_R1.fastq.gz"), full.names = T), collapse = " ")
      R2.list = paste(list.files(path = input.folder, pattern = paste0(samples[i], ".*_R2.fastq.gz"), full.names = T), collapse = " ")
      
      system(paste("zcat", R1.list, ">", paste0(output.folder, "/", samples[i], "_ATAC_mergedReps_R1.fastq")))
      system(paste("zcat", R2.list, ">", paste0(output.folder, "/", samples[i], "_ATAC_mergedReps_R2.fastq")))
      
      # R.utils::gzip(paste0(output.folder, "/", samples[i], "_mergedReps_R1.fastq"), overwrite = T, remove = T)
      # R.utils::gzip(paste0(output.folder, "/", samples[i], "_mergedReps_R2.fastq"), overwrite = T, remove = T)
    }
    system(paste("gzip", paste0(output.folder, "/*_ATAC_mergedReps_R*.fastq")))
  }
