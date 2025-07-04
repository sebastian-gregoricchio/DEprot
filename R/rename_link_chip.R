rename_link_chip =
  function(metadata.table,
           files_suffix = "_mdup.bam",
           output.dir,
           sample.file) {
    
    # Convert variables to absolute paths
    output.dir = tools::file_path_as_absolute(output.dir)
    metadata.tb = read.delim(metadata.table, sep = "\t", header = T)
    
    # Create the output.dir
    dir.create(output.dir, showWarnings = F, recursive = T)
    
    # Create the sample.config.file
    write(x = "chip_dict:", file = sample.file)
    
    for (i in 1:nrow(metadata.tb)){
      # Make symbolic link to rename the file in the output.dir
      R.utils::createLink(target = paste0(tools::file_path_as_absolute(metadata.tb$input_dir[i]), "/", metadata.tb$file_basename[i]),
                          link = paste0(output.dir, "/", metadata.tb$sample_ID[i], files_suffix),
                          overwrite = T)
      
      
      # Append the data to sample.config.file, if not an input file
      if (!((tolower(metadata.tb$corresponding_input_ID[i]) == "input") | (metadata.tb$corresponding_input_ID[i] == metadata.tb$sample_ID[i]))) {
        write(file = sample.file,
              append = T,
              paste0("  ", paste0(metadata.tb$sample_ID[i]), ":\n",
                     "    control: ", metadata.tb$corresponding_input_ID[i],"\n",
                     "    broad: ", stringr::str_to_sentence(metadata.tb$broad[i])))
      }
    }
  }