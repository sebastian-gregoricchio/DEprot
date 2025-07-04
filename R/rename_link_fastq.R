rename.link.fastq = 
  function(input.dir,
           output.dir,
           fastq.extension = ".fastq.gz",
           pattern.input.files = ".gz") {
    
    # Clean directory names
    input.dir = gsub("/$", "", tools::file_path_as_absolute(input.dir))
    output.dir = gsub("/$", "", tools::file_path_as_absolute(output.dir))
    fastq.extension.final = gsub("^[.]","",fastq.extension)
    fastq.extension.sub = gsub(".","[.]",fastq.extension.final, fixed = T)
    
    # List files
    files = list.files(path = input.dir,
                       pattern = pattern.input.files,
                       recursive = F,
                       full.names = F,
                       include.dirs = F)
    
    renamed.files = sapply(files,
                           function(x){
                             x = gsub("[0-9]*_[0-9]*_","",x)
                             x = gsub("_[A-Z]*-[A-Z]*_S[0-9]*_R", "_R", x)
                             x = gsub(paste0("_[0-9]*[.]",fastq.extension.sub, "$"),
                                      paste0(".",fastq.extension.final),
                                      x)},
                           USE.NAMES = F)
    
    # Make the renamed_link
    for (i in 1:length(renamed.files)) {
      system(paste0("ln -s ",
                    input.dir,"/",files[i]," ",
                    output.dir,"/",renamed.files[i]))
    }
  }