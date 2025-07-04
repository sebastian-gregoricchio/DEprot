sample.sheet = read.delim("/home/s.gregoricchio/DATA_seb/Organoids/ChIP/sample_sheets/7089_sample_sheet_ChIP.organoids.txt")


rename.chip.cram.to.bam(sample.sheet = sample.sheet,
                        output.folder = "/home/s.gregoricchio/DATA_seb/Organoids/ChIP/H3K27ac/01_BAM",
                        sample.config.file.name = "/home/s.gregoricchio/DATA_seb/Organoids/ChIP/sample_sheets/7089_sample_config.yaml")



rename.chip.cram.to.bam = 
  function(sample.sheet = NULL,
           output.folder,
           sample.ids = if(("character" %in% class(sample.sheet)) & !is.null(sample.sheet)) {read.delim(sample.sheet)[,1]} else {sample.sheet[,1]},
           cram.paths = if(("character" %in% class(sample.sheet)) & !is.null(sample.sheet)) {read.delim(sample.sheet)[,2]} else {sample.sheet[,2]},
           corresponding.input = if(("character" %in% class(sample.sheet)) & !is.null(sample.sheet)) {read.delim(sample.sheet)[,3]} else {sample.sheet[,3]},
           broad = if(("character" %in% class(sample.sheet)) & !is.null(sample.sheet)) {read.delim(sample.sheet)[,4]} else {sample.sheet[,4]},
           reference.genome = "/home/s.gregoricchio/annotations/genomes/Hg38/one_file_fasta/Homo_sapiens.GRCh38_v102.dna.primary_assembly.fa",
           bam.suffix = "_mdup.bam",
           sample.config.file.name = paste0(gsub("/$", "", dirname(output.folder)), "/sample_config.yaml"),
           threads = 32,
           samtools.command = "/home/s.gregoricchio/.conda/envs/NGS/bin/samtools") {
    
    # Check output.folder
    output.folder = gsub("/$", "", output.folder)
    dir.create(output.folder, showWarnings = F)
    parent.folder = gsub("/$", "", dirname(output.folder))
    
    # Compose the data.frame
    metadata = data.frame(id = sample.ids,
                          cram = cram.paths,
                          input = corresponding.input,
                          broad = broad,
                          renamed.bam = paste0(output.folder,"/",sample.ids,bam.suffix))
    
    input.list = sample.ids[grep(paste0(unique(metadata$input), collapse = "|"), sample.ids)]
    names(input.list) = unique(metadata$input)
    
    
    # carm-to-bam conversion
    for (i in 1:nrow(metadata)) {
      message(paste0("\n>>>>>>>>> [",i,"/",nrow(metadata),"]: ",metadata$id[i]))
      
      # Conversion cram to bam
      system(paste(samtools.command, "view",
                   "-@", threads,
                   "-b -t", reference.genome,
                   "-o", metadata$renamed.bam[i],
                   metadata$cram[i]))
      
      # Indexing
      system(paste(samtools.command, "index",
                   "-@", threads,
                   metadata$renamed.bam[i]))
    }
    
    
    # generate sample config file for snakepipes
    write(x = "chip_dict:", file = sample.config.file.name)
    
    for (i in 1:nrow(metadata)) {
      if (!(metadata$id[i] %in% input.list)) {
        write(x = paste0("  ", metadata$id[i],":"),
              file = sample.config.file.name,
              append = T)
        write(x = paste0("    broad: ", tolower(as.character(metadata$broad)[i])),
              file = sample.config.file.name,
              append = T)
        write(x = paste0("    control: '", input.list[grep(metadata$input[i], input.list)], "'"),
              file = sample.config.file.name,
              append = T) 
      }
    }
    
    message(paste0("*** Sample configuration file written in ->> ", tools::file_path_as_absolute(sample.config.file.name), " <<- ***"))
  }