#' @title Multi-bigWig average tool
#' 
#' @description This tools uses \code{bigwigCompare} (see details) to perform the average between multiple bigWig files in a multi-step process.
#'
#' @param sample.config A data.frame or a string indicating the path to a table. The table/data.frame must contain in the first column the file list, while the second column the group to which belong each file. Files belonging to the same group will be merged together. The group name corresponds to the name of the corresponding output file.
#' @param output.dir String of the path to the output directory. Default: \code{file.path(getwd(), "merged_bigWigs")}.
#' @param bin.size Size of the bins in base-pairs (bp). Default: \code{50}.
#' @param merged.suffix Suffix to use for the resulting output files: <group.id><suffix><extension>. Default: \code{paste0("_merged_bs", format(ceiling(bin.size), scientific = FALSE))}.
#' @param out.file.format Output file format. Possible choices: "bigWig", "bw", "bedGraph", "bdg". Default: \code{"bigwig"}.
#' @param number.of.processors Number of CPUs to use. Default: \code{4}.
#' @param bigwigCompare String of the path to the \code{bigwigCompare} tool. Default: \code{"~/.conda/envs/snakepipes/b9364eb954bd13fe4e5f104fa8a286e6/bin/bigwigCompare"}.
#' @param return.log Logic value to define whether the log list of the \code{bigwigCompare} steps performed by 'group.id'. Default: \code{TRUE}.
#' @param verbose Logic value to define whether the function should print messages. Default: \code{TRUE}.
#'
#' @return A list of vectors with the commands used to run \code{bigwigCompare} at each step. Each element is one group ID.
#' 
#' @details For details on \code{bigwigCompare} visit: \url{https://deeptools.readthedocs.io/en/develop/content/tools/bigwigCompare.html}
#' 
#' @export multi.bigwig.mean

multi.bigwig.mean =
  function(sample.config,
           output.dir = file.path(getwd(), "merged_bigWigs"),
           bin.size = 50,
           merged.suffix = paste0("_merged_bs", format(ceiling(bin.size), scientific = FALSE)),
           out.file.format = "bigwig",
           number.of.processors = 4,
           bigwigCompare = "~/.conda/envs/chip_zwart/bin/bigwigCompare",
           return.log = FALSE,
           verbose = TRUE) {
    
    start.time = Sys.time()
    
    # Load sample.config
    if ("character" %in% class(sample.config)) {
      config = data.frame(data.table::fread(sample.config))
    } else if ("data.frame" %in% class(sample.config)) {
      config = data.frame(sample.config)
    } else {
      return(warning("The 'sample.config' must be either a data.frame or the path to a table."))
    }
    colnames(config)[1:2] = c("file", "group")
    config = dplyr::mutate(config, file = sapply(config$file, tools::file_path_as_absolute))
    
    
    # Create output.folder if does not exits => `mkdir -p`
    dir.create(path = output.dir, recursive = T, showWarnings = F)
    output.dir = gsub("/$", "", tools::file_path_as_absolute(output.dir))
    if (isTRUE(verbose)) {message(paste(">>> output folder:", output.dir))}
    
    
    # generate log file
    log.file = file.path(output.dir, "multi.bigwig.mean.log")
    write(x = paste("----------------------------------------------\n Multi-bigWig average ///", Sys.time(), "\n by Sebastian Gregoricchio\n"), file = log.file, append = F)
    system(paste0(bigwigCompare, " --version >> ", file.path(output.dir, "multi.bigwig.mean.log")))                
    write(x = paste("----------------------------------------------\n\n\n*** Sample-group configuration table:\n"), file = log.file, append = T)
    suppressWarnings(write.table(x = config, file = log.file, append = T, col.names = T, row.names = F, quote = F, sep = "\t"))
    
    
    # Collect the group names
    group.ID = unique(config$group)
    
    
    # Decide file extension
    file.ext = ifelse(test = tolower(out.file.format) %in% c("bigwig", "bw"),
                      yes = ".bw", no = ".bdg")
    out.format = ifelse(test = tolower(out.file.format) %in% c("bigwig", "bw"),
                        yes = "bigwig", no = "bedgraph")
    
    
    # Perform the average of the bigwigs and generate log
    log.vector = list()
    write(x = paste("\n\n\n*** Command line for bigwigCompare:\n"), file = log.file, append = T)
    
    for (i in 1:length(group.ID)) {
      write(x = paste0("\n>>> [",i,"/",length(group.ID),"]: ", group.ID[i]), file = log.file, append = T)
      files.to.merge = (dplyr::filter(config, group == group.ID[i]))$file
      
      # 1 file
      if (length(files.to.merge) == 1) {
        if (isTRUE(verbose)) {message(paste0("[",i,"/",length(group.ID),"]: ",group.ID[i], " (only ", length(files.to.merge), " file -> averaging with itself)"))}
        command = paste(bigwigCompare,
                        "--bigwig1", files.to.merge[1],
                        "--bigwig2", files.to.merge[1],
                        "--operation mean",
                        "--skipZeroOverZero",
                        "--binSize", format(ceiling(bin.size), scientific = FALSE),
                        "--numberOfProcessors", number.of.processors,
                        "--outFileFormat", out.format,
                        "--outFileName", file.path(output.dir, paste0(group.ID[i], merged.suffix, file.ext)))
        log.vector[[i]] = command
        write(x = paste0("  * (step.1) --> ", command), file = log.file, append = T)
        system(command)
      }
      
      # 2 files
      else if (length(files.to.merge) == 2) {
        if (isTRUE(verbose)) {message(paste0("[",i,"/",length(group.ID),"]: ",group.ID[i], " (merging ", length(files.to.merge), " files)"))}
        command = paste(bigwigCompare,
                        "--bigwig1", files.to.merge[1],
                        "--bigwig2", files.to.merge[2],
                        "--operation mean",
                        "--skipZeroOverZero",
                        "--binSize", format(ceiling(bin.size), scientific = FALSE),
                        "--numberOfProcessors", number.of.processors,
                        "--outFileFormat", out.format,
                        "--outFileName", file.path(output.dir, paste0(group.ID[i], merged.suffix, file.ext)))
        log.vector[[i]] = command
        write(x = paste0("  * (setp.1) --> ", command), file = log.file, append = T)
        system(command) 
      }
      
      # 3+ files
      else {
        if (isTRUE(verbose)) {message(paste0("[",i,"/",length(group.ID),"]: ",group.ID[i], " (merging ", length(files.to.merge), " files)"))}
        scaling.factor = 1/length(files.to.merge)
        n.loops = length(files.to.merge)-1
        params.tb = data.frame(bw1 = files.to.merge[1:n.loops],
                               bw2 = c(files.to.merge[n.loops+1],
                                       file.path(output.dir, paste0("temp_",group.ID[i],"_loop",1:(n.loops-1), file.ext))),
                               out = c(file.path(output.dir, paste0("temp_",group.ID[i],"_loop",1:(n.loops-1), file.ext)),
                                       file.path(output.dir, paste0(group.ID[i], merged.suffix, file.ext))),
                               sf.bw1 = rep(scaling.factor, n.loops),
                               sf.bw2 = c(scaling.factor, rep(1, (n.loops-1))))
        
        multi.log = c()
        
        for (k in 1:n.loops) {
          if (isTRUE(verbose)) {message(paste0("   > ",k," of ", n.loops," steps --> ", basename(params.tb$out[k])))}
          command = paste(bigwigCompare,
                          "--bigwig1", params.tb$bw1[k],
                          "--bigwig2", params.tb$bw2[k],
                          "--operation add",
                          "--scaleFactors", paste0(params.tb$sf.bw1[k],":",params.tb$sf.bw2[k]),
                          "--skipZeroOverZero",
                          "--binSize", format(ceiling(bin.size), scientific = FALSE),
                          "--numberOfProcessors", number.of.processors,
                          "--outFileFormat", out.format,
                          "--outFileName", params.tb$out[k])
          multi.log = c(multi.log, command)
          write(x = paste0("  * (step.",k,") --> ", command, "\n"), file = log.file, append = T)
          system(command)
        }
        
        log.vector[[i]] = multi.log
        write(x = "\n", file = log.file, append = T)
        
        # removing temp files
        if (isTRUE(verbose)) {message(paste0("   > removing temp files"))}
        system(paste("rm", file.path(output.dir, paste0("temp_",group.ID[i],"*"))))
      }
    }
    
    names(log.vector) = group.ID
    write(x = paste0("\n\n\nRun time: ", round((Sys.time()-start.time)[[1]],3), " min"), file = log.file, append = T)
    
    if (isTRUE(return.log)) {return(log.vector)}
  } #END function

