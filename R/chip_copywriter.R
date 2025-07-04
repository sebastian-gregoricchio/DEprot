chip_copywriter =
  function(config_tb,  # data.frame: sample_id, target_id, peaks
           output.directory = paste0(getwd(), "/CNA_profiles_ChIP"),
           kb.resolution = 20,
           CNV.threshold = 2,
           genome = "hg38",
           excluded.chr = "",
           chr.prefix = "",
           cpu = 10,
           bedtools = "$HOME/.conda/envs/chip_zwart/bin/bedtools"){
    
    # Load libraries
    require(dplyr)
    require(ggplot2)
    
    
    # read config table
    if ("character" %in% class(config_tb)) {
      config = data.frame(data.table::fread(config_tb))
    } else if ("data.frame" %in% class(config_tb)) {
      config = data.frame(config_tb)
    } else {
      return(warning("The 'config_tb' must be either string with the path to a table or a data.frame with the following columns: sample_id target_id peaks."))
    }
    
    
    if (ncol(config) < 3 | colnames(config)[1:3] != c("sample_id", "target_id", "peaks")) {
      return(warning("The 'config_tb' must be either string with the path to a table or a data.frame with the following columns: sample_id target_id peaks."))
    }
    
    
    # Remove / at the end of the path if present
    data.folder = gsub("/$", "", output.directory)
    
    
    ### Set-up copywriteR
    bp.param = BiocParallel::SnowParam(workers = cpu, type = "SOCK")
    
    
    ### Generate genome index
    CopywriteR::preCopywriteR(output.folder = tools::file_path_as_absolute(file.path(data.folder)),
                              bin.size = kb.resolution*1000, #bp
                              ref.genome = genome,
                              prefix = chr.prefix)
    
    
    ### Run analyses
    CNA.plot_list = 
      purrr::pmap(.l = list(sample_id = config$sample_id,
                            target_id = config$target_id,
                            peaks = config$peaks),
                  .f = function(sample_id = sample_id, target_id = target_id, peaks = peaks){
                    
                    # Create output directory
                    sample_dir = paste0(data.folder, "/", sample_id)
                    if (dir.exists(sample_dir)) {system(paste0("rm -r ", sample_dir))}
                    dir.create(sample_dir, recursive = T, showWarnings = F)
                    
                    # Merge the peaks
                    system(paste0(bedtools, " merge -i ", peaks, " > ",
                                  sample_dir, "/", sample_id, "_collapsed.bed"))
                    
                    message(paste0("Dir: ", sample_dir))
                    
                    # Run CopyWriteR
                    CopywriteR::CopywriteR(sample.control = data.frame(target_id, target_id),
                                           destination.folder = sample_dir,
                                           reference.folder = paste0(data.folder, "/",genome,"_",kb.resolution,"kb"),
                                           capture.regions.file = paste0(sample_dir, "/", sample_id, "_collapsed.bed"),
                                           bp.param = bp.param,
                                           keep.intermediary.files = F)
                    
                    message(paste0(sample_id, ": computation done!"))
                    
                    
                    # Read log2_readCounts and plot it
                    log2.tb =
                      read.table(file = file.path(sample_dir, "CNAprofiles", "log2_read_counts.igv"), header = TRUE) %>%
                      dplyr::filter(!(Chromosome %in% paste0(chr.prefix, excluded.chr))) %>%
                      dplyr::mutate(mid.point = (Start+End)/2,
                                    Chromosome = factor(Chromosome, levels = unique(Chromosome)))
                    names(log2.tb)[5] = "log2.value"
                    
                    CNA.plot =
                      ggplot(log2.tb,
                             aes(x = mid.point,
                                 y = log2.value)) +
                      geom_point(size = 0.001) +
                      theme_classic() +
                      theme(axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.ticks.y = element_line(color = "black"),
                            axis.text.y = element_text(color = "black"),
                            panel.spacing = unit(0,'lines'),
                            strip.background = element_blank(),
                            strip.placement = "outside",
                            panel.grid = element_blank(),
                            panel.border = element_blank()) +
                      geom_hline(yintercept = c(-1,1)*log2(CNV.threshold), linetype = 2, color = "indianred") +
                      geom_vline(xintercept = +Inf, linetype = 1, color = "gray70") +
                      ylim(c(-1,1)*4) +
                      xlab("Chromosome") +
                      ylab("log2(copy number)") +
                      ggtitle(sample_id) +
                      facet_grid(~ Chromosome,
                                 space = "free_x",
                                 scales = "free_x",
                                 switch = "x")
                    
                    pdf(paste0(sample_dir, "/CNA.plot_", sample_id, "_all.chr.pdf"), width = 16, height = 5)
                    print(CNA.plot)
                    invisible(dev.off())
                    
                    
                    # Plot by copyWriteR
                    CopywriteR::plotCNA(destination.folder = sample_dir)
                    
                    
                    # Export filtered table
                    log2.tb_filtered =
                      read.table(file = file.path(sample_dir, "CNAprofiles", "log2_read_counts.igv"), header = TRUE) %>%
                      dplyr::filter(abs(.[[5]]) >= log2(CNV.threshold))
                    
                    write.table(x = log2.tb_filtered,
                                file = file.path(sample_dir, "CNAprofiles", "log2_read_counts_filtered.igv"),
                                sep = "\t", quote = F, row.names = F, col.names = T)
                    
                    
                    return(CNA.plot)
                  })
    
    return(CNA.plot_list)
  } # END function