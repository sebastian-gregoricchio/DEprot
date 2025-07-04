flowcell.1 = "/shared/gcf/s.gregoricchio/7083/NXT289_fastq_files/"
flowcell.2 = "/shared/gcf/s.gregoricchio/7083/NXT290_fastq_files/"
output.folder = "/home/s.gregoricchio/DATA_seb/Endometrium_vs_Breast_project/ATAC/00_fastq/"
pattern = "Ishikawa|MCF7"



merge.flowcells.nki =
  function(flowcell.1,
           flowcell.2,
           output.folder,
           pattern = "") {
    
    fastq.1 = list.files(flowcell.1, pattern = pattern)
    fastq.2 = list.files(flowcell.2, pattern = pattern)
    if (length(unique(fastq.1 %in% fastq.2)) > 1) {
      return(warning("The fastq names in the two flowcells do not correspond"))
    }
    
    for (i in 1:length(fastq.1)) {
      name = fastq.1[i]
      file.1 = paste0(gsub("/$","/",flowcell.1), name)
      file.2 = paste0(gsub("/$","/",flowcell.2), name)
      file.out = paste0(gsub("/$","/",output.folder), name)
      
      #system(paste("cat", file.1, file.2, ">", file.out))
      
      system(paste("zcat", file.1, file.2, "| gzip --fast >", file.out))
      
    }
  }


merge.flowcells.nki(flowcell.1, flowcell.2, output.folder, pattern)