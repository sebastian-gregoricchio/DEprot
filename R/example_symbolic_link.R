#ATAC
fqd = "/home/s.gregoricchio/DATA_seb/Alex_endomitrium_project/HiC_tissues/fastq_files_allRuns/PE_GCF68002_NextSeq75bp"
link_dir = "/home/s.gregoricchio/DATA_seb/Alex_endomitrium_project/HiC_tissues/analyses_SG/00_fastQ"
GCF_num = 6812
fastq_list = list.files(fqd, full.names = T)

for (i in 1:length(fastq_list)) {
  new_name = gsub(pattern = paste0("(",GCF_num,"_[0-9]*_)|(_[A-Z]*-[A-Z]*_S[0-9]*)|_001"), replacement = "", x = fastq_list[i])
  new_name = gsub(pattern = fqd, replacement = link_dir, x = new_name)
  
  R.utils::createLink(link = new_name,
                      target = fastq_list[i],
                      overwrite=T)
}



# wo / at the end
fqd = "/home/s.gregoricchio/DATA_seb/Emma_project/CNT/00_fastq"
link_dir = "/home/s.gregoricchio/DATA_seb/Emma_project/CNT/00_runs"
GCF_num = 6927
fastq_list = list.files(fqd, full.names = T)
