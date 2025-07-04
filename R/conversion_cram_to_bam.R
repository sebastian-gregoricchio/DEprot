cram_dir = "/shared/gcf/s.gregoricchio/6811/cram_files/"
bam_dir = "/home/s.gregoricchio/DATA_seb/Organoids/ChIP/H3K27ac/01_BAM/"
ref_genome = "/home/s.gregoricchio/annotations/genomes/Hg38/one_file_fasta/Homo_sapiens.GRCh38_v102.dna.primary_assembly.fa"

cram_list = list.files(cram_dir, pattern = ".cram$")

bam_basename = gsub(pattern = "^[0-9]*_[0-9]*_wz[0-9]*_ITO[-]", replacement = "ITO", x = cram_list)
bam_basename = gsub(pattern = "_[A-Z]*[-][A-Z]*_S[0-9]*[-]mdup.cram", replacement = "_mdup", x = bam_basename)

dir.create(bam_dir, showWarnings = F, recursive = T)

for (i in 1:length(cram_list)) {
  # conversion to bam
  system(paste0("samtools view -@ 4 -b -T ", ref_genome,
                " -o ", bam_dir,bam_basename[i],".bam ",
                cram_dir,cram_list[i]))
  
  system(paste0("samtools index -@ 4 -b ", bam_dir,bam_basename[i],".bam "))
}


# ----------------------------------------------------------
# not organoids

cram_dir = "/shared/gcf/n.eickhoff/7166/cram_files/"
bam_dir = "/home/s.gregoricchio/DATA_seb/Nils_project/GCF-7166/bam_files/"
ref_genome = "/home/s.gregoricchio/annotations/genomes/Hg38/one_file_fasta/Homo_sapiens.GRCh38_v102.dna.primary_assembly.fa"

cram_list = list.files(cram_dir, pattern = ".cram$")

bam_basename = gsub(pattern = "^[0-9]*_[0-9]*_", replacement = "", x = cram_list)
bam_basename = gsub(pattern = "_[A-Z]*[-][A-Z]*_S[0-9]*[-]mdup.cram", replacement = "_mdup", x = bam_basename)

dir.create(bam_dir, showWarnings = F, recursive = T)

for (i in 1:length(cram_list)) {
  # conversion to bam
  system(paste0("samtools view -@ 4 -b -T ", ref_genome,
                " -o ", bam_dir,bam_basename[i],".bam ",
                cram_dir,cram_list[i]))
  
  system(paste0("samtools index -@ 4 -b ", bam_dir,bam_basename[i],".bam "))
}