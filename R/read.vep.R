read.vep =
  function(vep.vcf,
           output.file = NULL,
           return.df = TRUE) {
    
    if ((!is.null(output.file)) & (return.df == T)) {return()}
    
    require(vcfR)
    require(ensemblVEP)
    require(dplyr)
    
    vep_file = "VCF_vep.vcf"
    
    vep_vcfr = read.vcfR(vep.vcf, verbose = FALSE)
    vep_header = data.frame(vep_vcfr@meta)
    vep_variants = data.frame(vep_vcfr@fix)
    vep_gt = data.frame(vep_vcfr@gt)
    
    # Parse into a GRanges and include the 'VCFRowID' column.
    vep_ens = readVcf(vep_file, "hg19")
    csq_vep = parseCSQToGRanges(vep_ens)
    csq_vep = data.frame(csq_vep)
    
    
    VEP = cbind.data.frame(vep_variants,csq_vep,vep_gt)
    
    if (!is.null(output.file)) {
      write.table(x = VEP,
                  file = output.file,
                  col.names = T, row.names = F, sep = "\t", quote = F)
    }
    
    if (return.df == T) {return(VEP)}
  }