get.tss =
  function(bed,
           upstream = 1,
           downstream = 1) {
    
    bed_tb = as.data.frame(data.table::fread(bed))[,1:6]
    colnames(bed_tb) = c("chr", "start", "end", "name", "score", "strand")
    
    plus = dplyr::filter(bed_tb, strand == "+")
    minus = dplyr::filter(bed_tb, strand == "-")
    
    TSS_pos = rbind(data.frame(chr = plus$chr,
                               start = plus$start-upstream,
                               end = plus$start+downstream,
                               name = plus$name,
                               score = plus$score,
                               strand = plus$strand),
                    data.frame(chr = minus$chr,
                               start = minus$start-upstream,
                               end = minus$start+downstream,
                               name = minus$name,
                               score = minus$score,
                               strand = minus$strand))
    
    return(Rseb::sort.bed(TSS_pos))
  }
