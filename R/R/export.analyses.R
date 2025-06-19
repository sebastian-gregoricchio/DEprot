#' @title export.analyses
#'
#' @description Export all the results from a DEprot.analyses object into a user-defined folder.
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param output.folder Path (string) to the output directory. If the folder is not already existing it will be created. Default: \code{"./export"}.
#' @param contrast.subset Numeric vector indicating the contrasts to use. Default: \code{NULL} (all contrasts are exported).
#' @param verbose Logical value to indicate whether progress messages should be printed. Default: \code{TRUE}.
#'
#' @export export.analyses


export.analyses =
  function(DEprot.analyses.object,
           output.folder = "./export",
           contrast.subset = NULL,
           verbose = TRUE) {


    ### libraries
    require(dplyr)


    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      warning("The input must be an object of class 'DEprot.analyses'.")
      return()
    }


    ### Subset contrasts
    if (!is.null(contrast.subset)) {
      if (all(contrast.subset %in% 1:length(DEprot.analyses.object@analyses.result.list))) {
        contrasts = contrast.subset
      } else {
        warning("Not all the contrasts indicated in the subset are present in the 'analyses.result.list' of the object provided.")
        return()
      }
    } else {
      contrasts = 1:length(DEprot.analyses.object@analyses.result.list)
    }


    ### make output dir
    outdir = gsub("/$", "", output.folder)
    dir.create(path = outdir, showWarnings = F, recursive = T)
    unlink(outdir, recursive = TRUE)
    dir.create(path = paste0(outdir,"/differential_analyses"), showWarnings = F, recursive = T)



    #################

    if (isTRUE(verbose)) {message("Exporting DEprot.analyses object...")}

    ## Export the full object
    saveRDS(object = DEprot.analyses.object,
            file = paste0(outdir,"/differential_analyses/differential_analyses.Rds"))



    ### Export metadata
    if (isTRUE(verbose)) {message("Exporting metadata...")}
    write.table(x = as.data.frame(DEprot.analyses.object@metadata),
                file = paste0(outdir,"/metadata.tsv"),
                sep = "\t", quote = F, col.names = T, row.names = F)


    ### Export counts
    if (isTRUE(verbose)) {message("Exporting all counts...")}
    dir.create(path = paste0(outdir,"/counts/boxplots"), showWarnings = F, recursive = T)

    if (!is.null(DEprot.analyses.object@raw.counts)) {
      write.table(x = as.data.frame(DEprot.analyses.object@raw.counts) %>% dplyr::mutate(prot.id = rownames(DEprot.analyses.object@raw.counts)) %>% dplyr::relocate(prot.id),
                  file = paste0(outdir,"/counts/raw_counts.tsv"),
                  sep = "\t", quote = F, col.names = T, row.names = F)

      pdf(file = paste0(outdir,"/counts/boxplots/raw_counts_boxplot.pdf"), height = 5, width = nrow(DEprot.analyses.object@metadata)*0.375)
      DEprot.analyses.object@boxplot.raw
      invisible(dev.off())
    }


    if (!is.null(DEprot.analyses.object@norm.counts)) {
      write.table(x = as.data.frame(DEprot.analyses.object@norm.counts) %>% dplyr::mutate(prot.id = rownames(DEprot.analyses.object@norm.counts)) %>% dplyr::relocate(prot.id),
                  file = paste0(outdir,"/counts/normalized_counts.tsv"),
                  sep = "\t", quote = F, col.names = T, row.names = F)

      pdf(file = paste0(outdir,"/counts/boxplots/normalized_counts_boxplot.pdf"), height = 5, width = nrow(DEprot.analyses.object@metadata)*0.375)
      DEprot.analyses.object@boxplot.norm
      invisible(dev.off())
    }


    if (!is.null(DEprot.analyses.object@imputed.counts)) {
      write.table(x = as.data.frame(DEprot.analyses.object@imputed.counts) %>% dplyr::mutate(prot.id = rownames(DEprot.analyses.object@imputed.counts)) %>% dplyr::relocate(prot.id),
                  file = paste0(outdir,"/counts/imputed_counts.tsv"),
                  sep = "\t", quote = F, col.names = T, row.names = F)

      pdf(file = paste0(outdir,"/counts/boxplots/imputed_counts_boxplot.pdf"), height = 5, width = nrow(DEprot.analyses.object@metadata)*0.375)
      DEprot.analyses.object@boxplot.imputed
      invisible(dev.off())
    }




    ########################

    # Export diff analyses parameters
    if (isTRUE(verbose)) {message("Exporting differential analyses parameters...")}
    params = DEprot.analyses.object@differential.analyses.params
    params$linear.FC.unresp.range = paste0("[",params$linear.FC.unresp.range[1]," ; ",params$linear.FC.unresp.range[2],"]")
    params.tb = data.frame(t(data.frame(params)))
    colnames(params.tb) = "value"
    params.tb =
      params.tb %>%
      dplyr::mutate(parameter = rownames(params.tb)) %>%
      dplyr::relocate(parameter)

    write.table(x = params.tb,
                file = paste0(outdir,"/differential_analyses/differential_analyses_parameters.tsv"),
                sep = "\t", quote = F, col.names = T, row.names = F)



    # Export diff analyses summary
    if (isTRUE(verbose)) {message("Exporting differential analyses summary...")}
    write.table(x = summary(DEprot.analyses.object),
                file = paste0(outdir,"/differential_analyses/differential_analyses_summary.n.differential.proteins.tsv"),
                sep = "\t", quote = F, col.names = T, row.names = F)


    # Export upset-plot
    if (isTRUE(verbose)) {message("Exporting an upset plot for all the contrasts requested...")}
    pdf(file = paste0(outdir,"/differential_analyses/differential_analyses_upsetplot.pdf"),
        width = 40, height = 10)
    print(DEprot::plot.upset(DEprot.analyses.object = DEprot.analyses.object, contrast.subset = contrasts))
    invisible(dev.off())



    ########################

    # Export diff analyses contrasts
    ## Set-up progress bar
    if (isTRUE(verbose)) {
      pb = progress::progress_bar$new(format = "(:spin) :bar :percent [Elapsed time: :elapsedfull || ETA: :eta]",
                                      total = length(contrasts),
                                      complete = "\u2588",   # Completion bar character
                                      incomplete = "\u2591", # Incomplete bar character
                                      current = "\u2592",    # Current bar character
                                      clear = T,        # If TRUE, clears the bar when finish
                                      width = 120)       # Width of the progress bar
    }


    if (isTRUE(verbose)) {message("Exporting contrasts...")}

    for (i in contrasts) {
      ## get contrast ID
      contrast.id = names(DEprot.analyses.object@contrasts)[i]

      ## make folder
      dir.create(path = paste0(outdir,"/differential_analyses/",contrast.id), showWarnings = F, recursive = T)


      ## Export diff table and n.diff
      write.table(x = DEprot.analyses.object@analyses.result.list[[i]]$results,
                  file = paste0(outdir,"/differential_analyses/",contrast.id,"/RESULTS_",contrast.id,".tsv"),
                  sep = "\t", quote = F, col.names = T, row.names = F)

      write.table(x = DEprot.analyses.object@analyses.result.list[[i]]$n.diff,
                  file = paste0(outdir,"/differential_analyses/",contrast.id,"/N.DIFF.PROTEINS_",contrast.id,".tsv"),
                  sep = "\t", quote = F, col.names = T, row.names = F)




      ## Export PCA and Correlation
      pdf(file = paste0(outdir,"/differential_analyses/",contrast.id,"/PCA.123_",contrast.id,".pdf"),
          width = 8, height = 8)
      print(DEprot.analyses.object@analyses.result.list[[i]]$PCA.plots)
      invisible(dev.off())

      saveRDS(object = DEprot.analyses.object@analyses.result.list[[i]]$PCA.data,
              file = paste0(outdir,"/differential_analyses/",contrast.id,"/PCA.DATA_",contrast.id,".Rds"))


      pdf(file = paste0(outdir,"/differential_analyses/",contrast.id,"/CORRELATION_",contrast.id,".pdf"),
          width = 16, height = 8)
      print(DEprot.analyses.object@analyses.result.list[[i]]$correlations)
      invisible(dev.off())




      ## Export volcano/MA plots
      pdf(file = paste0(outdir,"/differential_analyses/",contrast.id,"/VOLCANO_",contrast.id,".pdf"),
          width = 8, height = 5)
      print(DEprot.analyses.object@analyses.result.list[[i]]$volcano)
      invisible(dev.off())


      pdf(file = paste0(outdir,"/differential_analyses/",contrast.id,"/MA.PLOT_",contrast.id,".pdf"),
          width = 8, height = 5)
      print(DEprot.analyses.object@analyses.result.list[[i]]$MA.plot)
      invisible(dev.off())



      ## add tick to progress bar
      if (isTRUE(verbose)) {pb$tick()}
      } #end contrasts



    if (isTRUE(verbose)) {message("\nExport completed!")}


  } #END function
