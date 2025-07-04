# TOBIAS ATAC aggregates analyses

analyze.TF.activity =
  function(sample.groups,  # sample.ID | group
           aggregate.scores.dir,
           comparison.list = NULL,
           scale.zero.one = FALSE,
           paired.stats = FALSE,
           ignore.statistics = FALSE,
           pseudocount = 0.00000000001,
           compute.confidence.interval = T,
           confidence.level = 0.9,
           foldChange.threshold = 1.3,
           percentage.top.hits = 5,
           basemean.AUC.threshold = 0.1,
           p.adjust.method = "bonferroni",  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
           heatmap.cluster.rows = TRUE,
           heatmap.cluster.columns = TRUE,
           heatmap.scaling.mode = "row", # row, column, none
           heatmap.palette = "Spectral",
           heatmap.nColors = 101,
           heatmap.showRowNames = FALSE,
           heatmap.showColNames = TRUE,
           heatmap.n.row.kmeans = NA,
           heatmap.border.color = NA,
           aggregate.table.suffix = "_tobias_aggregate_score_table.txt",
           group.colors = NULL,
           group.order = NULL,
           error.transparency = 0.15,
           MAplot.show.names = TRUE,
           TF_whitelist = NULL,
           TF_blacklist = NULL) {
    
    #############################################
    # Load libraries
    require(dplyr)
    require(plyr)
    require(ggplot2)
    require(pryr)
    
    #############################################
    
    
    # Clean input data
    aggregate.scores.dir = gsub("/$", "", aggregate.scores.dir)
    
    colnames(sample.groups)[1:2] = c("sample.ID", "group")
    if (is.null(group.order)) {levels = unique(sample.groups$group)} else {levels = group.order}
    sample.groups =
      sample.groups %>%
      dplyr::mutate(group = factor(group, levels = levels)) %>%
      dplyr::arrange(group, sample.ID)
    
    
    # get the list of tables and corresponding TF name
    score_file_list = list.files(path = aggregate.scores.dir,
                                 pattern = aggregate.table.suffix,
                                 full.names = T)
    
    TF_list = gsub(paste0(aggregate.scores.dir,"/|",aggregate.table.suffix), "", score_file_list)
    
    
    # Filtering user-defined TFs
    if (!is.null(TF_whitelist)) {TF_list = TF_list[TF_list %in% TF_whitelist]}
    if (!is.null(TF_blacklist)) {TF_list = TF_list[!(TF_list %in% TF_blacklist)]}
    
    
    # Define colors for groups in footprint
    if (is.null(group.colors)) {
      group.colors = rainbow(length(levels))
    }
    
    
    # Define comparisons if not user-defined
    if (is.null(comparison.list)) {
      combos = combn(x = unique(sample.groups$group), m = 2)
      comparison.list = list()
      for (i in 1:ncol(combos)) {
        comparison.list[[i]] = c(as.vector(combos[1,i]), as.vector(combos[2,i]))
      }
    }
    
    
    # Check if groups in comparison
    all_groups_in_comparison = c()
    for (i in 1:length(comparison.list)) {
      all_groups_in_comparison = c(all_groups_in_comparison, comparison.list[[i]])
    }
    
    if (FALSE %in% unique(unique(all_groups_in_comparison) %in% unique(sample.groups$group))) {
      return(warning("The groups defined in the 'comparison.list' do not match the ones in the 'sample.groups'."))
    }
      
    
    
    # ------------------------------------------------------------------------------------------------------
    ### Read score tables and run foorprint analyses
    message("Read score tables...")
    
    score_list = list()
    
    
    score_list =
      purrr::pmap(.l = list(TF = TF_list),
                .f = function(TF) {
                  score_list =
                    data.frame(data.table::fread(paste0(aggregate.scores.dir,"/",TF,aggregate.table.suffix), skip = 2) %>%
                                 tidyr::separate(col = V1, sep = "\t", into = c("sample.ID", "TF", "V1")))
                  
                  score_list =
                    dplyr::left_join(score_list,
                                     sample.groups,
                                     by = "sample.ID") %>%
                    dplyr::mutate(group = factor(group, levels = levels))
                  
                  return(score_list)
                })
    
    
    
    
    for (i in 1:length(TF_list)) {
      score_list[[i]] =
        data.frame(data.table::fread(paste0(aggregate.scores.dir,"/",TF_list[i],aggregate.table.suffix), skip = 2) %>%
                     tidyr::separate(col = V1, sep = "\t", into = c("sample.ID", "TF", "V1")))
      
      score_list[[i]] =
        dplyr::left_join(score_list[[i]],
                         sample.groups,
                         by = "sample.ID") %>%
        dplyr::mutate(group = factor(group, levels = levels))
    }
    names(score_list) = TF_list
    
    
    
    ### Run footprint analyses and plots by group
    reshaped_score_list = list()
    reshaped_score_mean_list = list()
    footprint_list = list()
    
    for (i in 1:length(TF_list)) {
      # Reshaping score table
      reshaped_score_list[[i]] = data.frame()
      flanking = ceiling((ncol(score_list[[i]]) - 3)/2)
      
      for (k in 1:nrow(score_list[[i]])) {
        reshaped_score_list[[i]] = rbind(reshaped_score_list[[i]],
                                         data.frame(sample.ID = score_list[[i]]$sample[k],
                                                    group = as.vector(score_list[[i]]$group[k]),
                                                    TF = TF_list[i],
                                                    distance_center_bp = c(-flanking:flanking)[-(flanking+1)],
                                                    score = as.numeric(t(score_list[[i]][k,3:(ncol(score_list[[i]])-1)])[,1]) + pseudocount))
      }
      
      # bring to zero if required
      message("Scale between 0 and 1 the scores per TF...")
      if (scale.zero.one == TRUE) {
        reshaped_score_list[[i]] =
          reshaped_score_list[[i]] %>%
          dplyr::mutate(score = score - min(score)) %>%
          dplyr::mutate(score = score / max(score))
      }
      
      # Compute mean+stats by group at each position
      message("Compute mean+stats by group at each position...")
      
      reshaped_score_mean_list[[i]] =
        data.frame(reshaped_score_list[[i]] %>%
                     dplyr::group_by(distance_center_bp, group, TF) %>%
                     dplyr::summarise(n = n(),
                                      score_mean = mean(score, na.rm = T),
                                      SD = sd(score, na.rm = T),
                                      SEM = sd(score, na.rm = T)/sqrt(n()),
                                      .groups = "keep") %>%
                     dplyr::mutate(group = factor(group, levels = levels)) %>%
                     dplyr::arrange(group, distance_center_bp))
      
      
      ### Generate footprint plot
      message("Generate individual TF plots...")
      
      footprint_list[[i]] =
        ggplot(reshaped_score_mean_list[[i]],
               aes(x = distance_center_bp,
                   y = score_mean,
                   color = group,
                   fill = group)) +
        geom_ribbon(data = reshaped_score_mean_list[[i]],
                    aes(x = distance_center_bp,
                        ymin = score_mean-SEM,
                        ymax = score_mean+SEM),
                    alpha = error.transparency,
                    color = NA) +
        geom_line() +
        scale_color_manual(values = group.colors) +
        scale_fill_manual(values = group.colors) +
        ggtitle(TF_list[i]) +
        xlab("Distance from motif center [bp]") +
        ylab("mean Footprint score \u00B1 SEM") +
        theme_classic() +
        theme(axis.text = element_text(color = "black"),
              axis.ticks = element_line(color = "black"),
              plot.title = element_text(hjust = 0.5))
    }
    names(reshaped_score_list) = TF_list
    names(reshaped_score_mean_list) = TF_list
    names(footprint_list) = TF_list
    
    
    # ----------------------------------------------------------------------------------------------------
    ### Compute area under the curve (AUC)
    message("Computing AUC...")
    
    AUC_bySample = list()
    AUC_byGroup = list()
    
    for (i in 1:length(TF_list)) {
      # table by sample
      AUC_bySample[[i]] =
        reshaped_score_list[[i]] %>%
        dplyr::group_by(sample.ID, group, TF) %>%
        dplyr::summarise(AUC = DescTools::AUC(x = distance_center_bp,
                                              y = score,
                                              na.rm = T),
                         .groups = "keep") %>%
        dplyr::mutate(group = factor(group, levels = levels)) %>%
        dplyr::arrange(group)
      
      
      # Heatmap table (bySample)
      heatmap_mat_bySample = matrix(AUC_bySample[[i]]$AUC, nrow = 1)
      colnames(heatmap_mat_bySample) = AUC_bySample[[i]]$sample.ID
      rownames(heatmap_mat_bySample) = unique(AUC_bySample[[i]]$TF)
      
      if (i == 1) {
        heatmap_matrix_bySample = heatmap_mat_bySample
      } else {
        heatmap_matrix_bySample = rbind(heatmap_matrix_bySample, heatmap_mat_bySample)
      }
      
      
      
      # Table by Group
      AUC_byGroup[[i]] =
        AUC_bySample[[i]] %>%
        dplyr::group_by(group, TF) %>%
        dplyr::summarise(n = n(),
                         AUC_mean = mean(AUC, na.rm = T),
                         AUC.SD = sd(AUC, na.rm = T),
                         AUC.SEM = sd(AUC, na.rm = T)/sqrt(n()),
                         .groups = "keep") %>%
        dplyr::mutate(group = factor(group, levels = levels)) %>%
        dplyr::arrange(group)
      
      
      # Heatmap table (byGroup)
      heatmap_mat_byGroup = matrix(AUC_byGroup[[i]]$AUC_mean, nrow = 1)
      colnames(heatmap_mat_byGroup) = as.vector(AUC_byGroup[[i]]$group)
      rownames(heatmap_mat_byGroup) = unique(AUC_byGroup[[i]]$TF)
      
      if (i == 1) {
        heatmap_matrix_byGroup = heatmap_mat_byGroup
      } else {
        heatmap_matrix_byGroup = rbind(heatmap_matrix_byGroup, heatmap_mat_byGroup)
      }
    }
    names(AUC_byGroup) = TF_list
    names(AUC_bySample) = TF_list
    
    
    # Computing the Fold_Changes + stats of the required comparisons
    message("Compute Fold Changes...")
    
    FC_table_list = list()
    FC_MAplot_list = list()
    Diff_MAplot_list = list()
    
    for (j in 1:length(comparison.list)) {
      FC_table_list[[j]] = data.frame()
      
      for (i in 1:length(TF_list)) {
        if (ignore.statistics == F) {
        wilcox_test = wilcox.test(x = dplyr::filter(AUC_bySample[[i]], group == comparison.list[[j]][1])$AUC,
                                  y = dplyr::filter(AUC_bySample[[i]], group == comparison.list[[j]][2])$AUC,
                                  paired = paired.stats,
                                  conf.level = confidence.level,
                                  conf.int = compute.confidence.interval)
        } else {
          wilcox_test = list(p.value = 0, estimate = 0)
        }
        
        FC_table_list[[j]] = rbind(FC_table_list[[j]],
                                   data.frame(TF = TF_list[[i]],
                                              group.A = comparison.list[[j]][1],
                                              group.B = comparison.list[[j]][2],
                                              basemean.AUC = mean(dplyr::filter(AUC_bySample[[i]], group %in% comparison.list[[j]])$AUC, na.rm=T),
                                              mean.AUC.A = mean(dplyr::filter(AUC_bySample[[i]], group == comparison.list[[j]][1])$AUC, na.rm=T),
                                              mean.AUC.B = mean(dplyr::filter(AUC_bySample[[i]], group == comparison.list[[j]][2])$AUC, na.rm=T)) %>%
                                     dplyr::mutate(`difference.mean.AUC.A-B` = mean.AUC.A - mean.AUC.B,
                                                   Foldchange_A.vs.B = abs(mean.AUC.A / mean.AUC.B),
                                                   p.value = wilcox_test$p.value,
                                                   p.adjusted = wilcox_test$p.value,
                                                   difference.in.location = as.vector(wilcox_test$estimate)))
      }
      
      # Adjust p.values
      FC_table_list[[j]] =
        dplyr::mutate(FC_table_list[[j]],
                      p.adjusted = p.adjust(p.value,
                                            method = p.adjust.method,
                                            n = length(p.value))) %>%
        dplyr::mutate(activity = ifelse((Foldchange_A.vs.B >= foldChange.threshold) & (abs(basemean.AUC) > basemean.AUC.threshold),
                                        yes = "differential", no = "unresponsive")) %>%
        dplyr::mutate(activity = factor(activity, levels = c("differential", "unresponsive")))
      
      names(FC_table_list)[j] = paste0(comparison.list[[j]][1], ".vs.", comparison.list[[j]][2])
      
      
      #*** plot MA-plots on the foldchanges
      FC_MAplot_list[[j]] =
        ggplot(data = FC_table_list[[j]],
               aes(x = basemean.AUC,
                   y = log2(Foldchange_A.vs.B),
                   color = activity)) +
        geom_point(size = 2, alpha = 0.5) +
        scale_color_manual(values = c("differential" = "indianred",
                                      "unresponsive" = "gray50"),
                           drop = F) +
        ggtitle(paste0(comparison.list[[j]][1], " vs ", comparison.list[[j]][2])) +
        xlab(paste0("mean AUC (", comparison.list[[j]][1], ":", comparison.list[[j]][2], ")")) +
        ylab(paste0("log2(AUC Foldchange): ", comparison.list[[j]][1], " / ", comparison.list[[j]][2])) +
        geom_hline(yintercept = log2(1), linetype = "dotted", color = "gray50") +
        geom_hline(yintercept = c(-1,1)*log2(foldChange.threshold), linetype = "dashed", color = "steelblue") +
        geom_vline(xintercept = basemean.AUC.threshold, linetype = "dashed", color = "steelblue") +
        theme_classic() +
        theme(axis.text = element_text(color = "black"),
              axis.ticks = element_line(color = "black"),
              plot.title = element_text(hjust = 0.5))
      
      # Add significant TF.names if required
      if (MAplot.show.names == TRUE) {
        FC_MAplot_list[[j]] =
          FC_MAplot_list[[j]] + 
          ggrepel::geom_text_repel(data = FC_table_list[[j]] %>% dplyr::filter(activity == "differential"),
                                   aes(x = basemean.AUC,
                                       y = log2(Foldchange_A.vs.B),
                                       color = activity,
                                       label = TF),
                                   max.overlaps = 100,
                                   force = 5,
                                   size = 3,
                                   show.legend = F)
      }
      
      # Make y-axis symmetric
      ggbuild = ggplot_build(FC_MAplot_list[[j]])
      y.max = max(abs(ggbuild$layout$panel_params[[1]]$y.range))
      FC_MAplot_list[[j]] = FC_MAplot_list[[j]] + ylim(c(-y.max, y.max))
      
      names(FC_MAplot_list)[j] = paste0(comparison.list[[j]][1], ".vs.", comparison.list[[j]][2])
      
      
      
      #*** plot MA-plots on the _differences_
      # get quantiles fo threshold
      quantiles = quantile(log2(FC_table_list[[j]]$`difference.mean.AUC.A-B`+1),
                           c(percentage.top.hits/100, 1-(percentage.top.hits/100)),
                           na.rm = T)
      
      Diff_MAplot_list[[j]] =
        ggplot(data = FC_table_list[[j]] %>% dplyr::filter(log2(`difference.mean.AUC.A-B`+1)>quantiles[1] & log2(`difference.mean.AUC.A-B`+1)<quantiles[2]),
               aes(x = basemean.AUC,
                   y = log2(`difference.mean.AUC.A-B`+1))) +
        geom_point(size = 2, alpha = 0.5, show.legend = F, color = "gray50") +
        geom_point(data = FC_table_list[[j]] %>% dplyr::filter((log2(`difference.mean.AUC.A-B`+1)<=quantiles[1] | log2(`difference.mean.AUC.A-B`+1)>=quantiles[2]) & (abs(basemean.AUC) > basemean.AUC.threshold)),
                   size = 2, alpha = 0.5, show.legend = F, color = "#07c22c") +
        ggtitle(paste0(comparison.list[[j]][1], " vs ", comparison.list[[j]][2], " (top ", percentage.top.hits, "%)")) +
        xlab(paste0("mean AUC (", comparison.list[[j]][1], ":", comparison.list[[j]][2], ")")) +
        ylab(paste0("log2(AUC difference + 1): ", comparison.list[[j]][1], " - ", comparison.list[[j]][2])) +
        geom_hline(yintercept = log2(1), linetype = "dotted", color = "gray50") +
        geom_hline(yintercept = quantiles, linetype = "dashed", color = "steelblue") +
        geom_vline(xintercept = basemean.AUC.threshold, linetype = "dashed", color = "steelblue") +
        theme_classic() +
        theme(axis.text = element_text(color = "black"),
              axis.ticks = element_line(color = "black"),
              plot.title = element_text(hjust = 0.5))
      
      # Add significant TF.names if required
      if (MAplot.show.names == TRUE) {
        Diff_MAplot_list[[j]] =
          Diff_MAplot_list[[j]] + 
          ggrepel::geom_text_repel(data = FC_table_list[[j]] %>% dplyr::filter((log2(`difference.mean.AUC.A-B`+1)<=quantiles[1] | log2(`difference.mean.AUC.A-B`+1)>=quantiles[2]) & (abs(basemean.AUC) > basemean.AUC.threshold)),
                                   aes(x = basemean.AUC,
                                       y = log2(`difference.mean.AUC.A-B`+1),
                                       label = TF),
                                   color = "#07c22c",
                                   max.overlaps = 100,
                                   force = 5,
                                   size = 3,
                                   show.legend = F)
      }
      
      # Make y-axis symmetric
      ggbuild = ggplot_build(Diff_MAplot_list[[j]])
      y.max = max(abs(ggbuild$layout$panel_params[[1]]$y.range))
      Diff_MAplot_list[[j]] = Diff_MAplot_list[[j]] + ylim(c(-y.max, y.max))
      
      names(Diff_MAplot_list)[j] = paste0(comparison.list[[j]][1], ".vs.", comparison.list[[j]][2])
      
    }
    
    
    
    
    # ----------------------------------------------------------------------------------------------------
    ### Plot heatmap for all conditions (by Sample)
    palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = heatmap.palette))(heatmap.nColors)
    
    # build sample annotations (bySample)
    col_anno_bySample = data.frame(group = sample.groups$group)
    rownames(col_anno_bySample) = sample.groups$sample.ID
    
    # build sample annotations (byGroup)
    col_anno_byGroup = data.frame(group = levels)
    rownames(col_anno_byGroup) = levels
    
    # build sample color annotations
    col_anno_colors = list(group.colors)
    names(col_anno_colors[[1]]) = levels
    names(col_anno_colors) = "group"
    
    
    # Remove NaN from matrices
    heatmap_matrix_bySample[is.nan(heatmap_matrix_bySample)] = NA
    heatmap_matrix_byGroup[is.nan(heatmap_matrix_byGroup)] = NA
    
    
    # Plot linear and z-score heatmap (bySample)
    heatmap_linear_bySample %<a-%
      pheatmap::pheatmap(mat = heatmap_matrix_bySample[!is.nan(rowMeans(heatmap_matrix_bySample, na.rm = T)),],
                         cluster_rows = heatmap.cluster.rows,
                         cluster_cols = heatmap.cluster.columns,
                         show_rownames = heatmap.showRowNames,
                         show_colnames = heatmap.showColNames,
                         scale = "none",
                         kmeans_k = heatmap.n.row.kmeans,
                         border_color = heatmap.border.color,
                         annotation_col = col_anno_bySample,
                         annotation_colors = col_anno_colors,
                         color = palette,
                         na_col = palette[1])
    
    heatmap_scaled_bySample %<a-%
      pheatmap::pheatmap(mat = heatmap_matrix_bySample[!is.nan(rowMeans(heatmap_matrix_bySample, na.rm = T)),],
                         cluster_rows = heatmap.cluster.rows,
                         cluster_cols = heatmap.cluster.columns,
                         show_rownames = heatmap.showRowNames,
                         show_colnames = heatmap.showColNames,
                         scale = heatmap.scaling.mode,
                         kmeans_k = heatmap.n.row.kmeans,
                         border_color = heatmap.border.color,
                         annotation_col = col_anno_bySample,
                         annotation_colors = col_anno_colors,
                         color = palette,
                         na_col = palette[1])
    
    
    # Plot linear and z-score heatmap (byGroup)
    heatmap_linear_byGroup %<a-%
      pheatmap::pheatmap(mat = heatmap_matrix_byGroup[!is.nan(rowMeans(heatmap_matrix_byGroup, na.rm = T)),],
                         cluster_rows = heatmap.cluster.rows,
                         cluster_cols = heatmap.cluster.columns,
                         show_rownames = heatmap.showRowNames,
                         show_colnames = heatmap.showColNames,
                         scale = "none",
                         kmeans_k = heatmap.n.row.kmeans,
                         border_color = heatmap.border.color,
                         annotation_col = col_anno_byGroup,
                         annotation_colors = col_anno_colors,
                         color = palette,
                         na_col = palette[1])
    
    heatmap_scaled_byGroup %<a-%
      pheatmap::pheatmap(mat = heatmap_matrix_byGroup[!is.nan(rowMeans(heatmap_matrix_byGroup, na.rm = T)),],
                         cluster_rows = heatmap.cluster.rows,
                         cluster_cols = heatmap.cluster.columns,
                         show_rownames = heatmap.showRowNames,
                         show_colnames = heatmap.showColNames,
                         scale = heatmap.scaling.mode,
                         kmeans_k = heatmap.n.row.kmeans,
                         border_color = heatmap.border.color,
                         annotation_col = col_anno_byGroup,
                         annotation_colors = col_anno_colors,
                         color = palette,
                         na_col = palette[1])
    
    
    # ----------------------------------------------------------------------------------------------------
    ### Building the output list
    output_list =
      list(metadata = list(sample.configuration = sample.groups,
                           comparisons = comparison.list,
                           TF.analyzed = TF_list,
                           group.colors = col_anno_colors[[1]]),
           source.data = list(source.tables = score_list,
                              ggplot.format.source.tables = reshaped_score_list),
           footprint.analyses = list(ggplot.format.group.means = reshaped_score_mean_list,
                                     footprint.plots = footprint_list),
           AUC.analyses = list(AUC.table.bySample = AUC_bySample,
                               AUC.table.byGroup = AUC_byGroup,
                               AUC.foldchange.comparisons = FC_table_list,
                               AUC.comparisons.FC.MAplot = FC_MAplot_list,
                               AUC.comparisons.delta.MAplot = Diff_MAplot_list,
                               heatmaps = list(matrix.AUC.bySample = heatmap_matrix_bySample,
                                               matrix.AUC.byGroup = heatmap_matrix_byGroup,
                                               heatmap.AUC.bySample.raw = heatmap_linear_bySample,
                                               heatmap.AUC.byGroup.raw = heatmap_linear_byGroup,
                                               heatmap.AUC.bySample.zScore = heatmap_scaled_bySample,
                                               heatmap.AUC.byGroup.zScore = heatmap_scaled_byGroup,
                                               hetamap.annotations = list(sample.annotations = col_anno_bySample,
                                                                          group.annotations = col_anno_byGroup,
                                                                          group.annotation.colors = col_anno_colors))))
    
    return(output_list)
  } # END function
