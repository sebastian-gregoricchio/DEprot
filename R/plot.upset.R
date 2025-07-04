#' @title plot.upset
#'
#' @description Plots an upset-plot for differential expression results
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast.subset Numeric vector indicating the contrasts to use. Default: \code{NULL} (all contrasts are shown).
#' @param title String indicating the title to display. Markdown-mode supported. Default: \code{NULL} (no title).
#' @param subtitle String indicating the title to display. Markdown-mode supported. Default: \code{NULL} (subtitle).
#' @param sort.intersections String indicating which method to use for the intersection sorting: 'cardinality' or 'degree'. Default: \code{"cardinality"}.
#' @param sort.sets String indicating which method to use for the sets sorting: 'descending' or 'ascending'. Default: \code{"descending"}.
#' @param intersection.bar.color String indicating the color to use for the intersection size bar plot. Default: \code{"black"}.
#' @param setsize.bar.color String indicating the color to use for the set size bar plot. Default: \code{"black"}.
#' @param show.counts Logical value indicating whether the counts in the intersection bar pot should be shown. Default: \code{TRUE}.
#' @param min.size Numeric value indicating the minimal number of interactions to show. Default: \code{1}.
#' @param height.ratio Numeric value indicating the ratio of the intersection matrix to intersection size height. Default: \code{0.5}.
#' @param width.ratio Numeric value indicating the ratio of the overall set size width to intersection matrix width. Default: \code{0.3}.
#' @param use.uncorrected.pvalue Logical value indicating whether it should be used the normal p-value instead of the adjusted one (differential proteins numbers are recomputed). Default: \code{FALSE}, padj is used.
#'
#' @return A \code{DEprot.upset} with the upset plot and the matrix of the observed enriched proteins.
#'
#' @name plot.upset
#'
#' @export plot.upset

plot.upset =
  function(DEprot.analyses.object,
           contrast.subset = NULL,
           title = NULL,
           subtitle = NULL,
           sort.intersections = "cardinality",  # 'cardinality' / 'degree'
           sort.sets = "descending",
           intersection.bar.color = "black",
           setsize.bar.color = "black",
           show.counts = TRUE,
           min.size = 1,
           height.ratio = 0.5,
           width.ratio = 0.3,
           use.uncorrected.pvalue = FALSE) # 'ascending', 'descending', FALSE
    {

    ### Libraries
    require(dplyr)
    require(ggplot2)
    require(ComplexUpset)


    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      warning("The input must be an object of class 'DEprot.analyses'.")
      return()
    }


    ### Subset contrasts
    if (!is.null(contrast.subset)) {
      if (all(contrast.subset %in% 1:length(DEprot.analyses.object@analyses.result.list))) {
        analyses.result.list = DEprot.analyses.object@analyses.result.list[contrast.subset]
        contrasts = DEprot.analyses.object@contrasts[contrast.subset]
      } else {
        warning("Not all the contrasts indicated in the subset are present in the 'analyses.result.list' of the object provided.")
        return()
      }
    } else {
      analyses.result.list = DEprot.analyses.object@analyses.result.list
      contrasts = DEprot.analyses.object@contrasts
    }


    ### Redefine up and downs if required
    if (use.uncorrected.pvalue == TRUE) {

      de.status =
        function(FC, p, contr) {
          ifelse(p < DEprot.analyses.object@differential.analyses.params$padj.th,
                 yes = ifelse(abs(FC) >= log2(DEprot.analyses.object@differential.analyses.params$linear.FC.th),
                              yes = ifelse(sign(FC) == 1,
                                           yes = contr$var.1,
                                           no = contr$var.2),
                              no = "null"),
                 no = ifelse(FC >= log2(DEprot.analyses.object@differential.analyses.params$linear.FC.unresp.range[1]) & FC <= log2(DEprot.analyses.object@differential.analyses.params$linear.FC.unresp.range[2]),
                             yes = "unresponsive",
                             no = "null"))
        }

      for (i in 1:length(analyses.result.list)) {
        analyses.result.list[[i]]$results =
          analyses.result.list[[i]]$results %>%
          dplyr::mutate(diff.status = de.status(FC = analyses.result.list[[i]]$results[,5],
                                                p = analyses.result.list[[i]]$results[,6],
                                                contr = contrasts[[i]]))
      }
    }




    ### merge results
    overlaps.tb = data.frame(analyses.result.list[[1]]$results[,"prot.id"])
    colnames(overlaps.tb)[1] = "prot.id"

    for (i in 1:length(analyses.result.list)) {
      tb =
        dplyr::mutate(.data = analyses.result.list[[i]]$results,
                      diff.status = paste0(contrasts[[i]]$metadata.column, ": ",
                                           contrasts[[i]]$var.1, " *vs* ",
                                           contrasts[[i]]$var.2, " | **",
                                           diff.status, "**"))

      ### filter not differential
      tb = dplyr::filter(.data = tb, !grepl(x = diff.status, "unresponsive|null"))

      if (nrow(tb) > 0) {
        ## split groups
        groups = unique(tb$diff.status)
        group1 = dplyr::filter(.data = tb, diff.status == groups[1])
        group2 = dplyr::filter(.data = tb, diff.status == groups[2])

        ## check which proteins are in group1
        if (!(groups[1] %in% colnames(overlaps.tb))) {
          if (nrow(group1) > 0) {
            overlaps.tb = dplyr::mutate(.data = overlaps.tb, group1 = prot.id %in% group1$prot.id)
          } else {
            overlaps.tb$group1 = F
          }
          colnames(overlaps.tb)[ncol(overlaps.tb)] = groups[1]
        }


        ## check which proteins are in group2
        if (!(groups[2] %in% colnames(overlaps.tb))) {
          if (nrow(group2) > 0) {
            overlaps.tb = dplyr::mutate(.data = overlaps.tb, group2 = prot.id %in% group2$prot.id)
          } else {
            overlaps.tb$group2 = F
          }
          colnames(overlaps.tb)[ncol(overlaps.tb)] = groups[2]
        }
      }
    }


    if (ncol(overlaps.tb) > 1) {
      # Remove proteins missing everywhere
      overlaps.tb = overlaps.tb[rowSums(overlaps.tb[,-1]) > 0,]
    } else {
      warning("No differential proteins have been found. Upset plot cannot be generated")
      return()
    }


    ## Make upset plot
    if (nrow(overlaps.tb) > 0) {
      upset =
        upset(data = overlaps.tb,
              intersect = colnames(overlaps.tb)[-1],
              name = NULL,
              base_annotations = list("Intersection size" =
                                        intersection_size(fill = intersection.bar.color,
                                                          counts = show.counts,
                                                          width = 0.75) +
                                        scale_y_continuous(expand = c(0,0)) +
                                        ggtitle(label = title, subtitle = subtitle) +
                                        theme(axis.text = ggtext::element_markdown(color = "black"),
                                              panel.grid.major.x = element_blank(),
                                              axis.line.y = element_line(colour = "black"),
                                              axis.ticks.y = element_line(colour = "black"),
                                              plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
                                              plot.subtitle = ggtext::element_markdown(color = "black", hjust = 0.5))),
              set_sizes =
                upset_set_size(geom = geom_bar(fill = setsize.bar.color, width = 0.4, show.legend = F)) +
                scale_y_reverse(expand = c(0,0)) +
                theme(axis.line.x = element_line(colour = "black"),
                      axis.ticks.x = element_line(colour = "black")),
              themes = upset_modify_themes(list("intersections_matrix" = theme(axis.text = ggtext::element_markdown(color = "black"),
                                                                               panel.grid.major.x = element_blank()),
                                                "overall_sizes" = theme(axis.text = ggtext::element_markdown(color = "black")))),
              sort_intersections_by = sort.intersections,
              sort_sets = sort.sets,
              height_ratio = height.ratio,
              width_ratio = width.ratio,
              min_size = min.size)
    } else {
      warning("No differential proteins have been found. Upset plot cannot be generated")
      return()
    }


    ### Export upset object
    DEprot.upset.object =
      new(Class = "DEprot.upset",
          upset = upset,
          obs.matrix = overlaps.tb)

    return(DEprot.upset.object)

  } # END function

