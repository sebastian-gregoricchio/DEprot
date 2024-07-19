#' @title protein.summary
#'
#' @description Plots stacked-barplot recapitulating
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param group.column String indicating a column among the ones in the metadata table. Default: \code{"column.id"} (no groups).
#' @param n.lables String indicating the type of values to display on the barplot. One among: \code{NULL} (no labels), "frequency", "percentage", "counts". Default: \code{NULL} (no labels).
#' @param label.color String indicating the font-color to use for the barplot values labels. Default: \code{"white"}.
#' @param x.label.angle Numeric value indicating the rotation angle to use for the x-axis labels. Default: \code{30} (degrees).
#' @param show.frequency Logical value indicating whether the barplot should show the y-axis as frequency or as absolute counts. Default: \code{FALSE} (counts).
#' @param colors String-vector indicating the colors to use for the degrees of protein-presence. Default: \code{NULL} (colors will be assigned automatically).
#' @param title String indicating the title to use (markdown annotation supported). Default: \code{NULL}.
#' @param subtitle String indicating the subtitle to use (markdown annotation supported). Default: \code{NULL}.
#'
#' @return A ggplot object.
#'
#' @export protein.summary

protein.summary =
  function(DEprot.object,
           group.column = "column.id",
           n.lables = NULL, # frequency, percentage, counts
           label.color = "white",
           x.label.angle = 30,
           show.frequency = FALSE,
           colors = NULL,
           title = NULL,
           subtitle = NULL) {

    ### Libraries
    require(dplyr)
    require(ggplot2)


    ### check object and extract metadata table
    if (!("DEprot" %in% class(DEprot.object))) {
      if (!("DEprot.analyses" %in% class(DEprot.object))) {
        warning("The input must be an object of class 'DEprot' or 'DEprot.analyses'.")
        return(DEprot.object)
      }
    }



    ### Check and extract counts table
    if (!is.null(DEprot.object@raw.counts)) {
      mat = DEprot.object@raw.counts
    } else if (!is.null(DEprot.object@norm.counts)) {
      mat = DEprot.object@norm.counts
    } else {
      warning("Your object must contain either 'raw' or 'normalized' data.")
      return(DEprot.object)
    }


    ### check group column
    if (!(group.column %in% colnames(DEprot.object@metadata))) {
      warning(paste0("The group column '", group.column, "' is not present in the metadata table.\n",
                     "Available columns: ", paste0(colnames(DEprot.object@metadata), collapse = ", ")))
      return(DEprot.object)
    } else {
      meta = DEprot.object@metadata
    }


    ### Convert NaN to NA
    mat[is.nan(mat)] = NA
    mat[is.na(mat)] = 0

    ### Convert matrix to TRUE/FALSE (tf)
    mat.tf = mat > 0



    ### Summarize per group
    groups = unique(meta[,group.column])

    presence.table = data.frame()
    for (i in 1:length(groups)) {
      samples.in.group = meta[meta[,group.column] == groups[i], "column.id"]

      presence.table = rbind(presence.table,
                             data.frame(group = groups[i],
                                        shared = rowSums(data.frame(mat.tf[,samples.in.group]), na.rm = T)))
    }



    ### Count observations
    presence.recap =
      presence.table %>%
      dplyr::group_by(group, shared) %>%
      dplyr::summarise(.groups = "keep", n = n()) %>%
      dplyr::mutate(shared = ifelse(is.na(shared), yes = 0, no = shared),
                    n.total = nrow(mat)) %>%
      dplyr::mutate(fraction = n / n.total) %>%
      dplyr::mutate(percentage = fraction*100) %>%
      dplyr::select(-n.total)



    ### Make base plot
    # define colors
    n.categories = length(unique(presence.recap$shared))

    if (is.null(colors)) {
      shared.color = viridis::viridis(n = n.categories, direction = -1, end = 0.8)
    } else {
      if (length(colors) < n.categories) {
        message(paste0("The lenght of the 'colors' vector [n=", length(colors),
                       "] is lower than the number of colors required [n=", n.categories,"].\n",
                       "Colors will be re-defined automatically using the viridis palette."))
        shared.color = viridis::viridis(n = n.categories, direction = -1, end = 0.8)
      } else {
        shared.color = colors[1:n.categories]
      }
    }

    names(shared.color) = as.character(sort(unique(presence.recap$shared)))


    ### Barplot
    if (show.frequency == T) {
      barplot =
        ggplot(data = presence.recap,
               aes(x = group,
                   y = fraction,
                   fill = factor(as.character(shared), levels = as.character(names(shared.color))))) +
        geom_bar(stat = "identity", position = "fill")
    } else {
      barplot =
        ggplot(data = presence.recap,
               aes(x = group,
                   y = n,
                   fill = factor(as.character(shared), levels = as.character(names(shared.color))))) +
        geom_bar(stat = "identity")
    }

    barplot =
      barplot +
      ylab(ifelse(test = show.frequency == T, yes = "Protein frequency", no = "Protein count")) +
      xlab(NULL) +
      ggtitle(label = title, subtitle = subtitle) +
      scale_y_continuous(expand = c(0,0)) +
      theme_classic() +
      scale_fill_manual(values = shared.color,
                        name = "Detected in<br>*N* samples",
                        breaks = names(shared.color)) +
      theme(axis.line.x = element_blank(),
            axis.text.y = element_text(colour = "black"),
            axis.text.x = element_text(colour = "black", angle = x.label.angle, hjust = 1),
            axis.ticks.x = element_blank(),
            legend.title = ggtext::element_markdown(color = "black"),
            legend.text = element_text(colour = "black"),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(color = "black", hjust = 0.5))


    ### Add labels with numbers
    if (!is.null(n.lables)) {
      if (tolower(n.lables) %in% c("freq", "frequency", "frequence", "fraction", "frac")) {
        barplot =
          barplot +
          geom_text(aes(label = round(fraction, digits = 2)),
                    color = label.color,
                    show.legend = F,
                    position = position_stack(vjust = 0.5))
      } else if (tolower(n.lables) %in% c("percentage", "perc", "%")) {
        barplot =
          barplot +
          geom_text(aes(label = round(percentage, digits = 1)),
                    color = label.color,
                    show.legend = F,
                    position = position_stack(vjust = 0.5))
      } else {
        barplot =
          barplot +
          geom_text(aes(label = n),
                    color = label.color,
                    show.legend = F,
                    position = position_stack(vjust = 0.5))
      }
    }


    ### Export plot
    return(barplot)
  } # END of function
