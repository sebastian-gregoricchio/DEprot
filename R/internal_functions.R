##########################
### INTERNAL FUNCTIONS ###
##########################

# ----------------------------------------------------------------------------------------

#' @title corum.reshape
#'
#' @description
#' Reshapes the CORUM database downloads files for gene entrenchment analyses
#'
#' @param corum_table Either the path to the CORUM database file or a data.frame with the same structure.
#'
#' @return Data.frame with 4 columns: complex.id, complex.name, organism, protein.members.
#'
#' @importFrom data.table fread
#' @importFrom stringr str_split
#'
#' @keywords internal

.corum.reshape =
  function(corum_table) {

    ## import annotations
    if ("data.frame" %in% class(corum_table)) {
      corum = as.data.frame(corum_table) }
    else {
      corum = data.table::fread(corum_table, data.table = FALSE)
    }

    ## reshape database
    members = stringr::str_split(corum$subunits_gene_name, pattern = ";")
    n = lengths(members)

    corum_list = data.frame(complex.id = rep(corum$complex_id, n),
                            complex.name = rep(corum$complex_name, n),
                            organism = rep(corum$organism, n),
                            protein.members = unlist(members))

    return(corum_list)

  } # END corum.reshape


# ----------------------------------------------------------------------------------------


#' @title .ordination.slots
#'
#' @description Internal. Harmonizes the access to a \code{DEprot.PCA} or a \code{DEprot.PCoA} object,
#' returning the elements needed by the shared plotting helpers under a common set of names.
#'
#' @param object An object of class \code{DEprot.PCA} or \code{DEprot.PCoA}.
#'
#' @return A list with the elements: \code{coordinates}, \code{metadata}, \code{importance},
#' \code{axis.prefix}, \code{axis.title}, \code{positive.axis} (\code{NULL} for a PCA, a logical
#' vector for a PCoA) and \code{object.class}.
#'
#' @keywords internal

.ordination.slots =
  function(object) {

    if ("DEprot.PCA" %in% class(object)) {
      return(list(coordinates = object@PCs,
                  metadata = object@PCA.metadata,
                  importance = object@importance,
                  axis.prefix = "PC",
                  axis.title = "Principal Component (PC)",
                  positive.axis = NULL,
                  object.class = "DEprot.PCA"))

    } else if ("DEprot.PCoA" %in% class(object)) {
      return(list(coordinates = object@PCos,
                  metadata = object@PCoA.metadata,
                  importance = object@importance,
                  axis.prefix = "PCo",
                  axis.title = "Principal Coordinate (PCo)",
                  positive.axis = object@importance$Positive.eigenvalue,
                  object.class = "DEprot.PCoA"))

    } else {
      stop("The input must be an object of class 'DEprot.PCA' or 'DEprot.PCoA'.")
    }
  }


# ----------------------------------------------------------------------------------------


#' @title .check.aes.columns
#'
#' @description Internal. Validates the columns requested for the color/shape/label aesthetics of the
#' ordination plots and prepares the corresponding vectors and legend labels.
#'
#' @param coordinates data.frame combining the ordination coordinates and the metadata
#' (i.e. the \code{PCs} slot of a \code{DEprot.PCA} or the \code{PCos} slot of a \code{DEprot.PCoA}).
#' @param color.column String indicating the name of the column to use as factor for the dot colors.
#' @param shape.column String indicating the name of the column to use as factor for the dot shapes, or \code{NULL}.
#' @param label.column String indicating the name of the column to use as label of the dots, or \code{NULL}.
#' @param axis.prefix String indicating the prefix of the coordinate columns ('PC' or 'PCo').
#'
#' @return A list with the elements: \code{shape.scores}, \code{shape.label}, \code{label.column}, \code{show.labels}.
#'
#' @import ggplot2
#' @importFrom stringr str_to_title
#'
#' @keywords internal

.check.aes.columns =
  function(coordinates,
           color.column,
           shape.column = NULL,
           label.column = NULL,
           axis.prefix = "PC") {

    ## columns available for the aesthetics (i.e. all but the coordinates)
    axis.pattern = paste0("^", axis.prefix, "[0-9]+$")
    available.columns = colnames(coordinates)[!grepl(axis.pattern, colnames(coordinates))]

    .aes.error = function(what, value) {
      stop(paste0("The ", what, " column '", value, "' is not present in the ",
                  axis.prefix, " analyses table.\n",
                  "       Available columns: ", paste0(available.columns, collapse = ", ")))
    }


    ### color
    if (!(color.column %in% colnames(coordinates))) {.aes.error("color", color.column)}


    ### shape
    if (!is.null(shape.column)) {
      if (!(shape.column %in% colnames(coordinates))) {.aes.error("shape", shape.column)}
      shape.scores = coordinates[,shape.column]
      shape.label = guide_legend(title = stringr::str_to_title(shape.column))
    } else {
      shape.scores = "Samples"
      shape.label = "none"
    }


    ### label
    if (!is.null(label.column)) {
      if (!(label.column %in% colnames(coordinates))) {.aes.error("label", label.column)}
      show.labels = TRUE
    } else {
      label.column = color.column
      show.labels = FALSE
    }


    return(list(shape.scores = shape.scores,
                shape.label = shape.label,
                label.column = label.column,
                show.labels = show.labels))
  }


# ----------------------------------------------------------------------------------------


#' @title .check.ordination.axes
#'
#' @description Internal. Verifies that the two axes requested exist in the ordination and, for a PCoA,
#' that they do not correspond to a negative eigenvalue (which has no real representation).
#'
#' @param coordinates data.frame combining the ordination coordinates and the metadata.
#' @param axis.x Number of the axis displayed on the x-axis.
#' @param axis.y Number of the axis displayed on the y-axis.
#' @param axis.prefix String indicating the prefix of the coordinate columns ('PC' or 'PCo').
#' @param positive.axis Logical vector indicating which axes correspond to a positive eigenvalue, or \code{NULL} (PCA).
#'
#' @return Invisibly \code{TRUE}; it stops with an informative error otherwise.
#'
#' @keywords internal

.check.ordination.axes =
  function(coordinates,
           axis.x,
           axis.y,
           axis.prefix = "PC",
           positive.axis = NULL) {

    axis.pattern = paste0("^", axis.prefix, "[0-9]+$")
    available.axes = grep(axis.pattern, colnames(coordinates), value = TRUE)

    if (!(paste0(axis.prefix, axis.x) %in% available.axes) |
        !(paste0(axis.prefix, axis.y) %in% available.axes)) {
      stop(paste0("The axes requested (", axis.prefix, axis.x, ", ", axis.prefix, axis.y,
                  ") are not available.\n",
                  "       Available axes: ", paste0(available.axes, collapse = ", ")))
    }

    if (!is.null(positive.axis)) {
      if (FALSE %in% positive.axis[c(axis.x, axis.y)]) {
        stop(paste0("At least one of the axes requested (", axis.prefix, axis.x, ", ",
                    axis.prefix, axis.y, ") corresponds to a negative eigenvalue and cannot be plotted.\n",
                    "       Re-run 'perform.PCoA' using 'distance.transformation = \"sqrt\"' or a 'correction'."))
      }
    }

    return(invisible(TRUE))
  }


# ----------------------------------------------------------------------------------------


#' @title .ordination.scatter
#'
#' @description Internal. Builds the scatter plot of two ordination axes. This is the shared body of
#' \link{plot.PC.scatter} and \link{plot.PCoA.scatter}.
#'
#' @param coordinates data.frame combining the ordination coordinates and the metadata.
#' @param importance importance/summary table of the ordination.
#' @param axis.prefix String indicating the prefix of the coordinate columns ('PC' or 'PCo').
#' @param axis.x Number of the axis to display on the x-axis.
#' @param axis.y Number of the axis to display on the y-axis.
#' @param color.column String indicating the name of the column to use as factor for the dot colors.
#' @param shape.column String indicating the name of the column to use as factor for the dot shapes.
#' @param label.column String indicating the name of the column to use as label of the dots.
#' @param plot.zero.line.x Logical value indicating whether to plot a gray dashed line at x=0.
#' @param plot.zero.line.y Logical value indicating whether to plot a gray dashed line at y=0.
#' @param positive.axis Logical vector indicating which axes correspond to a positive eigenvalue, or \code{NULL}.
#' @param point.size Numeric value indicating the size of the dots.
#' @param extra.layers List of additional ggplot layers to be inserted between the zero-lines and the points (e.g. the loading arrows of a biplot). Default: \code{NULL}.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom stringr str_to_title
#' @import ggrepel
#'
#' @keywords internal

.ordination.scatter =
  function(coordinates,
           importance,
           axis.prefix,
           axis.x,
           axis.y,
           color.column = "column.id",
           shape.column = NULL,
           label.column = NULL,
           plot.zero.line.x = TRUE,
           plot.zero.line.y = TRUE,
           positive.axis = NULL,
           point.size = 3,
           extra.layers = NULL) {

    ### checks
    aes.list = .check.aes.columns(coordinates = coordinates,
                                  color.column = color.column,
                                  shape.column = shape.column,
                                  label.column = label.column,
                                  axis.prefix = axis.prefix)

    axis.x = round(axis.x, 0)
    axis.y = round(axis.y, 0)

    .check.ordination.axes(coordinates = coordinates,
                           axis.x = axis.x,
                           axis.y = axis.y,
                           axis.prefix = axis.prefix,
                           positive.axis = positive.axis)


    ### Build table for plot
    tb = data.frame(axis.x = coordinates[,paste0(axis.prefix, axis.x)],
                    axis.y = coordinates[,paste0(axis.prefix, axis.y)],
                    color = factor(coordinates[,color.column]),
                    shape = factor(aes.list$shape.scores),
                    label = coordinates[,aes.list$label.column])


    ### Make basic plot
    ordination.plot =
      ggplot(data = tb,
             aes(x = axis.x,
                 y = axis.y,
                 shape = shape,
                 color = color,
                 label = label))

    if (plot.zero.line.x == TRUE) {
      ordination.plot = ordination.plot + geom_vline(xintercept = 0, color = "gray", linetype = 2)
    }

    if (plot.zero.line.y == TRUE) {
      ordination.plot = ordination.plot + geom_hline(yintercept = 0, color = "gray", linetype = 2)
    }


    ## Extra layers (e.g. biplot arrows), added below the points
    if (!is.null(extra.layers)) {
      for (layer in extra.layers) {
        ordination.plot = ordination.plot + layer
      }
    }


    ordination.plot =
      ordination.plot +
      geom_point(size = point.size) +
      theme_classic() +
      xlab(paste0(axis.prefix, axis.x, " (", round(importance$Percentage.of.Variance[axis.x],1), "%)")) +
      ylab(paste0(axis.prefix, axis.y, " (", round(importance$Percentage.of.Variance[axis.y],1), "%)")) +
      theme(aspect.ratio = 1,
            axis.line = element_blank(),
            axis.ticks = element_line(color = "black"),
            axis.text = element_text(color = "black"),
            panel.border = element_rect(color = "black", fill = NA))


    ## Add labels
    if (aes.list$show.labels == TRUE) {
      ordination.plot = ordination.plot + ggrepel::geom_text_repel(max.overlaps = 100,
                                                                   show.legend = FALSE)
    }

    ordination.plot =
      ordination.plot +
      guides(color = guide_legend(title = ifelse(color.column == "column.id",
                                                 yes = "Samples",
                                                 no = stringr::str_to_title(color.column))),
             shape = aes.list$shape.label)


    return(ordination.plot)
  }


# ----------------------------------------------------------------------------------------


#' @title .ordination.scatter.123
#'
#' @description Internal. Builds the combined axis1-vs-axis2 / axis3-vs-axis2 scatters. This is the shared
#' body of \link{plot.PC.scatter.123} and \link{plot.PCoA.scatter.123}.
#'
#' @param coordinates data.frame combining the ordination coordinates and the metadata.
#' @param metadata metadata table of the samples used in the ordination.
#' @param importance importance/summary table of the ordination.
#' @param axis.prefix String indicating the prefix of the coordinate columns ('PC' or 'PCo').
#' @param color.column String indicating the name of the column to use as factor for the dot colors.
#' @param shape.column String indicating the name of the column to use as factor for the dot shapes.
#' @param label.column String indicating the name of the column to use as label of the dots.
#' @param dot.colors Color-vector to use for the points, or \code{NULL} for automatic colors.
#' @param title String indicating the title of the combined plot (markdown supported), or \code{NULL}.
#' @param plot.zero.line.y.12 Logical value for the y=0 line of the axis1-vs-axis2 plot.
#' @param plot.zero.line.x.12 Logical value for the x=0 line of the axis1-vs-axis2 plot.
#' @param plot.zero.line.y.23 Logical value for the y=0 line of the axis3-vs-axis2 plot.
#' @param plot.zero.line.x.23 Logical value for the x=0 line of the axis3-vs-axis2 plot.
#' @param positive.axis Logical vector indicating which axes correspond to a positive eigenvalue, or \code{NULL}.
#'
#' @return A patchwork object.
#'
#' @import ggplot2
#' @import ggtext
#' @import patchwork
#' @importFrom grDevices rainbow
#'
#' @keywords internal

.ordination.scatter.123 =
  function(coordinates,
           metadata,
           importance,
           axis.prefix,
           color.column = "column.id",
           shape.column = NULL,
           label.column = NULL,
           dot.colors = NULL,
           title = NULL,
           plot.zero.line.y.12 = TRUE,
           plot.zero.line.x.12 = TRUE,
           plot.zero.line.y.23 = TRUE,
           plot.zero.line.x.23 = TRUE,
           positive.axis = NULL) {

    ### checks (performed here as well to fail before generating any plot)
    invisible(.check.aes.columns(coordinates = coordinates,
                                 color.column = color.column,
                                 shape.column = shape.column,
                                 label.column = label.column,
                                 axis.prefix = axis.prefix))

    if (!(paste0(axis.prefix, "3") %in% colnames(coordinates))) {
      stop(paste0("At least 3 axes are required to generate the combined scatter.\n",
                  "       Available axes: ",
                  paste0(grep(paste0("^", axis.prefix, "[0-9]+$"), colnames(coordinates), value = TRUE),
                         collapse = ", ")))
    }


    ### Define/check colors
    if (is.factor(metadata[,color.column])) {
      color.labels = levels(metadata[,color.column])
    } else {
      color.labels = unique(metadata[,color.column])
    }

    if (is.null(dot.colors) | length(dot.colors) < length(color.labels)) {
      if (!is.null(dot.colors)) {
        message(paste0("The number of 'dot.colors' provided is lower the required values (",
                       length(color.labels), ").\nValues will be assigned automatically."))
      }
      scatter.colors = rainbow(n = length(color.labels), s = 0.5)
    } else {
      scatter.colors = dot.colors
    }

    if (is.null(names(scatter.colors))) {
      names(scatter.colors) = color.labels
    }


    ##### Generate plots
    plot.1.2 =
      .ordination.scatter(coordinates = coordinates,
                          importance = importance,
                          axis.prefix = axis.prefix,
                          axis.x = 1,
                          axis.y = 2,
                          color.column = color.column,
                          shape.column = shape.column,
                          label.column = label.column,
                          plot.zero.line.x = FALSE,
                          plot.zero.line.y = FALSE,
                          positive.axis = positive.axis) +
      scale_color_manual(values = scatter.colors) +
      theme(legend.position = "none")

    plot.2.3 =
      .ordination.scatter(coordinates = coordinates,
                          importance = importance,
                          axis.prefix = axis.prefix,
                          axis.x = 3,
                          axis.y = 2,
                          color.column = color.column,
                          shape.column = shape.column,
                          label.column = label.column,
                          plot.zero.line.x = FALSE,
                          plot.zero.line.y = FALSE,
                          positive.axis = positive.axis) +
      ylab(NULL) +
      scale_color_manual(values = scatter.colors)


    ## Add zero-lines
    if (plot.zero.line.y.12 == TRUE) {plot.1.2 = plot.1.2 + geom_hline(yintercept = 0, linetype = 2, color = "gray")}
    if (plot.zero.line.y.23 == TRUE) {plot.2.3 = plot.2.3 + geom_hline(yintercept = 0, linetype = 2, color = "gray")}

    if (plot.zero.line.x.12 == TRUE) {plot.1.2 = plot.1.2 + geom_vline(xintercept = 0, linetype = 2, color = "gray")}
    if (plot.zero.line.x.23 == TRUE) {plot.2.3 = plot.2.3 + geom_vline(xintercept = 0, linetype = 2, color = "gray")}


    ### Export combined plot
    if (is.null(title)) {
      return(patchwork::wrap_plots(plot.1.2, plot.2.3, nrow = 1))
    } else {
      return(patchwork::wrap_plots(plot.1.2, plot.2.3, nrow = 1) +
               patchwork::plot_annotation(title = title,
                                          theme = theme(plot.title = ggtext::element_markdown(hjust = 0.5))))
    }
  }


# ----------------------------------------------------------------------------------------


#' @title .cumulative.plot
#'
#' @description Internal. Builds the barplot of the proportion of variance of each ordination axis with the
#' cumulative curve on top. This is the shared body of \link{plot.PC.cumulative} and \link{plot.PCoA.cumulative}.
#'
#' @param importance importance/summary table of the ordination.
#' @param axis.prefix String indicating the prefix of the axes ('PC' or 'PCo').
#' @param axis.title String indicating the x-axis title.
#' @param bar.color String indicating the color of the bars.
#' @param line.color String indicating the color of the cumulative line and dots.
#' @param title String indicating the title of the plot (markdown supported), or \code{NULL}.
#' @param broken.stick Logical value indicating whether to add the broken-stick expectation as a dashed gray line. Requires a \code{Broken.stick} column in \code{importance}. Default: \code{FALSE}.
#' @param filter.positive Logical value indicating whether to keep only the axes with a positive eigenvalue. Requires a \code{Positive.eigenvalue} column in \code{importance}. Default: \code{FALSE}.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @import ggtext
#'
#' @keywords internal

.cumulative.plot =
  function(importance,
           axis.prefix = "PC",
           axis.title = "Principal Component (PC)",
           bar.color = "steelblue",
           line.color = "navyblue",
           title = NULL,
           broken.stick = FALSE,
           filter.positive = FALSE) {

    tb = importance

    ### keep only the axes with a positive eigenvalue (PCoA)
    if (filter.positive == TRUE & ("Positive.eigenvalue" %in% colnames(tb))) {
      tb = tb[tb$Positive.eigenvalue == TRUE, , drop = FALSE]
    }

    ### harmonize the name of the axis column ('PC' for a PCA, 'PCo' for a PCoA)
    if (!(axis.prefix %in% colnames(tb))) {
      stop(paste0("The importance table does not contain the axis column '", axis.prefix, "'."))
    }
    tb$axis = factor(as.character(tb[,axis.prefix]),
                     levels = as.character(tb[,axis.prefix]))


    ### Generate plot
    cumulative.plot =
      ggplot(data = tb,
             aes(x = axis)) +
      geom_bar(mapping = aes(y = Proportion.of.Variance),
               stat = "identity",
               fill = bar.color) +
      geom_line(mapping = aes(y = Cumulative.Proportion, group = 1),
                color = line.color,
                linetype = 1) +
      geom_point(mapping = aes(y = Cumulative.Proportion, group = 1),
                 stroke = NA,
                 size = 2,
                 color = line.color)


    ### Broken-stick expectation (PCoA)
    subtitle = NULL

    if (broken.stick == TRUE & ("Broken.stick" %in% colnames(tb))) {
      cumulative.plot =
        cumulative.plot +
        geom_line(mapping = aes(y = Broken.stick, group = 1),
                  color = "gray50",
                  linetype = 2)
      subtitle = "*dashed gray line: broken-stick expectation*"
    }


    cumulative.plot =
      cumulative.plot +
      ylab("Proportion of variance") +
      xlab(axis.title) +
      ggtitle(label = title, subtitle = subtitle) +
      theme_classic() +
      scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
      theme(axis.ticks.x = element_blank(),
            axis.text = element_text(color = "black"),
            axis.line.x = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.grid.major.y = element_line(color = "gray", linewidth = 0.1),
            panel.grid.minor.y = element_line(color = "gray", linewidth = 0.1),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(color = "black", hjust = 0.5),
            aspect.ratio = 0.5)

    return(cumulative.plot)
  }


# ----------------------------------------------------------------------------------------


#' @title .scale.loadings
#'
#' @description Internal. Ranks the variables by their euclidean distance from the origin in the selected
#' 2D plane and rescales them so that the arrows fit visually within the range of the sample scores.
#' Shared by \link{plot.PC.biplot} and \link{plot.PCoA.biplot}.
#'
#' @param loadings.df data.frame with at least the columns \code{variable}, \code{loading.x}, \code{loading.y}.
#' @param scores.x Numeric vector of the sample scores displayed on the x-axis.
#' @param scores.y Numeric vector of the sample scores displayed on the y-axis.
#' @param n.loadings Integer indicating the number of top variables to retain.
#' @param loading.scale Numeric multiplier to manually adjust the length of the arrows.
#'
#' @return A list with the elements \code{loadings} (the full, scaled and ranked table),
#' \code{top} (the top \code{n.loadings} rows) and \code{scale.factor}.
#'
#' @keywords internal

.scale.loadings =
  function(loadings.df,
           scores.x,
           scores.y,
           n.loadings = 10,
           loading.scale = 0.8) {

    ## Euclidean distance from the origin, used for the ranking
    loadings.df$distance = sqrt(loadings.df$loading.x^2 + loadings.df$loading.y^2)

    ## Scale the loadings to the range of the scores
    score.range = max(abs(c(scores.x, scores.y)), na.rm = TRUE)
    loading.range = max(abs(c(loadings.df$loading.x, loadings.df$loading.y)), na.rm = TRUE)

    if (loading.range > 0) {
      scale.factor = (score.range / loading.range) * loading.scale
    } else {
      scale.factor = 1
    }

    loadings.df$loading.x.scaled = loadings.df$loading.x * scale.factor
    loadings.df$loading.y.scaled = loadings.df$loading.y * scale.factor

    ## Select the top N
    n.loadings = min(n.loadings, nrow(loadings.df))
    loadings.top = loadings.df[order(loadings.df$distance, decreasing = TRUE),][seq_len(n.loadings),]

    return(list(loadings = loadings.df,
                top = loadings.top,
                scale.factor = scale.factor))
  }


# ----------------------------------------------------------------------------------------


#' @title .loading.layers
#'
#' @description Internal. Builds the ggplot layers (arrows and labels) used by the biplots.
#'
#' @param loadings.top data.frame of the top loadings, as returned by \link{.scale.loadings}.
#' @param loading.color String indicating the color of the arrows and labels.
#' @param loading.arrow.size Numeric value indicating the linewidth of the arrows.
#' @param loading.label.size Numeric value indicating the font size of the labels.
#' @param loading.alpha Numeric value (0-1) indicating the transparency of arrows and labels.
#'
#' @return A list of ggplot layers.
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @keywords internal

.loading.layers =
  function(loadings.top,
           loading.color = "turquoise4",
           loading.arrow.size = 0.6,
           loading.label.size = 3,
           loading.alpha = 0.7) {

    list(geom_segment(data = loadings.top,
                      aes(x = 0, y = 0,
                          xend = loading.x.scaled,
                          yend = loading.y.scaled),
                      inherit.aes = FALSE,
                      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                      color = loading.color,
                      linewidth = loading.arrow.size,
                      alpha = loading.alpha),

         ggrepel::geom_text_repel(data = loadings.top,
                                  aes(x = loading.x.scaled,
                                      y = loading.y.scaled,
                                      label = variable),
                                  inherit.aes = FALSE,
                                  color = loading.color,
                                  size = loading.label.size,
                                  alpha = loading.alpha,
                                  fontface = "italic",
                                  max.overlaps = 100,
                                  segment.color = loading.color,
                                  segment.alpha = loading.alpha * 0.5,
                                  segment.size = 0.3))
  }

# ----------------------------------------------------------------------------------------
