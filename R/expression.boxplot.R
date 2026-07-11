#' @title expression.boxplot
#'
#' @description Plots a boxplot of the expression of a specific protein. Samples can be groups depending on a metadata column. Optionally, the p-values of all the pairwise (2-by-2) comparisons between groups can be added on top of the plot.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param protein.id String indicating a protein for which plot the expression. The identifier must correspond to the full row.name of the counts table (equivalent to the \code{prot.id} column of the fold change table of \code{DEprot.analyses} object).
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'randomized', 'random', 'imp'. Default: \code{"imputed"}.
#' @param sample.subset Character vector indicating a subset of samples to display. The identifiers must correspond to a IDs in the \code{column.id} column of the object's metadata. Default: \code{NULL} (all samples are shown).
#' @param shape.column String indicating a column from the metadata table. This column will be used as factor for the shape of the points on the boxplot. Default: \code{NULL}: no different shapes.
#' @param group.by.metadata.column String indicating a column from the metadata table. This column will be used to define sample groups, and for each group it will be computed a mean of the counts. Default: \code{"column.id"} (no groups).
#' @param group.levels Ordered string vector indicating the order to use for the groups. Default: \code{NULL}, counts table order will be applied
#' @param scale.expression Logic value indicating whether Z-scores should be computed. Default: \code{FALSE} (no scaling).
#' @param x.label.angle Numeric value indicating the rotation angle to use for the x-axis labels. Default: \code{30}.
#' @param pairwise.comparisons Logical value indicating whether the p-values of all the pairwise (2-by-2) comparisons between the groups should be added on top of the boxplot using \code{ggpubr}. When \code{TRUE}, a comparison is computed for each possible pair of groups defined by \code{group.by.metadata.column}. Default: \code{FALSE}.
#' @param pairwise.test.type String indicating the statistical test to use for the pairwise comparisons. Any of the \code{ggpubr}-supported tests is accepted and the value is case/format-insensitive: capitalization, dots, spaces and hyphens are ignored (e.g. \code{"wilcox.test"}, \code{"Wilcoxon"}, \code{"WILCOX"}, \code{"mann-whitney"} are all equivalent). Supported families: \code{"t.test"} (Student's/Welch t-test), \code{"wilcox.test"} (Wilcoxon/Mann-Whitney), \code{"anova"} and \code{"kruskal.test"}. Since the comparisons are performed 2-by-2, the two latter are applied through their exact two-sample equivalents (\code{"anova"} -> pooled t-test, \code{"kruskal.test"} -> Wilcoxon rank-sum). Default: \code{"wilcox.test"}.
#' @param pairwise.include.ns Logical value indicating whether the non-significant comparisons (p > 0.05) should be displayed. If \code{FALSE}, only the significant comparisons are shown. Default: \code{TRUE}.
#' @param pairwise.p.label String indicating how the pairwise p-values should be displayed (case/format-insensitive). Use \code{"p.signif"} (aliases: \code{"stars"}, \code{"significance"}, \code{"symbol"}) to show the significance symbols (\code{ns}, \code{*}, \code{**}, \code{***}, \code{****}) drawn directly by \code{ggpubr::stat_compare_means}, or \code{"p.value"} (aliases: \code{"number"}, \code{"numeric"}, \code{"exact"}) to show the numeric p-value. Default: \code{"p.signif"}.
#' @param pairwise.p.decimals Numeric value indicating the number of decimals used to approximate the numeric p-values (used only when \code{pairwise.p.label} shows the numeric value). Values below 0.1 are rendered in scientific notation with a superscript exponent, e.g. 3.20\out{&times;10<sup>-2</sup>}. The actual (real) p-value is always displayed, so the uninformative \code{p < 2.2e-16} is never shown. Default: \code{2}.
#'
#' @return A boxplot of class ggplot2.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom ggpubr stat_compare_means stat_pvalue_manual
#' @import ggtext
#' @importFrom stats sd t.test wilcox.test
#' @importFrom utils combn
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' # Expression for all samples of protein 'protein.44'
#' expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                    protein.id = "protein.44",
#'                    shape.column = "replicate")
#'
#'
#' # Expression of protein 'protein.44' grouped by condition (combined.id)
#' expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                    protein.id = "protein.44",
#'                    group.by.metadata.column = "combined.id")
#'
#'
#' # Expression of protein 'protein.44' grouped by condition (combined.id) and Z-scored
#' expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                    protein.id = "protein.44",
#'                    group.by.metadata.column = "combined.id",
#'                    scale.expression = TRUE)
#'
#'
#' # Pairwise comparisons between conditions (significance symbols, Wilcoxon test)
#' expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                    protein.id = "protein.44",
#'                    group.by.metadata.column = "combined.id",
#'                    pairwise.comparisons = TRUE)
#'
#'
#' # Pairwise comparisons showing the exact numeric p-value
#' expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                    protein.id = "protein.44",
#'                    group.by.metadata.column = "combined.id",
#'                    pairwise.comparisons = TRUE,
#'                    pairwise.test.type = "t.test",
#'                    pairwise.p.label = "p.value",
#'                    pairwise.p.decimals = 3,
#'                    pairwise.include.ns = FALSE)
#'
#'
#' @export expression.boxplot




expression.boxplot =
  function(DEprot.object,
           protein.id,
           which.data = "imputed",
           sample.subset = NULL,
           shape.column = NULL,
           group.by.metadata.column = "column.id",
           group.levels = NULL,
           scale.expression = FALSE,
           x.label.angle = 30,
           pairwise.comparisons = FALSE,
           pairwise.test.type = "wilcox.test",
           pairwise.include.ns = TRUE,
           pairwise.p.label = "p.signif",
           pairwise.p.decimals = 2) {


    # ### Libraries
    # require(dplyr)
    # require(ggplot2)


    ### Internal functions
    check.matrix =
      function(m){
        warn = "Upon subsetting, no values to show are left."
        if (!is.logical(m)) {
          if (nrow(m) == 0 | ncol(m) == 0) {
            stop(warn)
            #return(return(invisible()))
          }
        } else {
          stop(warn)
          #return(return(invisible()))
        }
      }



    is.nan_df = function(data.frame) {do.call(cbind, lapply(data.frame, is.nan))}

    ######################################################################################

    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.object))) {
      if (!("DEprot" %in% class(DEprot.object))) {
        stop("The input must be an object of class 'DEprot' or 'DEprot.analyses'.")
        #return(invisible())
      }
    }


    ### check grouping column
    if (!is.null(group.by.metadata.column)) {
      if (!(group.by.metadata.column %in% colnames(DEprot.object@metadata))) {
        stop(paste0("The 'group.by.metadata.column' is not present in the metadata of the object provided.\n",
                    "       Available column IDs: ", paste0(colnames(DEprot.object@metadata), collapse = ", ")))
        #return(invisible())
      } else {
        meta = DEprot.object@metadata
      }
    }



    ### Check and extract table
    if (tolower(which.data) == "raw") {
      if (!is.null(DEprot.object@raw.counts)) {
        mat = DEprot.object@raw.counts
        data.used = "raw"
      } else {
        stop(paste0("Use of RAW counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
        #return(DEprot.object)
      }
    } else if (tolower(which.data) %in% c("norm", "normalized", "normal")) {
      if (!is.null(DEprot.object@norm.counts)) {
        mat = DEprot.object@norm.counts
        data.used = "normalized"
      } else {
        stop(paste0("Use of NORMALIZED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
        #return(DEprot.object)
      }
    } else if (tolower(which.data) %in% c("imputed", "imp", "impute")) {
      if (!is.null(DEprot.object@imputed.counts)) {
        mat = DEprot.object@imputed.counts
        data.used = "imputed"
      } else {
        stop(paste0("Use of IMPUTED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
        #return(DEprot.object)
      }
    } else if (tolower(which.data) %in% c("randomized", "random")) {
      if (!is.null(DEprot.object@random.counts)) {
        mat = DEprot.object@random.counts
        data.used = "randomized"
      } else {
        stop(paste0("Use of RANDOMIZED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
        #return(DEprot.object)
      }
    } else {
      stop(paste0("The 'which.data' value is not recognized.\n",
                  "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
      #return(DEprot.object)
    }



    ### Filter table of counts (samples and protein)
    if (!is.null(sample.subset)) {
      mat.filtered = mat[,which(colnames(mat) %in% sample.subset), drop=F]
    } else {
      mat.filtered = mat
    }
    check.matrix(mat.filtered)


    if (protein.id %in% rownames(mat.filtered)) {
      mat.filtered = mat.filtered[rownames(mat.filtered) == protein.id,,drop=F]
    } else {
      stop(paste0("The protein '", protein.id,"' is not present in the dataset."))
      #return(invisible())
    }



    ### reshape table
    exp.tb = as.data.frame(t(mat.filtered))
    colnames(exp.tb)[1] = "expression"
    exp.tb$column.id = rownames(exp.tb)

    # scale/center (z.score)
    if (scale.expression == TRUE) {
      exp.tb = dplyr::mutate(exp.tb, expression = (expression - mean(expression, na.rm = TRUE)) / sd(expression, na.rm = TRUE))
    }



    ### add group column
    if (is.null(group.by.metadata.column)) {group.by.metadata.column = "column.id"}

    if (group.by.metadata.column != "column.id") {
      exp.tb =
        dplyr::left_join(x = exp.tb,
                         y = meta[,c("column.id", group.by.metadata.column)],
                         by = "column.id")

      colnames(exp.tb)[ncol(exp.tb)] = "group"
    } else {
      exp.tb$group = exp.tb$column.id
    }



    ### add replicate column
    if (!is.null(shape.column)) {
      if (shape.column %in% colnames(meta)) {
        exp.tb =
          dplyr::left_join(x = exp.tb,
                           y = meta[,c("column.id", shape.column)],
                           by = "column.id")
        colnames(exp.tb)[ncol(exp.tb)] = "shape"
      } else {
        stop("The 'shape.column' provided is not present in the metadata table.")
        #return(invisible())
      }
    }



    ## add levels to group.column
    if (!is.null(group.levels)) {
      if (all(unique(exp.tb$group) %in% unique(group.levels))) {
        exp.tb = dplyr::mutate(.data = exp.tb, group = factor(group, levels = group.levels))
      } else {
        stop("The 'group.levels' do not include all the groups in the 'group.by.metadata.column'.")
        #return(invisible())
      }
    }




    ### Vertical offset for the global (Kruskal-Wallis/Wilcoxon) test label,
    ### so that it does not overlap with the pairwise p-value brackets
    kw.label.y = NULL
    if (isTRUE(pairwise.comparisons)) {
      if (is.factor(exp.tb$group)) {
        .pw.groups = levels(exp.tb$group)
      } else {
        .pw.groups = unique(as.character(exp.tb$group))
      }
      .pw.finite = exp.tb[is.finite(exp.tb$expression),,drop=F]
      .pw.sizes = table(as.character(.pw.finite$group))
      .pw.usable = .pw.groups[.pw.groups %in% names(.pw.sizes)[.pw.sizes >= 2]]

      if (length(.pw.usable) >= 2) {
        .pw.ncomparisons = choose(length(.pw.usable), 2)
        .pw.range = range(.pw.finite$expression, na.rm = TRUE)
        .pw.span = diff(.pw.range)
        if (!is.finite(.pw.span) || .pw.span == 0) {.pw.span = ifelse(.pw.range[2] == 0, 1, abs(.pw.range[2]))}
        # reserve room above the highest bracket (covers both stars and numeric layouts)
        kw.label.y = .pw.range[2] + (.pw.span * (0.10 + (0.13 * .pw.ncomparisons)))
      }
    }



    ### Generate boxplot
    boxplot =
      ggplot(data = exp.tb,
             aes(x = group,
                 y = expression,
                 fill = group,
                 color = group))

    if (scale.expression == TRUE) {
      boxplot = boxplot + geom_hline(yintercept = 0)
    }

    if (group.by.metadata.column != "column.id") {
      boxplot =
        boxplot +
        geom_boxplot(alpha = 0.25,
                     outliers = F,
                     show.legend = FALSE)
    }


    if (!is.null(shape.column)) {
      boxplot =
        boxplot +
        geom_point(aes(shape = factor(shape)),
                   #stroke = NA,
                   size = 3,
                   alpha = 0.5,
                   position = position_jitter(width = 0.15, height = 0),
                   show.legend = T) +
        guides(shape = guide_legend(title = shape.column))}
    else {
      boxplot =
        boxplot +
        geom_point(stroke = NA,
                   size = 3,
                   alpha = 0.5,
                   position = position_jitter(width = 0.15, height = 0),
                   show.legend = FALSE)
    }


    boxplot =
      boxplot +
      ggtitle(paste0("**",protein.id,"**")) +
      xlab(NULL) +
      ylab(ifelse(test = scale.expression == TRUE,
                  yes = paste0("centered log<sub>",DEprot.object@log.base,"</sub>(expression)"),
                  no = paste0("log<sub>",DEprot.object@log.base,"</sub>(expression)"))) +
      ggpubr::stat_compare_means(method = ifelse(test = length(unique(exp.tb$group)) == 2, yes = "wilcox", no = "kruskal"), label.y = kw.label.y, show.legend = FALSE) +
      guides(color = "none", fill = "none") +
      theme_classic() +
      theme(axis.title.x = ggtext::element_markdown(color = "black"),
            axis.title.y = ggtext::element_markdown(color = "black"),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            axis.text.x = element_text(color = "black", angle = x.label.angle, hjust = ifelse(x.label.angle %in% c(0), yes = 0.5, no = 1)),
            axis.text.y = element_text(color = "black"),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = "black"))



    ### Add pairwise (2-by-2) comparisons
    if (isTRUE(pairwise.comparisons)) {

      ## normalize the requested test type (case/format insensitive: dots, spaces, hyphens and apostrophes are ignored)
      .test.key = gsub("[[:space:]._'\u2019-]", "", tolower(trimws(pairwise.test.type)))

      requested.test =
        if (.test.key %in% c("ttest","t","student","students","studentt","studentst","welch","welcht","welchttest","unpairedttest","pairedttest","studentttest")) {
          "t.test"
        } else if (.test.key %in% c("wilcox","wilcoxon","wilcoxtest","wilcoxontest","mannwhitney","mannwhitneyu","mannwhitneyutest","mww","wmw","utest","u","ranksum","wilcoxonranksum")) {
          "wilcox.test"
        } else if (.test.key %in% c("anova","aov","onewayanova","oneway","ftest","f")) {
          "anova"
        } else if (.test.key %in% c("kruskal","kruskalwallis","kruskaltest","kruskalwallistest","kw")) {
          "kruskal.test"
        } else {
          stop(paste0("The 'pairwise.test.type' value ('", pairwise.test.type, "') is not recognized.\n",
                      "       Supported tests: 't.test', 'wilcox.test', 'anova', 'kruskal.test' (case/format insensitive)."))
        }

      ## two-sample equivalent applied to each 2-by-2 comparison (kept consistent between label styles).
      ## For exactly two groups: one-way ANOVA == pooled (var.equal) t-test, Kruskal-Wallis == Wilcoxon rank-sum.
      two.sample =
        switch(requested.test,
               "t.test"       = list(method = "t.test",      args = list()),
               "wilcox.test"  = list(method = "wilcox.test", args = list(exact = FALSE)),
               "anova"        = list(method = "t.test",      args = list(var.equal = TRUE)),
               "kruskal.test" = list(method = "wilcox.test", args = list(exact = FALSE, correct = FALSE)))

      ## normalize the label style (significance symbols vs numeric p-value)
      .label.key = gsub("[[:space:]._'\u2019-]", "", tolower(trimws(pairwise.p.label)))

      label.style =
        if (.label.key %in% c("psignif","signif","stars","star","asterisk","asterisks","significance","symbol","symbols","sign","star(s)")) {
          "stars"
        } else if (.label.key %in% c("pvalue","p","pval","value","number","numeric","num","exact","pformat")) {
          "pvalue"
        } else {
          stop(paste0("The 'pairwise.p.label' value ('", pairwise.p.label, "') is not recognized.\n",
                      "       Use 'p.signif' (significance symbols) or 'p.value' (numeric p-value)."))
        }

      ## number of decimals to approximate the numeric p-value
      p.dec = max(0L, as.integer(round(pairwise.p.decimals)))

      ## groups usable for a pairwise test = groups with at least two finite values (ordered as displayed)
      if (is.factor(exp.tb$group)) {
        ordered.groups = levels(exp.tb$group)
      } else {
        ordered.groups = unique(as.character(exp.tb$group))
      }

      finite.tb = exp.tb[is.finite(exp.tb$expression),,drop=F]
      group.sizes = table(as.character(finite.tb$group))
      usable.groups = ordered.groups[ordered.groups %in% names(group.sizes)[group.sizes >= 2]]


      if (length(usable.groups) < 2) {
        warning("Pairwise comparisons were requested but less than two groups with at least two values are available: no comparison is shown.")
      } else {

        ## all the possible 2-by-2 comparisons
        comparisons.list = utils::combn(usable.groups, 2, simplify = FALSE)


        if (label.style == "stars") {
          ##### Significance symbols: drawn directly through ggpubr::stat_compare_means (unadjusted pairwise p-values)
          boxplot =
            boxplot +
            ggpubr::stat_compare_means(data = finite.tb,
                                       comparisons = comparisons.list,
                                       method = two.sample$method,
                                       method.args = two.sample$args,
                                       label = "p.signif",
                                       hide.ns = !isTRUE(pairwise.include.ns),
                                       tip.length = 0.01,
                                       size = 3.3)

        } else {
          ##### Numeric p-value: computed exactly (the real value is always shown, never 'p < 2.2e-16'),
          #####                  formatted manually (custom decimals + scientific superscript when < 0.1),
          #####                  and drawn with ggpubr::stat_pvalue_manual.

          ## p-value formatter: returns a plain string using Unicode superscripts, e.g. "3.20\u00d710\u207b\u00b2"
          format.pairwise.p =
            function(p, decimals) {
              if (is.na(p)) {return("NA")}

              superscript = c("0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2", "3" = "\u00b3", "4" = "\u2074",
                              "5" = "\u2075", "6" = "\u2076", "7" = "\u2077", "8" = "\u2078", "9" = "\u2079",
                              "-" = "\u207b")
              to.superscript = function(n) {paste0(superscript[strsplit(as.character(n), "")[[1]]], collapse = "")}

              prefix = ""
              if (p <= 0) {p = .Machine$double.xmin; prefix = "< "} # numeric underflow safeguard

              exponent = floor(log10(p))

              if (exponent <= -2) {
                # scientific notation with superscript exponent (e.g. 3.20 x 10^-2)
                mantissa = round(p / (10^exponent), decimals)
                if (mantissa >= 10) {mantissa = mantissa / 10; exponent = exponent + 1}
                lab = paste0(formatC(mantissa, format = "f", digits = decimals), "\u00d7", "10", to.superscript(exponent))
              } else {
                # plain decimal notation (0.1 <= p <= 1)
                lab = formatC(round(p, decimals), format = "f", digits = decimals)
              }

              return(paste0(prefix, lab))
            }


          ## exact per-pair p-values
          pairwise.tb =
            do.call(rbind,
                    lapply(comparisons.list,
                           function(pair) {
                             d1 = finite.tb$expression[as.character(finite.tb$group) == pair[1]]
                             d2 = finite.tb$expression[as.character(finite.tb$group) == pair[2]]
                             pval = tryCatch(expr = suppressWarnings(do.call(two.sample$method, c(list(x = d1, y = d2), two.sample$args))$p.value),
                                             error = function(e){return(NA_real_)})
                             data.frame(group1 = pair[1], group2 = pair[2], p.value = pval, stringsAsFactors = FALSE)
                           }))

          ## keep only computable comparisons (and, if required, only the significant ones)
          pairwise.tb = pairwise.tb[!is.na(pairwise.tb$p.value),,drop=F]
          if (!isTRUE(pairwise.include.ns)) {
            pairwise.tb = pairwise.tb[pairwise.tb$p.value <= 0.05,,drop=F]
          }


          if (nrow(pairwise.tb) > 0) {
            ## y positions of the brackets (stacked above the data)
            y.range = range(finite.tb$expression, na.rm = TRUE)
            y.span = diff(y.range)
            if (!is.finite(y.span) || y.span == 0) {y.span = ifelse(y.range[2] == 0, 1, abs(y.range[2]))}
            pairwise.tb$y.position = y.range[2] + (0.08 * y.span) + ((seq_len(nrow(pairwise.tb)) - 1) * (0.09 * y.span))

            ## formatted labels
            pairwise.tb$p.label = vapply(X = pairwise.tb$p.value,
                                         FUN = function(x){format.pairwise.p(p = x, decimals = p.dec)},
                                         FUN.VALUE = character(1))

            boxplot =
              boxplot +
              ggpubr::stat_pvalue_manual(data = pairwise.tb[,c("group1", "group2", "y.position", "p.label")],
                                         label = "p.label",
                                         xmin = "group1",
                                         xmax = "group2",
                                         y.position = "y.position",
                                         tip.length = 0.01,
                                         size = 3.3,
                                         bracket.size = 0.3,
                                         inherit.aes = FALSE)
          }
        }
      }
    }



    ### return plot
    return(boxplot)
  } # END function
