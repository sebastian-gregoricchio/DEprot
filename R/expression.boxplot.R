#' @title expression.boxplot
#'
#' @description Plots a boxplot of the expression of one or more proteins. Samples can be grouped depending on a metadata column. Optionally, the p-values of all the pairwise (2-by-2) comparisons between groups can be added on top of the plot.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param protein.id String (or vector of strings) indicating the protein(s) for which plot the expression. The identifiers must correspond to the full row.names of the counts table (equivalent to the \code{prot.id} column of the fold change table of \code{DEprot.analyses} object). When several proteins are provided the plot is faceted, one panel per protein.
#' @param which.data String indicating which type of counts should be used. One among: 'raw', 'normalized', 'norm', 'imputed', 'randomized', 'random', 'imp'. Default: \code{"imputed"}.
#' @param sample.subset Character vector indicating a subset of samples to display. The identifiers must correspond to a IDs in the \code{column.id} column of the object's metadata. Default: \code{NULL} (all samples are shown).
#' @param shape.column String indicating a column from the metadata table. This column will be used as factor for the shape of the points on the boxplot. Default: \code{NULL}: no different shapes.
#' @param group.by.metadata.column String indicating a column from the metadata table. This column will be used to define sample groups, and for each group it will be computed a mean of the counts. Default: \code{"column.id"} (no groups).
#' @param group.levels Ordered string vector indicating the order to use for the groups. Default: \code{NULL}, counts table order will be applied
#' @param scale.expression Logic value indicating whether Z-scores should be computed. With several proteins the Z-score is computed protein by protein, so that the panels stay comparable. Default: \code{FALSE} (no scaling).
#' @param x.label.angle Numeric value indicating the rotation angle to use for the x-axis labels. Default: \code{30}.
#' @param ncol Numeric value indicating the number of columns of the facet grid, used only when several proteins are plotted. Default: \code{NULL} (automatic).
#' @param free.y Logic value indicating whether each panel should have its own y-axis, used only when several proteins are plotted. Proteins of different abundance are otherwise flattened onto a common scale. Default: \code{TRUE}.
#' @param pairwise.comparisons Logical value indicating whether the p-values of all the pairwise (2-by-2) comparisons between the groups should be added on top of the boxplot using \code{ggpubr}. When \code{TRUE}, a comparison is computed for each possible pair of groups defined by \code{group.by.metadata.column}, independently within each protein. Default: \code{FALSE}.
#' @param pairwise.test.type String indicating the statistical test to use for the pairwise comparisons. The value is case/format-insensitive: capitalization, dots, spaces and hyphens are ignored (e.g. \code{"wilcox.test"}, \code{"Wilcoxon"}, \code{"WILCOX"}, \code{"mann-whitney"} are all equivalent). Supported families: \code{"t.test"} (Student's/Welch t-test), \code{"wilcox.test"} (Wilcoxon/Mann-Whitney), \code{"anova"} and \code{"kruskal.test"}. Since the comparisons are performed 2-by-2, the two latter are applied through their two-sample equivalents (\code{"anova"} -> pooled t-test, \code{"kruskal.test"} -> Wilcoxon rank-sum). When the pairwise brackets are displayed, the global p-value shown at the top of the plot is computed with the same test and the same arguments: with two groups the two labels report therefore the same value, while with more groups the multi-sample version of the family is used (\code{"anova"} or \code{"kruskal.test"}). The Wilcoxon p-values follow the defaults of \code{wilcox.test}: exact when both groups have less than 50 values (older R versions require also the absence of ties), normal approximation with continuity correction otherwise. With few replicates the exact p-value has a floor: two groups of 3 samples cannot go below 0.1, even when they do not overlap at all. Names of paired tests are accepted too (e.g. \code{"paired t-test"}, \code{"signed-rank"}, \code{"Friedman"}) and switch on \code{paired.test}. Default: \code{"wilcox.test"}.
#' @param pairwise.include.ns Logical value indicating whether the non-significant comparisons (p > 0.05) should be displayed. If \code{FALSE}, only the significant comparisons are shown. Default: \code{TRUE}.
#' @param pairwise.p.label String indicating how the pairwise p-values should be displayed (case/format-insensitive). Use \code{"p.signif"} (aliases: \code{"stars"}, \code{"significance"}, \code{"symbol"}) to show the significance symbols (\code{ns}: p > 0.05, \code{*}: p <= 0.05, \code{**}: p <= 0.01, \code{***}: p <= 0.001, \code{****}: p <= 0.0001), or \code{"p.value"} (aliases: \code{"number"}, \code{"numeric"}, \code{"exact"}) to show the numeric p-value. The symbols are derived from the same p-values of the numeric labels. Default: \code{"p.signif"}.
#' @param pairwise.p.decimals Numeric value indicating the number of decimals used to approximate the numeric p-values (used only when \code{pairwise.p.label} shows the numeric value). Values below 0.1 are rendered in scientific notation with a superscript exponent, e.g. 3.20\out{&times;10<sup>-2</sup>}. The actual (real) p-value is always displayed, so the uninformative \code{p < 2.2e-16} is never shown. Default: \code{2}.
#' @param paired.test Logical value indicating whether paired statistical tests should be performed. The samples are matched between groups through the IDs stored in \code{replicate.column}, as in \code{diff.analyses}, and each comparison uses only the replicates measured in both groups. The t-tests become paired t-tests and the Wilcoxon rank-sum test a Wilcoxon signed-rank test. With more than two groups the global p-value comes from a repeated-measures ANOVA (\code{"t.test"}, \code{"anova"}) or from a Friedman test (\code{"wilcox.test"}, \code{"kruskal.test"}), computed on the replicates measured in all the groups. The pairing applies to the global p-value also when \code{pairwise.comparisons = FALSE}. Default: \code{FALSE}.
#' @param replicate.column String indicating the name of a column from the metadata table in which are stored the replicate IDs used to pair the samples. It is required when \code{paired.test = TRUE}, or when \code{pairwise.test.type} indicates a paired test. A replicate ID cannot be repeated within a group. Default: \code{NULL}.
#'
#' @return A boxplot of class ggplot2, faceted by protein when several proteins are provided.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom ggpubr stat_pvalue_manual
#' @import ggtext
#' @importFrom stats anova friedman.test kruskal.test lm oneway.test sd t.test wilcox.test
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
#' # Expression of several proteins: one panel per protein
#' expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                    protein.id = c("protein.44", "protein.45", "protein.46"),
#'                    group.by.metadata.column = "combined.id",
#'                    ncol = 3)
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
#'                    pairwise.p.label = "stars",
#'                    pairwise.p.decimals = 3,
#'                    pairwise.include.ns = FALSE)
#'
#'
#' # Paired t-test: the samples of the different conditions are matched through their replicate ID
#' expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                    protein.id = "protein.44",
#'                    group.by.metadata.column = "combined.id",
#'                    shape.column = "replicate",
#'                    pairwise.comparisons = TRUE,
#'                    pairwise.test.type = "t.test",
#'                    paired.test = TRUE,
#'                    replicate.column = "replicate")
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
           ncol = NULL,
           free.y = TRUE,
           pairwise.comparisons = FALSE,
           pairwise.test.type = "wilcox.test",
           pairwise.include.ns = TRUE,
           pairwise.p.label = "p.signif",
           pairwise.p.decimals = 2,
           paired.test = FALSE,
           replicate.column = NULL) {


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



    ## two-sample equivalent of each supported test family, with the arguments used to compute it.
    ## Both the global p-value and the 2-by-2 brackets go through this function, hence with two groups
    ## the two labels cannot disagree.
    ## For exactly two groups: one-way ANOVA == pooled (var.equal) t-test, Kruskal-Wallis == Wilcoxon rank-sum.
    ## The Wilcoxon test keeps the defaults of wilcox.test(), i.e. the exact p-value whenever it can be computed.
    ## With few replicates the normal approximation is anti-conservative: two groups of 3 samples that do not
    ## overlap give 0.081 (0.0495 without continuity correction), while the exact p-value cannot go below 0.1.
    ## In a paired design the t-tests become paired t-tests and the rank-sum test a signed-rank test.
    two.sample.equivalent =
      function(test, paired = FALSE) {
        switch(test,
               "t.test"       = list(method = "t.test",      args = list(paired = paired)),
               "wilcox.test"  = list(method = "wilcox.test", args = list(paired = paired)),
               "anova"        = list(method = "t.test",      args = list(paired = paired, var.equal = TRUE)),
               "kruskal.test" = list(method = "wilcox.test", args = list(paired = paired)))
      }



    ## values of two groups entering a two-sample test. In a paired design the two vectors are aligned on
    ## the replicate IDs, and a replicate missing in one of the groups is dropped from both.
    extract.pair =
      function(tb, group.1, group.2, paired = FALSE) {
        d1 = tb[as.character(tb$group) == group.1,,drop=F]
        d2 = tb[as.character(tb$group) == group.2,,drop=F]

        if (isTRUE(paired)) {
          shared.reps = intersect(d1$pair.id[!is.na(d1$pair.id)], d2$pair.id[!is.na(d2$pair.id)])
          return(list(x = d1$expression[match(shared.reps, d1$pair.id)],
                      y = d2$expression[match(shared.reps, d2$pair.id)]))
        } else {
          return(list(x = d1$expression, y = d2$expression))
        }
      }



    ## p-value of a two-sample test, NA when it cannot be computed (e.g. less than two values per group)
    two.sample.p =
      function(values, test) {
        if (length(values$x) < 2 | length(values$y) < 2) {return(NA_real_)}

        pval = tryCatch(expr = suppressWarnings(do.call(test$method, c(list(x = values$x, y = values$y), test$args))$p.value),
                        error = function(e){return(NA_real_)})
        return(pval)
      }



    ## p-value of the multi-sample version of the family, used with more than two groups.
    ## In a paired design only the replicates measured in all the groups are kept (complete blocks), and the
    ## repeated-measures ANOVA is the F-test of the group term once the replicate effect is removed.
    multi.sample.p =
      function(tb, groups, test, paired = FALSE) {
        tb = tb[as.character(tb$group) %in% groups,,drop=F]
        tb$group = factor(as.character(tb$group), levels = groups)

        if (isTRUE(paired)) {
          rep.counts = table(tb$pair.id)
          tb = tb[tb$pair.id %in% names(rep.counts)[rep.counts == length(groups)],,drop=F]
          if (length(unique(tb$pair.id)) < 2) {return(NA_real_)}
          tb$pair.id = factor(tb$pair.id)
        }

        pval =
          tryCatch(expr =
                     suppressWarnings(
                       switch(paste0(test, ifelse(test = isTRUE(paired), yes = ".paired", no = "")),
                              "anova"               = stats::oneway.test(expression ~ group, data = tb, var.equal = TRUE)$p.value,
                              "kruskal.test"        = stats::kruskal.test(expression ~ group, data = tb)$p.value,
                              "anova.paired"        = stats::anova(stats::lm(expression ~ pair.id + group, data = tb))["group", "Pr(>F)"],
                              "kruskal.test.paired" = stats::friedman.test(y = tb$expression, groups = tb$group, blocks = tb$pair.id)$p.value)),
                   error = function(e){return(NA_real_)})
        return(pval)
      }



    ## name of each test in the global label (for the unpaired tests the same names used by ggpubr)
    test.names = c("t.test"              = "T-test",
                   "wilcox.test"         = "Wilcoxon",
                   "anova"               = "Anova",
                   "kruskal.test"        = "Kruskal-Wallis",
                   "t.test.paired"       = "Paired t-test",
                   "wilcox.test.paired"  = "Wilcoxon signed-rank",
                   "anova.paired"        = "Repeated-measures Anova",
                   "kruskal.test.paired" = "Friedman")

    ######################################################################################

    ### check object
    if (!(methods::is(DEprot.object, "DEprot.analyses"))) {
      if (!(methods::is(DEprot.object, "DEprot"))) {
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


    ### check the column used to pair the samples
    if (!is.null(replicate.column)) {
      if (!(replicate.column %in% colnames(DEprot.object@metadata))) {
        stop(paste0("The 'replicate.column' is not present in the metadata of the object provided.\n",
                    "       Available column IDs: ", paste0(colnames(DEprot.object@metadata), collapse = ", ")))
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
      }
    } else if (tolower(which.data) %in% c("norm", "normalized", "normal")) {
      if (!is.null(DEprot.object@norm.counts)) {
        mat = DEprot.object@norm.counts
        data.used = "normalized"
      } else {
        stop(paste0("Use of NORMALIZED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
      }
    } else if (tolower(which.data) %in% c("imputed", "imp", "impute")) {
      if (!is.null(DEprot.object@imputed.counts)) {
        mat = DEprot.object@imputed.counts
        data.used = "imputed"
      } else {
        stop(paste0("Use of IMPUTED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
      }
    } else if (tolower(which.data) %in% c("randomized", "random")) {
      if (!is.null(DEprot.object@random.counts)) {
        mat = DEprot.object@random.counts
        data.used = "randomized"
      } else {
        stop(paste0("Use of RANDOMIZED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
      }
    } else {
      stop(paste0("The 'which.data' value is not recognized.\n",
                  "       Please indicated a count type among 'raw', 'normalized', 'randomized', 'imputed', using the option 'which.data'."))
    }



    ### Filter table of counts (samples and proteins)
    if (!is.null(sample.subset)) {
      mat.filtered = mat[,which(colnames(mat) %in% sample.subset), drop=FALSE]
    } else {
      mat.filtered = mat
    }
    check.matrix(mat.filtered)


    ## the proteins are validated all at once: a single missing ID must not hide the others
    protein.id = as.character(protein.id)

    if (length(protein.id) == 0) {
      stop("Provide at least one 'protein.id'.")
      #return(invisible())
    }

    missing.proteins = setdiff(protein.id, rownames(mat.filtered))

    if (length(missing.proteins) > 0) {
      stop(paste0("The following protein(s) are not present in the dataset: ",
                  paste0(missing.proteins, collapse = ", "), "."))
      #return(invisible())
    }

    mat.filtered = mat.filtered[protein.id, , drop = FALSE]

    ## a single protein keeps the historical layout (title, no facet), several proteins are
    ## displayed as a panel per protein
    multiple.proteins = (length(protein.id) > 1)



    ### reshape table
    ## the Z-score is computed protein by protein: scaling all the proteins together would
    ## simply re-center the panels on the average abundance of the set
    exp.tb =
      do.call(rbind,
              lapply(protein.id,
                     function(prot) {
                       values = as.numeric(mat.filtered[prot,])

                       if (scale.expression == TRUE) {
                         value.sd = sd(values, na.rm = TRUE)
                         values =
                           if (is.na(value.sd) | value.sd == 0) {
                             values - mean(values, na.rm = TRUE)
                           } else {
                             (values - mean(values, na.rm = TRUE)) / value.sd
                           }
                       }

                       data.frame(prot.id = prot,
                                  column.id = colnames(mat.filtered),
                                  expression = values,
                                  stringsAsFactors = FALSE)
                     }))

    ## the panels follow the order in which the proteins were requested
    exp.tb$prot.id = factor(exp.tb$prot.id, levels = protein.id)



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




    ### Normalize the requested test type (case/format insensitive: dots, spaces, hyphens and apostrophes are ignored).
    ### This happens before the plot is built because the global p-value shown at the top must be computed
    ### with the very same test, and the very same arguments, used for the 2-by-2 comparisons.
    .test.key = gsub("[[:space:]._'\u2019-]", "", tolower(trimws(pairwise.test.type)))

    ## the names of explicitly unpaired or paired tests are listed apart, since they define the design as well
    unpaired.keys = c("welch","welcht","welchttest","unpairedt","unpairedttest",
                      "mannwhitney","mannwhitneyu","mannwhitneyutest","mww","wmw","utest","u","ranksum","ranksumtest","wilcoxonranksum","wilcoxonranksumtest")

    paired.keys = c("pairedt","pairedttest","pairedstudentttest","pairedstudentsttest",
                    "pairedwilcox","pairedwilcoxon","pairedwilcoxtest","pairedwilcoxontest","signedrank","signedranktest","wilcoxonsignedrank","wilcoxonsignedranktest",
                    "rmanova","repeatedmeasuresanova","pairedanova",
                    "friedman","friedmantest")

    requested.test =
      if (.test.key %in% c("ttest","t","student","students","studentt","studentst","studentttest","studentsttest",
                           "welch","welcht","welchttest","unpairedt","unpairedttest",
                           "pairedt","pairedttest","pairedstudentttest","pairedstudentsttest")) {
        "t.test"
      } else if (.test.key %in% c("wilcox","wilcoxon","wilcoxtest","wilcoxontest",
                                  "mannwhitney","mannwhitneyu","mannwhitneyutest","mww","wmw","utest","u","ranksum","ranksumtest","wilcoxonranksum","wilcoxonranksumtest",
                                  "pairedwilcox","pairedwilcoxon","pairedwilcoxtest","pairedwilcoxontest","signedrank","signedranktest","wilcoxonsignedrank","wilcoxonsignedranktest")) {
        "wilcox.test"
      } else if (.test.key %in% c("anova","aov","onewayanova","oneway","ftest","f",
                                  "rmanova","repeatedmeasuresanova","pairedanova")) {
        "anova"
      } else if (.test.key %in% c("kruskal","kruskalwallis","kruskaltest","kruskalwallistest","kw",
                                  "friedman","friedmantest")) {
        "kruskal.test"
      } else {
        stop(paste0("The 'pairwise.test.type' value ('", pairwise.test.type, "') is not recognized.\n",
                    "       Supported tests: 't.test', 'wilcox.test', 'anova', 'kruskal.test' (case/format insensitive)."))
      }



    ### Paired design: the samples are matched between groups through their replicate ID, as in diff.analyses().
    ### A test named as paired (e.g. "paired t-test") switches it on, while a test named as unpaired (e.g. "Welch",
    ### "Mann-Whitney") cannot be combined with it. As for the family, the name given in 'pairwise.test.type'
    ### counts only when the brackets are shown.
    paired = isTRUE(paired.test)

    if (isTRUE(pairwise.comparisons)) {
      if (.test.key %in% paired.keys) {
        paired = TRUE
      } else if (.test.key %in% unpaired.keys & paired) {
        stop(paste0("The 'pairwise.test.type' value ('", pairwise.test.type, "') indicates an unpaired test, while 'paired.test = TRUE'."))
      }
    }

    if (paired) {
      if (is.null(replicate.column)) {
        stop("A paired test was required, but no 'replicate.column' was provided: the samples cannot be matched between groups.")
      }

      exp.tb$pair.id = as.character(DEprot.object@metadata[[replicate.column]][match(exp.tb$column.id, DEprot.object@metadata$column.id)])

      ## a replicate ID repeated within a group would make the pairing ambiguous
      sample.reps = unique(exp.tb[!is.na(exp.tb$pair.id), c("column.id", "group", "pair.id"), drop=F])

      if (any(duplicated(sample.reps[, c("group", "pair.id"), drop=F]))) {
        stop("At least one replicate ID in the 'replicate.column' is duplicated within a group: the samples cannot be paired.")
      }
    }



    ### Values actually usable by a test: non-finite values (NA, NaN, and the -Inf coming from the
    ### log of a zero count) are dropped once and the same table is used for the global p-value,
    ### for the pairwise ones and for the positioning of the brackets.
    finite.tb = exp.tb[is.finite(exp.tb$expression),,drop=F]

    ## groups usable for a test = groups with at least two finite values (ordered as displayed)
    if (is.factor(exp.tb$group)) {
      ordered.groups = levels(exp.tb$group)
    } else {
      ordered.groups = unique(as.character(exp.tb$group))
    }

    group.sizes = table(as.character(finite.tb$group))
    usable.groups = ordered.groups[ordered.groups %in% names(group.sizes)[group.sizes >= 2]]

    ## in a paired design, two groups sharing less than two replicates cannot be compared
    if (paired & length(usable.groups) >= 2) {
      unmatched.groups =
        unlist(lapply(utils::combn(usable.groups, 2, simplify = FALSE),
                      function(pair) {
                        shared.reps = intersect(exp.tb$pair.id[as.character(exp.tb$group) == pair[1]],
                                                exp.tb$pair.id[as.character(exp.tb$group) == pair[2]])
                        if (sum(!is.na(shared.reps)) < 2) {return(paste0(pair[1], " vs ", pair[2]))}
                        return(NULL)
                      }))

      if (length(unmatched.groups) > 0) {
        warning(paste0("The following groups share less than two replicate IDs, hence they cannot be compared by a paired test: ",
                       paste0(unmatched.groups, collapse = ", "), "."))
      }
    }



    ## two-sample equivalent applied to each 2-by-2 comparison (kept consistent between label styles)
    two.sample = two.sample.equivalent(requested.test, paired = paired)



    ### Global p-value displayed at the top of each panel. It is computed here, and not by ggpubr, so that
    ### it follows the same rules of the brackets (test, arguments, pairing and values used).
    ### With two groups it is exactly the pairwise test, arguments included, so that the two labels
    ### cannot disagree; with more groups its multi-sample version is used (ANOVA or Kruskal-Wallis, and
    ### in a paired design repeated-measures ANOVA or Friedman test).
    ### The family follows 'pairwise.test.type' only when the brackets are shown, otherwise the
    ### historical Wilcoxon/Kruskal-Wallis behaviour is kept.
    global.family = ifelse(test = isTRUE(pairwise.comparisons), yes = requested.test, no = "wilcox.test")
    global.two.sample = two.sample.equivalent(global.family, paired = paired)
    global.multi.sample = ifelse(test = global.family %in% c("t.test", "anova"), yes = "anova", no = "kruskal.test")

    global.tb =
      do.call(rbind,
              lapply(protein.id,
                     function(prot) {
                       prot.tb = finite.tb[as.character(finite.tb$prot.id) == prot,,drop=F]

                       ## groups with at least two values for this protein (ordered as displayed)
                       prot.sizes = table(as.character(prot.tb$group))
                       prot.groups = ordered.groups[ordered.groups %in% names(prot.sizes)[prot.sizes >= 2]]

                       if (length(prot.groups) < 2) {return(NULL)}

                       if (length(prot.groups) == 2) {
                         pval = two.sample.p(values = extract.pair(tb = prot.tb, group.1 = prot.groups[1], group.2 = prot.groups[2], paired = paired),
                                             test = global.two.sample)
                         test.id = global.two.sample$method
                       } else {
                         pval = multi.sample.p(tb = prot.tb, groups = prot.groups, test = global.multi.sample, paired = paired)
                         test.id = global.multi.sample
                       }

                       if (is.na(pval)) {return(NULL)}

                       ## vertical position: when the brackets are shown the label goes above them, reserving
                       ## room for all the possible comparisons (covers both stars and numeric layouts);
                       ## otherwise it sits at the top of the data of the panel, as done by ggpubr
                       if (isTRUE(pairwise.comparisons) & length(usable.groups) >= 2) {
                         y.range = range(prot.tb$expression)
                         y.span = diff(y.range)
                         if (!is.finite(y.span) || y.span == 0) {y.span = ifelse(y.range[2] == 0, 1, abs(y.range[2]))}
                         label.y = y.range[2] + (y.span * (0.10 + (0.13 * choose(length(usable.groups), 2))))
                       } else {
                         panel.values = if (multiple.proteins & isTRUE(free.y)) {prot.tb$expression} else {finite.tb$expression}
                         label.y = max(c(panel.values, ifelse(test = scale.expression == TRUE, yes = 0, no = -Inf)))
                       }

                       p.text = ifelse(test = pval < 2.2e-16, yes = "p < 2.2e-16", no = paste("p =", signif(pval, 2)))

                       data.frame(prot.id = prot,
                                  x = 1,
                                  y = label.y,
                                  label = paste0(test.names[[paste0(test.id, ifelse(test = paired, yes = ".paired", no = ""))]], ", ", p.text),
                                  stringsAsFactors = FALSE)
                     }))



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
      ggtitle(switch(multiple.proteins + 1, paste0("**",protein.id,"**"), NULL)) +
      xlab(NULL) +
      ylab(ifelse(test = scale.expression == TRUE,
                  yes = paste0("centered log<sub>",DEprot.object@log.base,"</sub>(expression)"),
                  no = paste0("log<sub>",DEprot.object@log.base,"</sub>(expression)"))) +
      guides(color = "none", fill = "none") +
      theme_classic() +
      theme(axis.title.x = ggtext::element_markdown(color = "black"),
            axis.title.y = ggtext::element_markdown(color = "black"),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            axis.text.x = element_text(color = "black", angle = x.label.angle, hjust = ifelse(x.label.angle %in% c(0), yes = 0.5, no = 1)),
            axis.text.y = element_text(color = "black"),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = "black"),
            ## naked facet labels, the protein name in bold
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"))


    ## one panel per protein: the proteins are rarely expressed in the same range, hence
    ## each panel gets its own y-axis unless required otherwise
    if (multiple.proteins) {
      boxplot =
        boxplot +
        facet_wrap(~ prot.id,
                   ncol = ncol,
                   scales = ifelse(isTRUE(free.y), yes = "free_y", no = "fixed"))
    }



    ### Add the global p-value: one label per panel, placed where ggpubr::stat_compare_means() put it
    ### (above the first group, left-aligned)
    if (!is.null(global.tb)) {
      global.tb$prot.id = factor(global.tb$prot.id, levels = protein.id)

      boxplot =
        boxplot +
        geom_text(data = global.tb,
                  mapping = aes(x = x, y = y, label = label),
                  hjust = 0.2,
                  vjust = 0,
                  inherit.aes = FALSE,
                  show.legend = FALSE)
    }



    ### Add pairwise (2-by-2) comparisons
    if (isTRUE(pairwise.comparisons)) {

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

      ## 'finite.tb' and 'usable.groups' come from the section above: the pairwise brackets
      ## and the global label must rely on the same values and on the same groups
      if (length(usable.groups) < 2) {
        warning("Pairwise comparisons were requested but less than two groups with at least two values are available: no comparison is shown.")
      } else {

        ## all the possible 2-by-2 comparisons
        comparisons.list = utils::combn(usable.groups, 2, simplify = FALSE)


        ##### The p-values are computed here for both label styles, protein by protein, through the same
        ##### functions used for the global p-value: the significance symbols are derived from the numbers
        ##### shown by the numeric labels. The real value is always available (never 'p < 2.2e-16'), the
        ##### numeric labels are formatted manually (custom decimals + scientific superscript when < 0.1),
        ##### and the brackets are drawn with ggpubr::stat_pvalue_manual.

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


        ## per-pair p-values, computed within each protein
        pairwise.tb =
          do.call(rbind,
                  lapply(protein.id,
                         function(prot) {
                           prot.tb = finite.tb[as.character(finite.tb$prot.id) == prot,,drop=F]

                           if (nrow(prot.tb) == 0) {return(NULL)}

                           prot.pairs =
                             do.call(rbind,
                                     lapply(comparisons.list,
                                            function(pair) {
                                              pval = two.sample.p(values = extract.pair(tb = prot.tb, group.1 = pair[1], group.2 = pair[2], paired = paired),
                                                                  test = two.sample)
                                              data.frame(prot.id = prot, group1 = pair[1], group2 = pair[2], p.value = pval, stringsAsFactors = FALSE)
                                            }))

                           ## keep only computable comparisons (and, if required, only the significant ones)
                           prot.pairs = prot.pairs[!is.na(prot.pairs$p.value),,drop=F]
                           if (!isTRUE(pairwise.include.ns)) {
                             prot.pairs = prot.pairs[prot.pairs$p.value <= 0.05,,drop=F]
                           }

                           if (nrow(prot.pairs) == 0) {return(NULL)}

                           ## y positions of the brackets (stacked above the data of THIS protein)
                           y.range = range(prot.tb$expression, na.rm = TRUE)
                           y.span = diff(y.range)
                           if (!is.finite(y.span) || y.span == 0) {y.span = ifelse(y.range[2] == 0, 1, abs(y.range[2]))}
                           prot.pairs$y.position = y.range[2] + (0.08 * y.span) + ((seq_len(nrow(prot.pairs)) - 1) * (0.09 * y.span))

                           return(prot.pairs)
                         }))


        if (!is.null(pairwise.tb)) {
          if (nrow(pairwise.tb) > 0) {
            ## labels: significance symbols (same thresholds of the 'pairwise.include.ns' filter) or formatted p-values
            if (label.style == "stars") {
              pairwise.tb$p.label = as.character(cut(x = pairwise.tb$p.value,
                                                     breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                     labels = c("****", "***", "**", "*", "ns")))
            } else {
              pairwise.tb$p.label = vapply(X = pairwise.tb$p.value,
                                           FUN = function(x){format.pairwise.p(p = x, decimals = p.dec)},
                                           FUN.VALUE = character(1))
            }

            ## the facetting variable must be carried over, otherwise every bracket would be
            ## drawn in every panel
            pairwise.tb$prot.id = factor(pairwise.tb$prot.id, levels = protein.id)

            bracket.columns = c(switch(multiple.proteins + 1, NULL, "prot.id"),
                                "group1", "group2", "y.position", "p.label")

            boxplot =
              boxplot +
              ggpubr::stat_pvalue_manual(data = pairwise.tb[,bracket.columns,drop=FALSE],
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
