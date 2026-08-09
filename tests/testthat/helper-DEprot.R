## ----------------------------------------------------------------------------------------
##  Shared fixtures
##  This file is sourced by testthat before every test file: the objects and the helpers
##  defined here are available everywhere, which avoids rebuilding them in each file.
## ----------------------------------------------------------------------------------------

## Real example objects
tb.dpo.raw    <- DEprot::test.toolbox$dpo.raw
tb.dpo.norm   <- DEprot::test.toolbox$dpo.norm
tb.dpo.random <- DEprot::test.toolbox$dpo.random
tb.dpo.imp    <- DEprot::test.toolbox$dpo.imp
tb.limma      <- DEprot::test.toolbox$diff.exp.limma


#' Temporary directory removed when the calling test exits
#'
#' 'withr' is not among the dependencies of the package, hence the deferred cleanup is
#' registered on the environment of the caller through testthat.

local.tmpdir <-
  function(envir = parent.frame()) {
    path <- file.path(tempdir(), paste0("DEprot.test.", paste0(sample(letters, 8), collapse = "")))
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    withr::defer(unlink(path, recursive = TRUE, force = TRUE), envir = envir)
    return(path)
  }


#' Small synthetic count matrix carrying missing values
#'
#' @param n.prot Number of proteins.
#' @param n.samples Number of samples (split evenly between two conditions).
#' @param n.missing Number of values replaced by NA.
#' @param seed Seed.

make.counts <-
  function(n.prot = 60,
           n.samples = 6,
           n.missing = 25,
           seed = 7) {

    set.seed(seed)

    m <- matrix(stats::rnorm(n.prot * n.samples, mean = 22, sd = 2),
                nrow = n.prot,
                dimnames = list(paste0("prot.", seq_len(n.prot)),
                                paste0("s", seq_len(n.samples))))

    ## the missing values are kept away from the first rows, so that at least a few
    ## proteins stay complete whatever the imputation method being tested
    if (n.missing > 0) {
      idx <- sample(seq(from = length(m) %/% 4, to = length(m)), n.missing)
      m[idx] <- NA
    }

    return(m)
  }


#' Minimal metadata matching make.counts()

make.metadata <-
  function(m) {
    data.frame(column.id = colnames(m),
               condition = rep(c("ctrl", "treat"), length.out = ncol(m)),
               replicate = rep(paste0("rep", seq_len(ncol(m) / 2)), each = 2)[seq_len(ncol(m))],
               stringsAsFactors = FALSE)
  }


#' Synthetic DEprot object built on make.counts()
#'
#' Counts are stored in the 'raw', 'normalized' and 'randomized' slots, so that the object
#' can be handed to any function whatever the 'which.data' being tested.

make.dpo <-
  function(n.prot = 60,
           n.samples = 6,
           n.missing = 25,
           seed = 7) {

    m <- make.counts(n.prot = n.prot, n.samples = n.samples, n.missing = n.missing, seed = seed)

    obj <- DEprot::load.counts2(counts = m,
                                metadata = make.metadata(m),
                                data.type = "raw",
                                log.base = 2)

    obj@norm.counts <- m
    obj@normalized <- TRUE
    obj@normalization.method <- "test"
    obj@random.counts <- m
    obj@randomized <- TRUE

    return(obj)
  }


## ----------------------------------------------------------------------------------------
##  Fixture writers: minimal reports mimicking the output of each search engine.
##  They are deliberately tiny (a handful of proteins), the point being to walk through the
##  readers rather than to test the quantification itself.
## ----------------------------------------------------------------------------------------

#' DIA-NN main report (long format)

write.diann.report <-
  function(dir,
           file = "report.tsv") {

    runs <- c("D:/data/20240101_sampleA.raw", "D:/data/20240101_sampleB.raw")

    tb <- expand.grid(Run = runs,
                      Protein.Group = c("P00001", "P00002", "P00003", "CON__P99999"),
                      Precursor.Id = c("PEPTIDEK2", "PEPTIDER2"),
                      stringsAsFactors = FALSE)

    tb$Protein.Ids <- tb$Protein.Group
    tb$Genes <- paste0("GENE", substr(tb$Protein.Group, 6, 6))
    tb$Modified.Sequence <- sub("[0-9]$", "", tb$Precursor.Id)
    tb$Precursor.Charge <- 2
    tb$Q.Value <- 0.001
    tb$PG.Q.Value <- 0.001
    tb$Precursor.Normalised <- seq_len(nrow(tb)) * 1000
    tb$PG.MaxLFQ <- seq_len(nrow(tb)) * 2000

    path <- file.path(dir, file)
    utils::write.table(tb, path, sep = "\t", row.names = FALSE, quote = FALSE)
    return(path)
  }


#' DIA-NN pg_matrix (wide format)

write.diann.matrix <-
  function(dir,
           file = "report.pg_matrix.tsv") {

    tb <- data.frame(Protein.Group = c("P00001", "P00002", "P00003"),
                     Protein.Ids = c("P00001", "P00002", "P00003"),
                     Genes = c("GENE1", "GENE2", "GENE3"),
                     First.Protein.Description = "desc",
                     `D:/data/sampleA.raw` = c(1000, 2000, 3000),
                     `D:/data/sampleB.raw` = c(1500, 2500, 3500),
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

    path <- file.path(dir, file)
    utils::write.table(tb, path, sep = "\t", row.names = FALSE, quote = FALSE)
    return(path)
  }


#' MaxQuant proteinGroups.txt

write.maxquant <-
  function(dir,
           file = "proteinGroups.txt") {

    tb <- data.frame(`Majority protein IDs` = c("P00001", "P00002", "P00003", "P00004"),
                     `Protein IDs` = c("P00001", "P00002", "P00003", "P00004"),
                     `Gene names` = c("GENE1", "GENE2", "GENE3", "GENE4"),
                     Reverse = c("", "", "+", ""),
                     `Potential contaminant` = c("", "", "", "+"),
                     `LFQ intensity sampleA` = c(1000, 2000, 3000, 4000),
                     `LFQ intensity sampleB` = c(1100, 2100, 3100, 4100),
                     `Intensity sampleA` = c(10, 20, 30, 40),
                     `Intensity sampleB` = c(11, 21, 31, 41),
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

    path <- file.path(dir, file)
    utils::write.table(tb, path, sep = "\t", row.names = FALSE, quote = FALSE)
    return(path)
  }


#' Spectronaut pivot report (wide format)

write.spectronaut.pivot <-
  function(dir,
           file = "spectronaut.pivot.tsv") {

    tb <- data.frame(PG.ProteinGroups = c("P00001", "P00002", "P00003"),
                     PG.Genes = c("GENE1", "GENE2", "GENE3"),
                     `[1] sampleA.raw.PG.Quantity` = c(1000, 2000, 3000),
                     `[2] sampleB.raw.PG.Quantity` = c(1200, 2200, 3200),
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

    path <- file.path(dir, file)
    utils::write.table(tb, path, sep = "\t", row.names = FALSE, quote = FALSE)
    return(path)
  }


#' Generic wide table, without any recognizable signature

write.generic <-
  function(dir,
           file = "generic.tsv") {

    tb <- data.frame(id = c("P00001", "P00002", "P00003"),
                     ctrl.1 = c(1000, 2000, 3000),
                     ctrl.2 = c(1100, 2100, 3100),
                     treat.1 = c(2000, 4000, 6000),
                     treat.2 = c(2100, 4100, 6100),
                     stringsAsFactors = FALSE)

    path <- file.path(dir, file)
    utils::write.table(tb, path, sep = "\t", row.names = FALSE, quote = FALSE)
    return(path)
  }


#' ProteinLevelData table as returned by MSstats::dataProcess()

make.msstats.lfq <-
  function(n.prot = 5,
           runs = c("run.A_1", "run.B_1", "run.A_2", "run.B_2")) {

    tb <- expand.grid(Protein = paste0("P0000", seq_len(n.prot)),
                      originalRUN = runs,
                      stringsAsFactors = FALSE)

    tb$RUN <- match(tb$originalRUN, runs)
    tb$GROUP <- ifelse(grepl("run.A", tb$originalRUN), "ctrl", "treat")
    tb$SUBJECT <- sub(".*_", "rep", tb$originalRUN)
    tb$LogIntensities <- 20 + seq_len(nrow(tb)) / 10
    tb$NumMeasuredFeature <- 3

    return(tb)
  }
