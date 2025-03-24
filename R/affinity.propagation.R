#' @title affinity.propagation
#'
#' @description
#' Use affinity propagation to cluster similar gene sets to reduce redundancy in report.
#'
#' @param idsInSet A list of set names and their member IDs.
#' @param score A vector of addible scores with the same length used to assign input preference;
#'  higher score has larger weight, i.e. -logP.
#'
#' @return A list of \code{clusters} and \code{representatives} for each cluster.
#' \describe{
#'  \item{clusters}{A list of character vectors of set IDs in each cluster.}
#'  \item{representatives}{A character vector of representatives for each cluster.}
#' }
#'
#' @export affinity.propagation
#' @importFrom apcluster apcluster
#' @importFrom stats rnorm
#' @author Zhiao Shi, Yuxing Liao

affinity.propagation <- function(idsInSet, score) {
  ## Libraries
  # require(apcluster)

  cat("Begin affinity propagation...\n")
  # compute the similiarity and input preference vector
  ret <- jaccardSim(idsInSet, as.numeric(score))

  sim.mat <- ret$sim.mat
  ip.vec <- ret$ip.vec

  apRes <- apcluster::apcluster(sim.mat,p=ip.vec)
  #sort clusters to make exemplar the first member
  clusters <- vector(mode="list", length(apRes@clusters))
  if(length(apRes@clusters) == 0){
    return(NULL)
  }
  for (i in 1:length(apRes@clusters)) {
    exemplar <- apRes@exemplars[[i]]
    clusters[[i]] <- apRes@clusters[[i]][order(apRes@clusters[[i]] == exemplar, decreasing=TRUE)]
  }
  cat("End affinity propagation...\n")
  return(list(clusters=sapply(clusters, names), representatives=names(apRes@exemplars)))
}

