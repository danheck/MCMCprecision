
#' Estimate stationary distribution of discrete MCMC output
#'
#' Transdimensional MCMC methods include a discrete model-indicator variable \eqn{Z} with a fixed but unknown stationary distribution \eqn{\pi} (i.e., the model posterior probabiltiies). This function provides posterior samples quantiying the estimation uncertainty of \eqn{\pi}.
#'
#' @param z MCMC output for the discrete indicator variable (can be a \code{\link[coda]{mcmc.list}} and include numerical or character labels)
#' @param labels a-prior vector of labels for all models to include models not sampled in chain \code{z}. This does not affect inferences due to the improper Dirichlet(0,..,0) prior.
#' @param add.diag constant added to the diagonal to ensure that the observed transition matrix is not singular. The special value \code{add.diag="1/M"} adds 1/M with M=number of models.
#' @param sample number of samples
#' @param summary whether to summarize results and get effective sample size (otherwise, the raw samples are r' @parameturned)
#' @param logBF whether to summarize log(BF) instead of BF
#' @param cpu number of CPUs used. Will only speed up computations for large numbers of models.
#' @param method how to compute eigenvectors:
#' \itemize{
#'   \item \code{"cpp"}: This default uses \code{RcppArmadillo::eig_gen} and provides a compromise of speed and stability.
#'   \item \code{"cpps"}: Uses sparse matrices in \code{RcppArmadillo::eigs_gen} and can be faster (but also more unstable) for large number of models.
#'   \item \code{"base"}: Uses \code{base::\link[base]{eigen}}, which is slower but most stable.
#   \item \code{"sparse"} (=> rARPACK:eigs; uses sparse matrices)
#'   \item \code{"iid"}: Assume independent sampling of model indicator variable and is only implemented as a benchmark. Results cannot be trusted if the samples \code{z} are correlated (which is usually the case of transdimensional MCMC output)
#'  }
#'
#' @details The mathod draws independent posterior samples of the transition matrix P for the discrete-valued indicator variable Z (which indexes the model). For each row of the transition matrix, a Dirichlet(0,...,0) prior is assumed, resulting in a conjugate Dirichlet posterior. For each sample, the eigenvector with eigenvalue 1 is computed and normalized to sum up to one. The resulting samples can be used to assess the estimation uncertainty in the true stationary distribution of interest (e.g., the model posterior probabilities).
#' @return posterior samples for the model posterior probabilities (i.e., the fixed but unknown stationary distribution of \code{z}). If \code{summary=TRUE}, a \code{\link{summary.stationary}} is returned.
#' @examples
#' P <- LaplacesDemon::rdirichlet(3, c(2,1,.5))
#' z <- sim.mc(1000, P)
#' stationary(z)
#' @export
stationary <- function(z,
                       labels=NULL,
                       add.diag=0,
                       sample=1000,
                       method ="cpp",
                       cpu=1,
                       summary=FALSE,
                       logBF=FALSE){
  tab <- table.mc(z, labels=labels)
  labels <- rownames(tab)
  M <- nrow(tab)
  if(det(tab) == 0 & add.diag != 0){
    if(add.diag == "1/M") add.diag <- 1/M
    diag(tab) <- diag(tab) + add.diag
    warning("Matrix 'tab' is singular! add.diag=", round(add.diag, 5),
            " is added to the diagonal.")
  }

  if(method == "cpp"){
    samp <- stationaryCpp(tab, sample)
  }else if(method == "cpps"){
    samp <- stationaryCppSparse(Matrix(tab, sparse = TRUE), sample)
  }else if(method == "iid"){
    samp <- rdirichlet(sample, rowSums(tab))
  }else if(cpu == 1){
    samp <- t(sapply(1:sample, posterior.sample,
                     tab=tab, method=method))
  }else{
    cl <- makeCluster(cpu)
    tmp <- clusterEvalQ(cl, {library(Matrix)}) #; library(rARPACK)})
    clusterExport(cl, c("posterior.sample"))
    samp <- t(parSapply(cl, 1:sample, posterior.sample,
                        tab=tab, method=method))
    stopCluster(cl)
  }
  # stability: rerun with base
  sel <- apply(is.na(samp) | samp == Inf, 1, any)
  if(sum(sel) > 0)
    samp[sel,] <-  t(sapply(1:sum(sel), posterior.sample,
                            tab=tab, method="base"))

  colnames(samp) <- colnames(tab)
  if(summary){
    summ <- summary.stationary(samp, logBF =logBF)
    if(method == "iid")
      summ$neff <- sum(tab)+1
    return(summ)
  }else{
    class(samp) <- "stationary"
    return(samp)
  }
}


# #' Print summary for posterior samples of stationary distribution
# #'
# #' See \code{\link{summary.stationary}}.
# #'
# #' @param ... ignored
# #' @param digits digits for rounding
# # #' @describeIn summary print summary for posterior samples
#' @export
print.stationary <- function(x, digits=NULL, ...){
  if(!is.null(digits))
    res <- lapply(summary(x), round, digits=digits)
  else
    res <- summary(x)
  print(res,...)
}

