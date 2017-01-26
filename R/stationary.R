#' Estimate stationary distribution of discrete MCMC output
#'
#' Transdimensional MCMC methods include a discrete model-indicator variable \eqn{Z} with a fixed but unknown stationary distribution \eqn{\pi} (i.e., the model posterior probabiltiies). This function provides posterior samples quantiying the estimation uncertainty of \eqn{\pi}.
#'
#' @param z MCMC output for the discrete indicator variable with numerical, character, or factor labels (can be a \code{\link[coda]{mcmc.list}} or a matrix with one chain per column). If \code{z} is a quadratic matrix, it will be interpreted as the observed transition matrix.
#' @param labels a-prior vector of labels for all models to include models not sampled in chain \code{z}. This does not affect inferences due to the improper Dirichlet(0,..,0) prior.
# DEPRECATED: @param add.diag constant added to the diagonal to ensure that the observed transition matrix is not singular. The special value \code{add.diag="1/M"} adds 1/M with M=number of models.
#' @param sample number of samples
#' @param summary whether to summarize results and get effective sample size (otherwise, the raw samples are r' @parameturned)
#' @param logBF whether to summarize log(BF) instead of BF
# @param cpu number of CPUs used. Will only speed up computations for large numbers of models.
#' @param progress whether to show the progress (not for \code{cpu>4} and \code{method="base"})
#' @param method how to compute eigenvectors:
#' \itemize{
#'   \item \code{"cpp"}: This default uses \code{RcppArmadillo::eig_gen} and provides a compromise of speed and stability.
#'   \item \code{"cpps"}: Uses sparse matrices in \code{RcppArmadillo::eigs_gen} and can be faster for large number of models (but also more unstable).
#'   \item \code{"base"}: Uses \code{base::\link[base]{eigen}}, which is slower but might be more stable.
#   \item \code{"sparse"} (=> rARPACK:eigs; uses sparse matrices)
#'   \item \code{"iid"}: Assume independent sampling of model indicator variable and is only implemented as a benchmark. Results cannot be trusted if the samples \code{z} are correlated (which is usually the case of transdimensional MCMC output)
#'  }
#'
#' @details The mathod draws independent posterior samples of the transition matrix P for the discrete-valued indicator variable Z (which indexes the model). For each row of the transition matrix, a Dirichlet(0,...,0) prior is assumed, resulting in a conjugate Dirichlet posterior. For each sample, the eigenvector with eigenvalue 1 is computed and normalized to sum up to one. The resulting samples can be used to assess the estimation uncertainty in the true stationary distribution of interest (e.g., the model posterior probabilities).
#' @return By default, a summary for the posterior distribution of the model posterior probabilities (i.e., the fixed but unknown stationary distribution of \code{z}) and Bayes factors is returned. Use \code{summary=FALSE} to get the raw samples.
#' @examples
#' P <- matrix(c(.9,.1,0,
#'               .1,.6,.3,
#'               .2,.3,.5), 3, byrow=TRUE)
#' z <- sim.mc(1000, P)
#' stationary(z)
#' @export
stationary <- function(z,
                       labels=NULL,
                       sample=1000,
                       # cpu=4,
                       method ="cpp",
                       progress=TRUE,
                       summary=TRUE,
                       logBF=FALSE){
  if(is.matrix(z) && ncol(z) == nrow(z)){
    if(any(z<0) || any(z != round(z)))
      stop("The transition matrix 'z' has negative or non-integer values.")
    tab <- z
  }
  else{
    tab <- table.mc(z, labels=labels)
  }
  labels <- rownames(tab)
  M <- nrow(tab)

  if(method == "cpp"){
    samp <- stationaryCpp(tab, sample = sample,
                          display_progress=progress)
  }else if(method == "cpps"){
    samp <- stationaryCppSparse(Matrix(tab, sparse = TRUE), sample = sample,
                                display_progress=progress)
  }else if(method == "iid"){
    samp <- rdirichlet(sample, rowSums(tab))
  }else if(method == "base"){
    samp <- matrix(NA, sample, M)
    if(progress)
      prog <- txtProgressBar(0, sample, style=3, width=50)
    for(i in 1:sample){
      if(progress) setTxtProgressBar(prog, i)
      samp[i,] <- posterior.sample(1, tab=tab, method=method)
    }
    if(progress) close(prog)
  }else{
    stop("method not supported.")
    # cl <- makeCluster(cpu)
    # tmp <- clusterEvalQ(cl, {library(Matrix)}) #; library(rARPACK)})
    # clusterExport(cl, c("posterior.sample"),
    #               envir = environment())
    # samp <- t(parSapply(cl, 1:sample, posterior.sample,
    #                     tab=tab, method=method))
    # stopCluster(cl)
  }

  colnames(samp) <- colnames(tab)
  if(summary){
    summ <- summary.stationary(samp, logBF =logBF)
    if(method == "iid")
      summ$neff <- sum(tab)+1
    return(summ)
  }else{
    class(samp) <- c("stationary", "matrix")
    return(samp)
  }
}

