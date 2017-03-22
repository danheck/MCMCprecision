#' Estimate stationary distribution for discrete MCMC variables
#'
#' Transdimensional MCMC methods include a discrete model-indicator variable \eqn{Z} with a fixed but unknown stationary distribution \eqn{\pi} (i.e., the model posterior probabiltiies). The function \code{stationary} draws posterior samples to assess the estimation uncertainty for \eqn{\pi}.
#'
#' @param z MCMC output for the discrete indicator variable with numerical, character, or factor labels (can also be a \code{\link[coda]{mcmc.list}} or a matrix with one MCMC chain per column).
#' @param N the observed transition matrix (if supplied, \code{z} is ignored). A quadratic matrix with sampled transition frequencies \eqn{N[i,j] = #{ z[t]=i & z[t+1]=j }}
#' @param labels optional: vector of labels for complete set of models (e.g., models not sampled in the chain \code{z}). If \code{epsilon=0}, this does not affect inferences due to the improper Dirichlet(0,..,0) prior.
#'
#' @param sample number of samples to be drawn from the posterior of the stationary distribution \eqn{\pi}
#' @param epsilon prior parameter for the rows of the estimated transition matrix \eqn{P}: \eqn{P[i,]} ~ Dirichlet\eqn{(\epsilon, ..., \epsilon)}. The default \code{epsilon="1/M"} (with M = number of sampled models) provides results close to the i.i.d. estimates and is numerically stable. The alternative \code{epsilon=0} minimizes the impact of the prior and also renders non-sampled models irrelevant. If \code{method="iid"}, a Dirichlet prior is assumed on the stationary distribution \eqn{\pi} instead of the rows of the transition matrix \eqn{P}.
#' @param summary whether the output should be summarized. If \code{FALSE}, posterior samples are returned
#' @param cpu number of CPUs used for parallel sampling. Will only speed up computations for large numbers of models (i.e., for large transition matrices).
#' @param progress whether to show a progress bar (not functional for \code{cpu>1})
#'
#' @param method how to compute eigenvectors:
#' \itemize{
#'   \item \code{"base"}: Uses \code{base::\link[base]{eigen}}, which is most stable, but can be slower than \code{"arma"} for small transition matrices
#'   \item \code{"arma"}: Uses \code{RcppArmadillo::eig_gen}
#'   \item \code{"armas"}: Uses sparse matrices with \code{RcppArmadillo::eigs_gen}, which can be faster for very large number of models if \code{epsilon=0} (but also numerically unstable)
#'   \item \code{"eigen"}: Uses package \code{RcppEigen::EigenSolver}
#'   \item \code{"iid"}: Assumes i.i.d. sampling of the model indicator variable \code{z}. This is only implemented as a benchmark, because results cannot be trusted if the samples \code{z} are correlated (which is usually the case for transdimensional MCMC output)
#'  }
#' @param digits number of digits that are used for checking whether the first eigenvalue is equal to 1 (any difference must be due to low numerical precision)
#'
#' @details The method draws independent posterior samples of the transition matrix \eqn{P} for the discrete-valued indicator variable \code{z} (usually, a sequence of sampled models). For each row of the transition matrix, a Dirichlet\eqn{(\epsilon,...,\epsilon)} prior is assumed, resulting in a conjugate Dirichlet posterior. For each sample, the eigenvector with eigenvalue 1 is computed and normalized. The (independent) posterior samples can be used to assess the estimation uncertainty in the 'true' stationary distribution of interest (e.g., the model posterior probabilities) and to estimate the effective sample size (see \code{\link{summary.stationary}}).
#'
#' @return default: a summary for the posterior distribution of the model posterior probabilities (i.e., the fixed but unknown stationary distribution of \code{z}) and Bayes factors. If \code{summary=FALSE}: posterior samples.
#' @seealso \code{\link{best.k}}, \code{\link{summary.stationary}}
#'
#' @examples
#' # data-generating transition matrix
#' P <- matrix(c(.1,.5,.4,
#'               0, .5,.5,
#'               .9,.1,0), ncol = 3, byrow=TRUE)
#'
#' # input: sequence of sampled models
#' z <- sim.mc(500, P)
#' stationary(z)
#'
#' # input: transition frequency
#' N <- table.mc(z)
#' samples <- stationary(N = N, summary = FALSE)
#'
#' # summaries:
#' best.k(samples, k = 3)
#' summary(samples)
#' @export
stationary <- function (z,
                        N,
                        labels,
                        sample = 1000,
                        epsilon = "1/M",
                        cpu = 1,
                        method = "base",
                        digits = 6,
                        progress = TRUE,
                        summary = TRUE){
  if (missing(labels))
    labels <- NULL

  if (!missing(N) && !is.null(N)){
    if (ncol(N) != nrow(N) || any(N<0))
      stop ("The transition matrix 'N' has negative values.")
    tab <- as.matrix(N)
  } else {
    tab <- table.mc(z, labels=labels)
  }

  M <- nrow(tab)
  if (epsilon == "1/M"){
    epsilon <- 1/M
  } else if (epsilon <0){
    stop ("'epsilon' must be zero or positive.")
  }

  labels <- rownames(tab)
  if (cpu == 1){
    if (method == "arma"){
      samp <- stationaryCpp(tab, sample = sample, epsilon = epsilon,
                            digits = digits, display_progress=progress)

    } else if (method == "armas") {
      if (epsilon != 0)
        warning ("'epsilon=0' is used if method='armas' (sparse matrices)")
      samp <- stationaryCppSparse(Matrix(tab, sparse = TRUE), sample = sample,
                                  digits = digits, display_progress=progress)

    } else if (method == "iid"){
      samp <- rdirichlet(sample, rowSums(tab) + epsilon)

    } else if (method == "base"){
      samp <- matrix(NA, sample, M)
      if (progress)
        prog <- txtProgressBar(0, sample, style=3, width=50)
      for (i in 1:sample){
        if (progress) setTxtProgressBar(prog, i)
        samp[i,] <- posterior.sample(1, tab=tab, epsilon=epsilon,
                                     digits = digits)
      }
      if (progress) close(prog)

    } else if (method == "eigen") {
      samp <- stationaryEigen(tab, sample = sample, epsilon = epsilon,
                              digits = digits, display_progress=progress)
      samp[samp == -99] <- NA

    } else {
      stop ("method not supported.")
    }

  } else {
    ################################ MULTICORE SUPPORT
    cl <- makeCluster(cpu)
    if (method == "arma"){
      samp <- do.call("rbind",
                      parSapply(cl, X = rep(ceiling(sample/cpu), cpu), simplify = FALSE,
                                FUN = function (ss)
                                  stationaryCpp(tab, sample = ss, epsilon = epsilon,
                                                digits = digits, display_progress = FALSE)))
    } else if (method == "armas"){
      if (epsilon != 0)
        warning ("'epsilon=0' is used if method='armas' (sparse matrices)")
      samp <- do.call("rbind",
                      parSapply(cl, X = rep(ceiling(sample/cpu), cpu), simplify = FALSE,
                                FUN = function (ss)
                                  stationaryCppSparse(Matrix(tab, sparse = TRUE), sample = ss,
                                                      digits = digits, display_progress=FALSE)))
    } else if (method == "eigen"){
      samp <- do.call("rbind",
                      parSapply(cl, X = rep(ceiling(sample/cpu), cpu), simplify = FALSE,
                                FUN = function (ss)
                                  stationaryEigen(tab, sample = ss, epsilon = epsilon,
                                                  digits = digits, display_progress=FALSE)))
    } else if (method == "iid"){
      samp <-  do.call("rbind", parSapply(cl = cl, X = rep(ceiling(sample/cpu), cpu),
                                          FUN = rdirichlet, a = rowSums(tab) + epsilon,
                                          simplify = FALSE))
    } else if (method == "base"){
      samp <- t(parSapply(cl, 1:sample, posterior.sample, epsilon = epsilon,
                          tab=tab, digits = digits))
    } else {
      stop ("method not supported.")
    }
    stopCluster(cl)
  }

  samp[samp < 0] <- 0
  if (anyNA(samp))
    warning (
      "Sampled posterior model probabilities contain NAs/NaNs.\n",
      "  This might be due to low precision of the eigenvalue decomposition.\n",
      "  As a remedy, the Dirichlet prior can be changed, e.g., 'epsilon = .01'.")

  if (!all(is.na(samp)) && (max(samp, na.rm = TRUE) == 1 ))
    warning (
      "Some models are assigned posterior probability of one.\n",
      "  This indicates numerical issues with the eigenvalue decomposition.\n",
      "  As a remedy, the Dirichlet prior can be changed, e.g., 'epsilon = .01'.")

  colnames(samp) <- colnames(tab)
  class(samp) <- c("stationary", "matrix")
  attr(samp, "epsilon") <- epsilon
  attr(samp, "method") <- method

  if (summary){
    samp <- summary.stationary(samp)
    if (method == "iid")
      samp$n.eff <- sum(tab) + epsilon*M
  }
  return (samp)
}

