#' Estimate stationary distribution for discrete MCMC variables
#'
#' Transdimensional MCMC methods include a discrete model-indicator variable \eqn{Z} with a fixed but unknown stationary distribution \eqn{\pi} (i.e., the model posterior probabiltiies). This function provides posterior samples quantiying the estimation uncertainty of \eqn{\pi}.
#'
#' @param z MCMC output for the discrete indicator variable with numerical, character, or factor labels (can also be a \code{\link[coda]{mcmc.list}} or a matrix with one MCMC chain per column).
#' @param N the observed transition matrix (if supplied, \code{z} is ignored). A quadratic matrix with switching frequencies \eqn{N[i,j] = #{ z[t]=i & z[t+1]=j }}
#' @param labels optional: vector of labels for complete set of models (e.g., models not sampled in the chain \code{z}). If \code{epsilon=0}, this does not affect inferences due to the improper Dirichlet(0,..,0) prior.
#'
#' @param sample number of samples to be drawn from the posterior of the stationary distribution \eqn{\pi}
#' @param N.min minimum frequency for sampled models to be included in the transition matrix \eqn{N}. The default \code{N.min = 1} includes all sampled models, which can lead to numerical issues when computing eigenvectors (i.e., all but one of the posterior model probabilities are estimated to be zero).
#' @param epsilon prior parameter for the rows of the estimated transition matrix P, that is, P[i,] ~ Dirichlet\eqn{(\epsilon, ..., \epsilon)}. The default \code{epsilon=0} minimizes the impact of the prior and also renders non-sampled models irrelevant. If \code{method="iid"}, the Dirichlet prior is assumed directly on the stationary distribution \eqn{\pi} instead of the rows of the transition matrix.
#' @param summary whether the output should be a summary or posterior samples
# ' @param logBF whether to summarize log(BF) instead of BF (see \code{\link{summary.stationary}})
#' @param cpu number of processing units used for parallel sampling. Will only speed up computations for large numbers of models.
#' @param progress whether to show a progress bar (not functional for \code{cpu>1} and \code{method="base"})
#'
#' @param method function used to compute eigenvectors:
#' \itemize{
#'   \item \code{"base"}: Uses \code{base::\link[base]{eigen}}, which is most stable, but can be slower than \code{"arma"}
#'   \item \code{"arma"}: Uses \code{RcppArmadillo::eig_gen}
#'   \item \code{"armas"}: Uses sparse matrices and \code{RcppArmadillo::eigs_gen}, which can be faster for very large number of models (but also numerically unstable)
#'   \item \code{"eigen"}: Uses package\code{RcppEigen}
#'   \item \code{"iid"}: Assume i.i.d. sampling of the model indicator variable \code{z}. This is only implemented as a benchmark, because results cannot be trusted if the samples \code{z} are correlated (which is usually the case for transdimensional MCMC output)
#'  }
#' @param digits number of digits that are used for checking whether the first eigenvalue is equal to 1 (any difference must be due to numerical precision)
#'
#' @details The method draws independent posterior samples of the transition matrix \eqn{P} for the discrete-valued indicator variable \code{z} (usually, a sequence of sampled models). For each row of the transition matrix, a Dirichlet(0,...,0) prior is assumed, resulting in a conjugate Dirichlet posterior. For each sample, the eigenvector with eigenvalue 1 is computed and normalized (i.e., devided by its sum). The (independent) posterior samples can be used to assess the estimation uncertainty in the 'true' stationary distribution of interest (e.g., the model posterior probabilities).
#'
#' @return default: a summary for the posterior distribution of the model posterior probabilities (i.e., the fixed but unknown stationary distribution of \code{z}) and Bayes factors. otherwise, if \code{summary=FALSE}: posterior samples.
#' @seealso \code{\link{best.k}}, \code{\link{summary.stationary}}
#'
#' @examples
#' P <- matrix(c(.9,.1,0,
#'               .1,.6,.3,
#'               .2,.3,.5), 3, byrow=TRUE)
#' z <- sim.mc(1000, P)
#' stationary(z)
#'
#' # input: transition frequency
#' tab <- table.mc(z)
#' stationary(N = tab)
#' @export
stationary <- function (z,
                        N,
                        labels,
                        sample = 1000,
                        N.min = 1,
                        epsilon = 0,
                        cpu = 1,
                        method = "base",
                        digits = 6,
                        progress = TRUE,
                        summary = TRUE){
  if (epsilon <0)
    stop ("'epsilon' must be zero or positive.")
  if (missing(labels))
    labels <- NULL

  if (!missing(N) && !is.null(N)){
    if (ncol(N) != nrow(N) || any(N<0) || any(N != round(N)))
      stop ("The transition matrix 'N' has negative or non-integer values.")
    tab <- as.matrix(N)
  } else {
    tab <- table.mc(z, labels=labels)
  }

  # numerical issues for large condition number
  # (i.e., when maximum and minimum frequency in N differ a lot)
  if (N.min  > 1){
    sel.tab <- colSums(tab) >= N.min
    tab <- tab[sel.tab, sel.tab, drop = FALSE]
  }

  labels <- rownames(tab)
  if (cpu == 1){
    if (method == "arma"){
      samp <- stationaryCpp(tab, sample = sample, epsilon = epsilon,
                            digits = digits, display_progress=progress)

    } else if (method == "armas") {
      if (epsilon != 0)
        warning ("'epsilon is ignored if method='armas' (sparse matrices) is used.")
      samp <- stationaryCppSparse(Matrix(tab, sparse = TRUE), sample = sample,
                                  digits = digits, display_progress=progress)

    } else if (method == "iid"){
      samp <- rdirichlet(sample, rowSums(tab) + epsilon)

    } else if (method == "base"){
      samp <- matrix(NA, sample, nrow(tab))
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
        warning ("'epsilon is ignored if method='armas' (sparse matrices) is used.")
      samp <- do.call("rbind",
                      parSapply(cl, X = rep(ceiling(sample/cpu), cpu), simplify = FALSE,
                                FUN = function (ss)
                                  stationaryCppSparse(Matrix(tab, sparse = TRUE), sample = ss,
                                                      digits = digits, display_progress=FALSE)))
    } else if (method == "iid"){
      samp <-  do.call("rbind", parSapply(cl = cl, X = rep(ceiling(sample/cpu), cpu),
                                          FUN = rdirichlet, a = rowSums(tab) + epsilon,
                                          simplify = FALSE))
    } else if (method == "base"){
      samp <- t(parSapply(cl, 1:sample, posterior.sample,
                          tab=tab, digits = digits))
    } else {
      stop ("method not supported.")
    }
    stopCluster(cl)
  }

  samp[samp < 0] <- 0
  if (anyNA(samp))
    warning ("Sampled posterior model probabilities contain missing values.\n",
             "  This might be due to a low precision of the eigenvalue decomposition.")
  if (kappa(tab) == Inf)
    warning ("Condition number of transition matrix is infinite.\n",
             "  The eigenvalue decomposition might be numerically unstable.\n",
             "  Use 'N.min=2' (or larger) to check for robustness.")
  if (!all(is.na(samp)) && (max(samp, na.rm = TRUE) == 1 ))
    warning ("Some models are assigned posterior probability of one.\n",
             "  This indicates numerical issues with the eigenvalue decomposition.")

  colnames(samp) <- colnames(tab)
  if (summary){
    summ <- summary.stationary(samp)
    if (method == "iid")
      summ$n.eff <- sum(tab) + 1
    return (summ)
  } else {
    class(samp) <- c("stationary", "matrix")
    return (samp)
  }
}

