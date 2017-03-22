#' Summary for Posterior Model Probabilities
#'
#' Summary for the estiamted posterior model probabilities (\code{\link{stationary}}), estimated effective sample size, and summaries for all pairwise Bayes factors.
#'
#' @param object posterior samples of the stationary distribution (rows = replications; columns = model probabilities)
#' @param BF whether to compute summaries for all pairwise Bayes factors
#' @param logBF whether to summarize log Bayes factors instead of Bayes factors
#' @param ... ignored
#'
#' @examples
#' pp <- matrix(runif(100*3),100)
#' pp <- pp/rowSums(pp)
#' class(pp) <- "stationary"
#' summary(pp)
#'
#' @return a list with estimates for \code{"pp"} = model posterior probabilities, \code{"n.eff"} = effective sample size, and \code{"bf"} = pairwise Bayes factors (optional)
#' @seealso \code{\link{stationary}}, \code{\link{dirichlet.mle}}
#' @export
summary.stationary <- function(object,
                               BF = FALSE,
                               logBF = FALSE,
                               ...){
  samples <- object
  pp <- t(apply(samples, 2, summ.samples))

  n.eff <- NA
  try (n.eff <- floor(dirichlet.mle(na.omit(samples))$sum))

  res <- list(pp=pp, n.eff=n.eff)

  if (BF){
    combs <- combn(ncol(samples),2)

    if (!is.null(colnames(samples))){
      rownames(pp) <- colnames(samples)
      c1 <- rownames(pp)[combs[1,]]
      c2 <- rownames(pp)[combs[2,]]
      bf.names <- paste0("BF_", apply(rbind(c1,c2),2,paste,collapse=""))

    } else {
      rownames(pp) <- paste0("M",1:ncol(samples))
      bf.names <- paste0("BF_", apply(combs,2,paste,collapse=""))
    }

    res$bf <- matrix(NA, ncol(combs), 5,
                     dimnames = list(BF = bf.names,
                                     Statistic = colnames(pp)))
    for (bb in 1:ncol(combs)){
      suppressWarnings(
        bf.tmp <- log(samples[,combs[1,bb]])-log(samples[,combs[2,bb]])
      )
      if (!logBF)
        bf.tmp <- exp(bf.tmp)
      res$bf[bb,] <- summ.samples(bf.tmp)
    }
  }
  return (res)
}

summ.samples <- function(x)
  c(Mean=mean(x, na.rm = TRUE),
    SD=sd(x, na.rm = TRUE),
    quantile(x, probs=c(.05,.5,.95), na.rm = TRUE))
