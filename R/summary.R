
#' Summarize posterior probabilities and BF estimates
#' @param samples posterior samples of the stationary distribution (rows = replications; columns = model probabilities)
#' @param logBF whether to summarize log Bayes factors instead of Bayes factors
#' @examples
#' n <- 100
#' pp <- LaplacesDemon::rdirichlet(200, n*c(.7,.2,.1))
#' class(pp) <- "stationary"
#' summary(pp)
#' @return a list with estimates for \code{"pp"} = model posterior probabilities, \code{"bf"} = Bayes factors, \code{"neff"} = effective sample size
#' @export
summary.stationary <- function(samples,
                               logBF=FALSE){
  pp <- t(apply(samples, 2, summ.samples))
  combs <- combn(ncol(samples),2)

  if(!is.null(colnames(samples))){
    rownames(pp) <- colnames(samples)
    c1 <- rownames(pp)[combs[1,]]
    c2 <- rownames(pp)[combs[2,]]
    bf.names <- paste0("BF_", apply(rbind(c1,c2),2,paste,collapse=""))
  }else{
    rownames(pp) <- paste0("M",1:ncol(samples))
    bf.names <- paste0("BF_", apply(combs,2,paste,collapse=""))
  }

  bf <- matrix(NA, ncol(combs), 5,
               dimnames = list(BF = bf.names,
                               Statistic = colnames(pp)))
  for(bb in 1:ncol(combs)){
    suppressWarnings(
      bf.tmp <- log(samples[,combs[1,bb]])-log(samples[,combs[2,bb]])
    )
    if(!logBF) bf.tmp <- exp(bf.tmp)
    bf[bb,] <- summ.samples(bf.tmp)
  }
  neff <- NA
  try(neff <- floor(dirichlet.mle(samples)$alpha0)) # sirt
  list(pp=pp, bf=bf, neff=neff)
}

summ.samples <- function(x)
  c(Mean=mean(x, na.rm = TRUE),
    SD=sd(x, na.rm = TRUE),
    quantile(x, probs=c(.05,.5,.95), na.rm = TRUE))
