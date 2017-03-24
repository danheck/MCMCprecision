#' Generate a sample of a discrete-state Markov chain
#'
#' Generates a sequence of discrete states from a discrete-time Markov chain with transition matrix P.
#'
#' @param n length
#' @param P transition matrix (rows are normalized to sum to 1)
#' @param start starting distribution (discrete uniform by default)
#' @examples
#' P <- matrix(c(.3,.5,.2,
#'               .05,.25,.7,
#'               0,.1,.9), 3, 3, byrow=TRUE)
#' sim.mc(50, P)
#' @export
sim.mc <- function(n, P, start=rep(1, ncol(P))){

  if(!is.matrix(P) || ncol(P) != nrow(P) ||
     any(P<0) || ncol(P) < 2)
      stop("'P' must be a transition matrix with positive values.")
  P <- P/rowSums(P)
  if(!is.numeric(n) || n != round(n) || n <= 1 || length(n) != 1)
    stop("'n' must be a positive integer")
  if(!is.numeric(start) || !is.vector(start) || length(start) != ncol(P))
    stop("Length of 'start' does not match size of 'P'.")

  c0 <- sample(1:ncol(P), 1, prob=start)
  z <- sim_mc(n, P, c0)
  c(z)
}
