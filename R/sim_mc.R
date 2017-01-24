
#' Generate a sample of a discrete-state Markov chain
#'
#' Generates a sequence of discrete states from a discrete-time Markov chain with transition matrix P.
#'
#' @param n length
#' @param P transition matrix
#' @param start starting distribution (discrete uniform by default)
#' @examples
#' P <- matrix(c(.3,.5,.2,
#'               .05,.25,.7,
#'               0,.1,.9), 3, 3, byrow=TRUE)
#' sim.mc(50, P)
#' @export
sim.mc <- function(n, P, start=rep(1, ncol(P))){

  if(!is.matrix(P) || ncol(P) != nrow(P) ||
     any(P<0) || any(rowSums(P) != 1) || ncol(P) <= 2)
      stop("'P' must be a transition matrix with positive values and rows summing to one.")
  if(!is.numeric(n) || n != round(n) || n <= 1 || length(n) != 1)
    stop("'n' must be a positive integer")
  if(!is.numeric(start) || !is.vector(start) || length(start) != ncol(P))
    stop("Length of 'start' does not match size of 'P'.")

  c0 <- sample(1:ncol(P), 1, prob=start)
  z <- sim_mc(n, P, c0)
  c(z)
}


# Generate transition matrix for a given stationary distribution
#
# Returns a transition matrix with a specific autocorrelation. Useful for simulations.
#
# @param pi: stationary distribution
# @param alpha: probability to stay with current model (=> autocorrelation)
# @expamples
# pi <- c(.1,.5,.01,.3, .09)
# P <- make.P(pi, alpha=.5)
# P
# @export
# make.P <- function(pi, alpha = 0){
#   pi <- pi/sum(pi)
#   M <- length(pi)
#   P <- alpha*diag(M)  + (1-alpha)*matrix(pi,M,M,byrow = TRUE)
#   P
# }



# Simulation of Markov chain estimation
#
# @param n: length of chain
# @param pi: true stationary distribution
# @param alpha: probability to switch to other model (distributed uniformly)
# @param chains: number of Markov chains
# @param R: number of simulation replications
# @param iter: number of (independent) Dirichlet posterior samples
# @examples
#
# @export
# sim.chain <- function(n, P, chains=1,
#                       R=100, iter=5000){
#
#   if(!is.matrix(P) || ncol(P) != nrow(P) || any(P<0) || any(rowSums(P) != 1))
#     stop("'P' must be a transition matrix with positive values and rows summing to one.")
#   if(!is.numeric(n) || n != round(n) || n < 1 || length(n) != 1)
#     stop("'n' must be a positive integer")
#
#   M <- ncol(P)
#   est <- array(NA, c(R, 2, ncol(P), 5),
#                dimnames=list(Replication=NULL, Method=c("Markov","iid"),
#                              Model = paste0("Model ",1:ncol(P)),
#                              Statistic=c("Mean","SD","5%","50%","95%")))
#   for(i in 1:R){
#     obs <- gen.chain(n, P, chains)
#     s1 <- stationary(obs, labels = 1:M, iter = iter, method = "cpp")
#     s2 <- stationary(obs, labels = 1:M, iter = iter, method = "iid")
#     est[i,1,,] <- s1$pp
#     est[i,2,,] <- s2$pp
#   }
#   est
# }
