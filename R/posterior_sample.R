# #' Get a sample from the posterior for the stationary distribution
# #'
# #' @param i ignored (just to use parSapply)
# #' @param tab matrix of transition frequencies
# #' @param M number of models
# #' @export
posterior.sample <- function (i,
                              tab,
                              epsilon = 0,
                              # method = "base",
                              digits = 8){
  M <- ncol(tab)

  # 1.) sample from conjugate posterior with prior: Dirichlet(0,..,0)
  P <- matrix(rgamma(M^2, tab + epsilon, 1),
                    nrow = M, ncol = M)
  sel <- rowSums(P) > 0
  P[sel,] <- P[sel,,drop=FALSE]/rowSums(P[sel,,drop=FALSE])
  # 2.) get estimate for stationary distribution (largest eigenvalue = 1)
  # if (method == "base"){
  decomp <- eigen(t(P))
  # # }else if (method == "sparse"){
  # #   decomp <- eigs(Matrix(t(P), sparse=TRUE), k = 1, which = "LR") ## rARPACK
  # } else {
  #   stop ("Method not supported.")
  # }
  idx <- which.max(Re(decomp$values))
  if (round(decomp$values[idx],digits) != 1){
    return (rep(NA, M))
  } else {
    ev <- Re(decomp$vectors[,idx])
    return (ev/sum(ev))
  }
}
