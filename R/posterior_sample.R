# #' Get a sample from the posterior for the stationary distribution
# #'
# #' @param i ignored (just to use parSapply below)
# #' @param tab matrix of transition frequencies
# #' @param M number of models
# #' @export
# @importFrom LaplacesDemon rdirichlet
#' @importFrom stats rgamma
posterior.sample <- function(i, tab, method = "sparse"){
  M <- ncol(tab)

  # 1.) sample from conjugate posterior with prior: Dirichlet(0,..,0)
  P <- aa <- matrix(rgamma(M^2, tab, 1),
                    nrow = M, ncol = M)
  sel <- rowSums(aa) > 0
  P[sel,] <- aa[sel,,drop=FALSE]/rowSums(aa[sel,,drop=FALSE])
  # 2.) get estimate for stationary distribution (largest eigenvalue = 1)
  if(method == "base"){
    decomp <- eigen(t(P))
  # }else if (method == "sparse"){
  #   decomp <- eigs(Matrix(t(P), sparse=TRUE), k = 1, which = "LR") ## rARPACK
  }else{
    stop("Method not supported.")
  }
  if(round(decomp$values[1],8) != 1){
    return( rep(NA, M) )
  }else{
    ev <- Re(decomp$vectors[,1])
    return( ev/sum(ev) )
  }
}
