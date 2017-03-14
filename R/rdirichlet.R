#' Random Sample from Dirichlet Distribution
#'
#' Random generation from the Dirichlet distribution.
#'
#' @param n number of samples
#' @param a vector or matrix of shape parameters
#' @examples
#' rdirichlet(2, c(1,5,3,8))
#'
#' @seealso \code{\link{dirichlet.est}}
#' @export
rdirichlet <- function (n, a){

  # parameters in columns:
  a <- rbind(a)
  M <- ncol(a)

  # for a single set of parameters, extend to matrix
  if (n > nrow(a))
    a <- matrix(a, n, M, byrow = TRUE)

  x <- matrix(rgamma(M*n, a), nrow=n, ncol=M)

  # normalize to sum up to one:
  x/rowSums(x)
}
