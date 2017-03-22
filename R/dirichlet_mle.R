#' Estimate Parameters of Dirichlet Distribution
#'
#' Fast C++ implementation of the fixed-point iteration algorithm by Minka (2000).
#'
#' @param x a matrix of Dirichlet samples, one row per observation
#' @param const constant that is added to zeros (to avoid problems in \code{log(x)})
#' @param maxit maximum number of iterations
#' @param abstol The absolute convergence tolerance: maximaum of absolute differences of Dirichlet parameters.
#' @examples
#' x <- rdirichlet(50, c(8,1,3,9))
#' dirichlet.mle(x)
#' @seealso \code{\link{rdirichlet}}
#' @references
#' Minka, T. (2000). Estimating a Dirichlet distribution. Technical Report.
#' @export
dirichlet.mle <- function (x,
                           const,
                           maxit = 1e5,
                           abstol = 1e-4)
{
  # adjust for x=0  (because of log(x) = -Inf)
  if (min(x) == 0){
    if (missing(const))
      const <- min(x[x>0])*.01
    x <- (x + const)/(1 + 2 * const)
  }
  x <- x/rowSums(x)
  logx.mean <- colMeans(log(x))
  N <- nrow(x)

  # heuristic for starting values:
  x.mean <- colMeans(x)
  x.squares <- colMeans(x^2)
  xi <- (x.mean - x.squares)/(x.squares - x.mean^2)
  alpha0 <- xi * x.mean
  alpha <- dirichlet_fp(pmax(.05, pmin(20, alpha0)), logx.mean,
                        maxit = maxit, abstol = abstol)
  # if this fails: random starting values
  if (anyNA(alpha))
    alpha <- dirichlet_fp(runif(length(alpha),0.5,1),
                          logx.mean,
                          maxit = maxit, abstol = abstol)

  res <- list(alpha = alpha, sum = sum(alpha))
  return(res)
}
