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
#' dirichlet.est(x)
#' @seealso \code{\link{rdirichlet}}
#' @references
#' Minka, T. (2000). Estimating a Dirichlet distribution. Technical Report.
#' @export
dirichlet.est <- function (x,
                           const = min(x)*.01,
                           maxit = 1000,
                           abstol = 1e-4)
{
  # adjust for x=0  (because of log(x) = -Inf)
  x <- (x + const)/(1 + 2 * const)
  x <- x/rowSums(x)

  # start values:
  x.mean <- colMeans(x)
  xlog.mean <- colMeans(log(x))
  x.squares <- colMeans(x^2)
  xi <- (x.mean - x.squares)/(x.squares - x.mean^2)
  alpha <- xi * x.mean

  x <- (x + const)/(1 + 2 * const)
  x <- x/rowSums(x)

  # start values:
  x.mean <- colMeans(x)
  x.squares <- colMeans(x^2)
  xi <- (x.mean - x.squares)/(x.squares - x.mean^2)
  alpha <- xi * x.mean

  logx.mean <- colMeans(log(x))
  N <- nrow(x)
  cnt <- diff <- 1
  min.x <- min(x)*.001
  alpha <- dirichlet_fp(alpha, logx.mean,
                        min = min.x, maxit = maxit, abstol = abstol)
  res <- list(alpha = alpha, sum = sum(alpha))
  return(res)
}
