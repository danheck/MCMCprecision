#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector stationary_mle (NumericVector pi,
                              NumericMatrix N,
                              double abstol = 1e-5,
                              int maxit = 1e5)
{
  NumericVector pi0(clone(pi));
  int M = N.cols();
  NumericVector N_row = rowSums(N);

  int cnt = 0;
  double diff = 1.;
  while ((diff > abstol) && (cnt < maxit))
  {
    pi0 = clone(pi);
    for (int i = 0; i < M; i++)
    {
      pi[i] = sum( (N.row(i) + N.column(i))/(N_row[i]/pi0[i] + N_row/pi0));
    }
    cnt += 1;
    diff = max(abs(pi0 - pi));
    // Rcout << cnt << " / ", diff << "::: " << pi << " // " << pi0 << "\n";
  }
  if (cnt == maxit)
    warning("Maximum number of iterations reached.");
  return (pi);
}
