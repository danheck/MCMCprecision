#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec sim_mc(int n,
                 arma::mat P,
                 int start)
{
  arma::vec chain(n);
  arma::vec p(P.n_cols);
  arma::vec idx = arma::linspace<arma::vec>(1,P.n_cols,P.n_cols);

  chain(0) = start;
  for(int i=1; i<n; i++)
  {
    p = P.row( chain(i-1)-1 ).t();
    chain(i) = arma::as_scalar(RcppArmadillo::sample(idx, 1, true, p) );
  }
  return chain;
}
