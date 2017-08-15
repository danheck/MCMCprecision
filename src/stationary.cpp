#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/rmultinom.h>
#include <RcppArmadilloExtensions/sample.h>
#include <progress.hpp>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace arma;

// [[Rcpp::export]]
arma::vec sim_mc(int n, arma::mat P, int start)
{
  vec chain(n);
  vec p(P.n_cols);
  vec idx = linspace<vec>(1,P.n_cols,P.n_cols);

  chain(0) = start;
  for(int i=1; i<n; i++)
  {
    p = P.row( chain(i-1)-1 ).t();
    chain(i) = as_scalar(sample(idx, 1, true, p) ); // Rcpp::RcppArmadillo::
  }
  return chain;
}



// sample posterior of (transposed) transition matrix P
// (conjugate Dirichlet(0,...,0) prior per row )
arma::sp_mat rdirichletPt(arma::sp_mat Pt)
{
  double colsum;
  for(uword j=0; j < Pt.n_cols; j++)
  {
    for(uword i=0; i < Pt.n_cols; i++)
    {
      if(Pt(i,j) != 0)
        Pt(i,j) = R::rgamma(Pt(i,j), 1);
    }
    colsum = accu(Pt.col(j));
    if(colsum > 0)
      Pt.col(j) /= colsum;
  }
  return (Pt);
}

// [[Rcpp::export]]
arma::mat rdirichletPt(arma::mat Pt)
{
  double colsum;
  for (uword j=0; j < Pt.n_cols; j++)
  {
    for (uword i=0; i < Pt.n_cols; i++)
    {
      if (Pt(i,j) != 0)
        Pt(i,j) = R::rgamma(Pt(i,j), 1);
    }
    colsum = accu(Pt.col(j));
    if(colsum > 0)
      Pt.col(j) /= colsum;
  }
  return (Pt);
}

// get expected transition probabilities of order 2: z[t-1, t, t+1]
// [[Rcpp::export]]
arma::cube getP2(arma::mat P, arma::vec pi)
{
  int M = P.n_cols;
  cube P2(M,M,M);
  for (int i = 0; i < M; i++)
    for (int j = 0; j < M; j++)
      for (int k = 0; k < M; k++)
        P2(i,j,k) = pi(i) * P(i,j) * P(j,k);
  return (P2);
}

double x2(arma::vec o, arma::vec e)
{
  return(accu(pow(o - e, 2) / e));
}

// [[Rcpp::export]]
arma::vec postpred(arma::mat P, arma::vec pi, arma::cube N2)
{
  vec statistic(2);
  statistic.fill(datum::nan);
  if (!N2.has_nan())
  {
    double N = accu(N2);
    vec P2vec = vectorise(getP2(P, pi)); // identical order as in R: c(P2)
    // observed:
    statistic(0) = x2(vectorise(N2), P2vec * N);
    // posterior predicted:
    ivec N2pred = rmultinom(N, as<NumericVector>(wrap(P2vec)));
    statistic(1) = x2(conv_to<vec>::from(N2pred), P2vec * N);
  }
  return(statistic);
}

arma::vec postpred(arma::sp_mat P, arma::vec pi, arma::cube N2)
{
  mat Pmat = conv_to<mat>::from(P);
  return (postpred(Pmat, pi, N2));
}

// Posterior distribution for stationary distribution
// N: matrix with transition frequencies
// a: prior vector for transition probabilities - Dirichlet(a[1],...,a[M])
// sample: number of (independent) posterior samples
// [[Rcpp::export]]
arma::mat stationaryArma(arma::mat N, arma::cube N2,
                         double epsilon = 0, int sample = 5000,
                         bool progress = true, double digits = 8.)
{
  int M = N.n_cols;
  int steps = round(1000/M);
  mat mcmc(M, sample), x2(2, sample);
  mcmc.fill(datum::nan);
  x2.fill(datum::nan);
  mat freqt = N.t() + epsilon;
  Progress p(sample, progress);

  uword maxIdx;
  cx_vec eigval;
  cx_mat eigvec;
  vec ev(M), pi(M);
  bool run = true;
  for(int i=0; i<sample; i++)
  {
    p.increment();   // update progress bar
    if(run && i % steps == 0)
      run = !Progress::check_abort(); // check whether to abort

    if (run)
    {
      // 1.) sample from conjugate posterior: Dirichlet
      mat Pt = rdirichletPt(freqt);
      try
      {
        // 2.) get estimate for stationary distribution
        //     (normalized left eigenvector for eigenvalue = 1)
        eig_gen(eigval, eigvec, Pt);
        maxIdx = index_max(real(eigval));
        // Rcout << "\n eigval: " << real(eigval.t()) << " maxIdx = " << maxIdx;
        if( abs(real(eigval(maxIdx)) - 1) < pow(10,-digits))
        {
          ev = real(eigvec.col(maxIdx));
          pi = ev / accu(ev);
          mcmc.col(i) = pi;
          // posterior predictive check:
          x2.col(i) = postpred(Pt.t(), pi, N2);
        }
      }
      catch(...)
      {
        Rcout << "# RcppArmadillo::eig_gen unstable: \n#" <<
          "method='base' or setting 'epsilon=.01' might provide more stable results#";
      }
    }
  }
  return join_vert(mcmc, x2).t();
}


// Posterior distribution for stationary distribution
// N: matrix with transition frequencies
// a: prior for transition probabilities - Dirichlet(a,...,a)
// sample: number of (independent) posterior samples
// [[Rcpp::export]]
arma::mat stationaryArmaSparse(arma::sp_mat N, arma::cube N2, double epsilon = 0,
                               int sample = 5000,
                               bool progress=true, double digits = 8.)
{
  int M = N.n_cols;
  mat mcmc(M, sample), x2(2, sample);
  mcmc.fill(datum::nan);
  sp_mat freqt = N.t();

  int steps = round(1000/M);
  Progress p(sample, progress);

  cx_vec eigval;
  cx_mat eigvec;
  vec ev, pi;
  bool run = true;
  for(int i=0; i<sample; i++)
  {
    p.increment();   // update progress bar
    if(run && i % steps == 0)
      run = !Progress::check_abort(); // check whether to abort

    if (run)
    {
      // 1.) sample from conjugate posterior: Dirichlet
      sp_mat Pt = rdirichletPt(freqt);
      // 2.) get estimate for stationary distribution
      //     (normalized left eigenvector for eigenvalue = 1)
      try
      {
        eigs_gen(eigval, eigvec, Pt, 1, "lr");
        if( abs(real(eigval(0)) - 1.) < pow(10,-digits))
        {
          ev = real(eigvec.col(0));
          pi = ev / accu(ev);
          mcmc.col(i) = pi;
        }
        // posterior predictive check
        x2.col(i) = postpred(Pt.t(), pi, N2);
      }
      catch(...)
      {
        Rcout << "# RcppArmadillo::eigs_gen unstable: \n#" <<
          "method='base' or 'epsilon=.01' might provide more stable results#";
      }
    }
  }
  return join_vert(mcmc, x2).t();
}


