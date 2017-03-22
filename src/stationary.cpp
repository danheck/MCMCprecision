#include <RcppArmadillo.h>
#include <progress.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;


// sample posterior of (transposed) transition matrix P
// (conjugate Dirichlet(0,...,0) prior per row )
arma::sp_mat rdirichletPt(arma::sp_mat Pt)
{
  double colsum;
  for(arma::uword j=0; j<Pt.n_cols; j++)
  {
    for(arma::uword i=0; i<Pt.n_cols; i++)
    {
      if(Pt(i,j) != 0)
        Pt(i,j) = R::rgamma(Pt(i,j), 1);
    }
    colsum = arma::accu(Pt.col(j));
    if(colsum > 0)
      Pt.col(j) /= colsum;
  }
  return (Pt);
}

arma::mat rdirichletPt(arma::mat Pt)
{
  double colsum;
  for (arma::uword j=0; j<Pt.n_cols; j++)
  {
    for (arma::uword i=0; i<Pt.n_cols; i++)
    {
      if (Pt(i,j) != 0)
        Pt(i,j) = R::rgamma(Pt(i,j), 1);
    }
    colsum = arma::accu(Pt.col(j));
    if(colsum > 0)
      Pt.col(j) /= colsum;
  }
  return (Pt);
}



// Posterior distribution for stationary distribution
// N: matrix with transition frequencies
// a: prior vector for transition probabilities - Dirichlet(a[1],...,a[M])
// sample: number of (independent) posterior samples
// [[Rcpp::export]]
arma::mat stationaryArma(arma::mat N,
                         double epsilon = 0,
                         int sample = 5000,
                         bool display_progress=true,
                         double digits = 8.)
{
  int M = N.n_cols;
  int steps = round(1000/M);
  arma::mat samp(sample, M);
  samp.fill(arma::datum::nan);
  arma::mat freqt = N.t() + epsilon;
  Progress p(sample, display_progress);

  arma::uword maxIdx;
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  bool run = true;
  for(int i=0; i<sample; i++)
  {
    p.increment();   // update progress bar
    if(run && i % steps == 0)
      run = !Progress::check_abort(); // check whether to abort

    if (run)
    {
      // 1.) sample from conjugate posterior: Dirichlet
      arma::mat Pt = rdirichletPt(freqt);
      // 2.) get estimate for stationary distribution
      //     (normalized left eigenvector for eigenvalue = 1)
      try
      {
        arma::eig_gen(eigval, eigvec, Pt);
        maxIdx = arma::index_max(real(eigval));
        // Rcout << "\n eigval: " << real(eigval.t()) << " maxIdx = " << maxIdx;
        if( abs(real(eigval(maxIdx)) - 1.) < pow(10,-digits))
        {
          arma::vec ev = real(eigvec.col(maxIdx));
          samp.row(i) = (ev / accu(ev)).t();
        }
      }
      catch(...)
      {
        Rcout << "# RcppArmadillo::eigs_gen unstable: \n#" <<
          "method='base' or 'epsilon=.01' might provide more stable results#";
      }
    }
  }
  return samp;
}


// Posterior distribution for stationary distribution
// N: matrix with transition frequencies
// a: prior for transition probabilities - Dirichlet(a,...,a)
// sample: number of (independent) posterior samples
// [[Rcpp::export]]
arma::mat stationaryArmaSparse(arma::sp_mat N,
                               int sample = 5000,
                               bool display_progress=true,
                               double digits = 8.)
{
  int M = N.n_cols;
  arma::mat samp(sample, M);
  samp.fill(arma::datum::nan);
  arma::sp_mat freqt = N.t();

  int steps = round(1000/M);
  Progress p(sample, display_progress);

  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  bool run = true;
  for(int i=0; i<sample; i++)
  {
    p.increment();   // update progress bar
    if(run && i % steps == 0)
      run = !Progress::check_abort(); // check whether to abort

    if (run)
    {
      // 1.) sample from conjugate posterior: Dirichlet
      arma::sp_mat Pt = rdirichletPt(freqt);
      // 2.) get estimate for stationary distribution
      //     (normalized left eigenvector for eigenvalue = 1)
      try
      {
        arma::eigs_gen(eigval, eigvec, Pt, 1, "lr");
        if( abs(real(eigval(0)) - 1.) < pow(10,-digits))
        {
          arma::vec ev = real(eigvec.col(0));
          samp.row(i) = (ev / accu(ev)).t();
        }
      }
      catch(...)
      {
        Rcout << "# RcppArmadillo::eigs_gen unstable: \n#" <<
          "method='base' or 'epsilon=.01' might provide more stable results#";
      }
    }
  }
  return samp;
}


