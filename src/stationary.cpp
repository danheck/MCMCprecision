#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// sample posterior of (transposed) transition matrix P
// (conjugate Dirichlet(0,...,0) prior per row )
// [[Rcpp::export]]
arma::sp_mat rdirichletPt(arma::sp_mat Pt) {
  double colsum;
  for(int j=0; j<Pt.n_cols; j++) {
    for(int i=0; i<Pt.n_cols; i++) {
      if(Pt(i,j) != 0)
        Pt(i,j) = R::rgamma(Pt(i,j), 1);
    }
    colsum = arma::accu(Pt.col(j));
    if(colsum > 0)
      Pt.col(j) /= colsum;
  }
  return( Pt );
}

arma::mat rdirichletPt(arma::mat Pt) {
  double colsum;
  for(int j=0; j<Pt.n_cols; j++) {
    for(int i=0; i<Pt.n_cols; i++) {
      if(Pt(i,j) != 0)
        Pt(i,j) = R::rgamma(Pt(i,j), 1);
    }
    colsum = arma::accu(Pt.col(j));
    if(colsum > 0)
      Pt.col(j) /= colsum;
  }
  return( Pt );
}



// Posterior distribution for stationary distribution
// freq: matrix with transition frequencies
// a: prior vector for transition probabilities - Dirichlet(a[1],...,a[M])
// iter: number of (independent) posterior samples
// [[Rcpp::export]]
arma::mat stationaryCpp(arma::mat freq,
                        int iter = 10000){
  int M = freq.n_cols;
  arma::mat samp(iter, M);
  arma::mat Pt,freqt = freq.t();
  arma::cx_vec eigval;
  arma::cx_mat eigvec;

  for(int i=0; i<iter; i++){
    // 1.) sample from conjugate posterior: Dirichlet
    Pt = rdirichletPt(freqt);
    // 2.) get estimate for stationary distribution
    //     (normalized left eigenvector for eigenvalue = 1)
    try{
      arma::eig_gen(eigval, eigvec, Pt);
      arma::uvec idx = find(round(eigval*1.e8)/1.e8 == 1.);
      arma::vec eigenvec1 = real(eigvec.cols(idx));
      samp.row(i) = (eigenvec1 / arma::accu(eigenvec1)).t();
    }catch(...) {
      samp.row(i).fill(arma::datum::nan);
      Rcout << "RcppArmadillo::eig_gen seems to be unstable. Use method='base'";
    }
  }

  return samp;
}


// Posterior distribution for stationary distribution
// freq: matrix with transition frequencies
// a: prior for transition probabilities - Dirichlet(a,...,a)
// iter: number of (independent) posterior samples
// [[Rcpp::export]]
arma::mat stationaryCppSparse(arma::sp_mat freq,
                              int iter = 10000){
  int M = freq.n_cols;
  arma::mat samp(iter, M);
  arma::sp_mat Pt,freqt = freq.t();
  arma::cx_vec eigval;
  arma::cx_mat eigvec;

  for(int i=0; i<iter; i++){
    // 1.) sample from conjugate posterior: Dirichlet
    Pt = rdirichletPt(freqt);
    // 2.) get estimate for stationary distribution
    //     (normalized left eigenvector for eigenvalue = 1)
    try{
      arma::eigs_gen(eigval, eigvec, Pt, 1, "lr");
      arma::vec ev = real(eigvec);
      samp.row(i) = (ev / accu(ev)).t();
    }catch(...) {
      samp.row(i).fill(arma::datum::nan);
      Rcout << "RcppArmadillo::eigs_gen seems to be unstable. Use method='base'";
    }
  }
  return samp;
}


