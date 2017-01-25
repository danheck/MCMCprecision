#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <progress.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;


// sample posterior of (transposed) transition matrix P
// (conjugate Dirichlet(0,...,0) prior per row )
// [[Rcpp::export]]
arma::sp_mat rdirichletPt(arma::sp_mat Pt) {
  double colsum;
  for(arma::uword j=0; j<Pt.n_cols; j++) {
    for(arma::uword i=0; i<Pt.n_cols; i++) {
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
  for(arma::uword j=0; j<Pt.n_cols; j++) {
    for(arma::uword i=0; i<Pt.n_cols; i++) {
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
// sample: number of (independent) posterior samples
// [[Rcpp::export]]
arma::mat stationaryCpp(arma::mat freq,
                        int sample = 5000,
                        int cpu=4,
                        bool display_progress=true){
  int M = freq.n_cols;
  arma::mat samp(sample, M);
  arma::mat freqt = freq.t();
  int steps = round(1000/M);
  Progress p(round(sample/steps), display_progress);

#ifdef _OPENMP
  if ( cpu > 0 )
    omp_set_num_threads( cpu );
#endif

  bool run = true;
#pragma omp parallel for schedule(dynamic) private(run)
  for(int i=0; i<sample; i++){
    if(i%steps==0){
      p.increment();                  // update progress
      run = !Progress::check_abort(); // check whether to abort
    }
    if (run){

      // 1.) sample from conjugate posterior: Dirichlet
      arma::mat Pt = rdirichletPt(freqt);
      // 2.) get estimate for stationary distribution
      //     (normalized left eigenvector for eigenvalue = 1)
      try{
        arma::cx_vec eigval(M);
        arma::cx_mat eigvec(M,M);
        arma::eig_gen(eigval, eigvec, Pt);
        arma::uvec idx = find(round(eigval*1.e10)/1.e10 == 1.);
        arma::vec eigenvec1 = real(eigvec.cols(idx));
        // arma::vec eigenvec1 = real(eigvec.cols(0));
        samp.row(i) = (eigenvec1 / arma::accu(eigenvec1)).t();
      }catch(...) {
        samp.row(i).fill(arma::datum::nan);
        Rcout << "### RcppArmadillo::eig_gen unstable; method='base' might be more stable ###";
      }
    }
  }

  return samp;
}


// Posterior distribution for stationary distribution
// freq: matrix with transition frequencies
// a: prior for transition probabilities - Dirichlet(a,...,a)
// sample: number of (independent) posterior samples
// [[Rcpp::export]]
arma::mat stationaryCppSparse(arma::sp_mat freq,
                              int sample = 5000,
                              int cpu=4,
                              bool display_progress=true){
  int M = freq.n_cols;
  arma::mat samp(sample, M);
  arma::sp_mat freqt = freq.t();

  int steps = round(1000/M);
  Progress p(round(sample/steps), display_progress);

#ifdef _OPENMP
  if ( cpu > 0 )
    omp_set_num_threads( cpu );
#endif

  bool run = true;
#pragma omp parallel for schedule(dynamic) private(run)
for(int i=0; i<sample; i++){
  if(i%steps==0){
    p.increment();                  // update progress
    run = !Progress::check_abort(); // check whether to abort
  }
  if (run){

      // 1.) sample from conjugate posterior: Dirichlet
      arma::sp_mat Pt = rdirichletPt(freqt);
      // 2.) get estimate for stationary distribution
      //     (normalized left eigenvector for eigenvalue = 1)
      try{
        arma::cx_vec eigval(M);
        arma::cx_mat eigvec(M,1);
        arma::eigs_gen(eigval, eigvec, Pt, 1, "lr");
        arma::vec ev = real(eigvec);
        samp.row(i) = (ev / accu(ev)).t();
      }catch(...) {
        samp.row(i).fill(arma::datum::nan);
        Rcout << "### RcppArmadillo::eig_gen unstable; method='base'/'cpp' might be more stable ###";
      }
    }
  }
  return samp;
}


