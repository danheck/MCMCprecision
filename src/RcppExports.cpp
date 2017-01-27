// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sim_mc
arma::vec sim_mc(int n, arma::mat P, int start);
RcppExport SEXP MCMCprec_sim_mc(SEXP nSEXP, SEXP PSEXP, SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_mc(n, P, start));
    return rcpp_result_gen;
END_RCPP
}
// stationaryCpp
arma::mat stationaryCpp(arma::mat N, double epsilon, int sample, bool display_progress, int digits);
RcppExport SEXP MCMCprec_stationaryCpp(SEXP NSEXP, SEXP epsilonSEXP, SEXP sampleSEXP, SEXP display_progressSEXP, SEXP digitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< int >::type digits(digitsSEXP);
    rcpp_result_gen = Rcpp::wrap(stationaryCpp(N, epsilon, sample, display_progress, digits));
    return rcpp_result_gen;
END_RCPP
}
// stationaryCppSparse
arma::mat stationaryCppSparse(arma::sp_mat N, int sample, bool display_progress, int digits);
RcppExport SEXP MCMCprec_stationaryCppSparse(SEXP NSEXP, SEXP sampleSEXP, SEXP display_progressSEXP, SEXP digitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< int >::type digits(digitsSEXP);
    rcpp_result_gen = Rcpp::wrap(stationaryCppSparse(N, sample, display_progress, digits));
    return rcpp_result_gen;
END_RCPP
}
