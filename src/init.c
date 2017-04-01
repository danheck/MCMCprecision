#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP MCMCprecision_dirichlet_fp(SEXP, SEXP, SEXP, SEXP);
extern SEXP MCMCprecision_inv_digamma(SEXP, SEXP);
extern SEXP MCMCprecision_sim_mc(SEXP, SEXP, SEXP);
extern SEXP MCMCprecision_stationary_mle(SEXP, SEXP, SEXP, SEXP);
extern SEXP MCMCprecision_stationaryArma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MCMCprecision_stationaryArmaSparse(SEXP, SEXP, SEXP, SEXP);
extern SEXP MCMCprecision_stationaryEigen(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"MCMCprecision_dirichlet_fp",         (DL_FUNC) &MCMCprecision_dirichlet_fp,         4},
  {"MCMCprecision_inv_digamma",          (DL_FUNC) &MCMCprecision_inv_digamma,          2},
  {"MCMCprecision_sim_mc",               (DL_FUNC) &MCMCprecision_sim_mc,               3},
  {"MCMCprecision_stationary_mle",       (DL_FUNC) &MCMCprecision_stationary_mle,       4},
  {"MCMCprecision_stationaryArma",       (DL_FUNC) &MCMCprecision_stationaryArma,       5},
  {"MCMCprecision_stationaryArmaSparse", (DL_FUNC) &MCMCprecision_stationaryArmaSparse, 4},
  {"MCMCprecision_stationaryEigen",      (DL_FUNC) &MCMCprecision_stationaryEigen,      5},
  {NULL, NULL, 0}
};

void R_init_MCMCprecision(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
