#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _MCMCprecision_dirichlet_fp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _MCMCprecision_inv_digamma(SEXP, SEXP);
extern SEXP _MCMCprecision_sim_mc(SEXP, SEXP, SEXP);
extern SEXP _MCMCprecision_stationary_mle(SEXP, SEXP, SEXP, SEXP);
extern SEXP _MCMCprecision_stationaryArma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MCMCprecision_stationaryArmaSparse(SEXP, SEXP, SEXP, SEXP);
extern SEXP _MCMCprecision_stationaryEigen(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_MCMCprecision_dirichlet_fp",         (DL_FUNC) &_MCMCprecision_dirichlet_fp,         4},
    {"_MCMCprecision_inv_digamma",          (DL_FUNC) &_MCMCprecision_inv_digamma,          2},
    {"_MCMCprecision_sim_mc",               (DL_FUNC) &_MCMCprecision_sim_mc,               3},
    {"_MCMCprecision_stationary_mle",       (DL_FUNC) &_MCMCprecision_stationary_mle,       4},
    {"_MCMCprecision_stationaryArma",       (DL_FUNC) &_MCMCprecision_stationaryArma,       5},
    {"_MCMCprecision_stationaryArmaSparse", (DL_FUNC) &_MCMCprecision_stationaryArmaSparse, 4},
    {"_MCMCprecision_stationaryEigen",      (DL_FUNC) &_MCMCprecision_stationaryEigen,      5},
    {NULL, NULL, 0}
};

void R_init_MCMCprecision(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}