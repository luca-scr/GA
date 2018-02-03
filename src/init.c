#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* 
Routines registration obtained with 

tools::package_native_routine_registration_skeleton(".", character_only = FALSE)
 
FIXME: Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _GA_c_double(SEXP, SEXP);
extern SEXP _GA_c_int(SEXP, SEXP);
extern SEXP _GA_ga_lrSelection_Rcpp(SEXP, SEXP, SEXP);
extern SEXP _GA_ga_nlrSelection_Rcpp(SEXP, SEXP);
extern SEXP _GA_ga_pmutation_Rcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GA_ga_rwSelection_Rcpp(SEXP);
extern SEXP _GA_ga_spCrossover_Rcpp(SEXP, SEXP);
extern SEXP _GA_ga_tourSelection_Rcpp(SEXP, SEXP);
extern SEXP _GA_gabin_Population_Rcpp(SEXP);
extern SEXP _GA_gabin_raMutation_Rcpp(SEXP, SEXP);
extern SEXP _GA_gabin_uCrossover_Rcpp(SEXP, SEXP);
extern SEXP _GA_gaperm_cxCrossover_Rcpp(SEXP, SEXP);
extern SEXP _GA_gaperm_dmMutation_Rcpp(SEXP, SEXP);
extern SEXP _GA_gaperm_ismMutation_Rcpp(SEXP, SEXP);
extern SEXP _GA_gaperm_oxCrossover_Rcpp(SEXP, SEXP);
extern SEXP _GA_gaperm_pbxCrossover_Rcpp(SEXP, SEXP);
extern SEXP _GA_gaperm_pmxCrossover_Rcpp(SEXP, SEXP);
extern SEXP _GA_gaperm_Population_Rcpp(SEXP);
extern SEXP _GA_gaperm_scrMutation_Rcpp(SEXP, SEXP);
extern SEXP _GA_gaperm_simMutation_Rcpp(SEXP, SEXP);
extern SEXP _GA_gaperm_swMutation_Rcpp(SEXP, SEXP);
extern SEXP _GA_gareal_blxCrossover_Rcpp(SEXP, SEXP);
extern SEXP _GA_gareal_laCrossover_Rcpp(SEXP, SEXP);
extern SEXP _GA_gareal_laplaceCrossover_Rcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GA_gareal_lsSelection_Rcpp(SEXP);
extern SEXP _GA_gareal_nraMutation_Rcpp(SEXP, SEXP);
extern SEXP _GA_gareal_Population_Rcpp(SEXP);
extern SEXP _GA_gareal_powMutation_Rcpp(SEXP, SEXP, SEXP);
extern SEXP _GA_gareal_raMutation_Rcpp(SEXP, SEXP);
extern SEXP _GA_gareal_rsMutation_Rcpp(SEXP, SEXP);
extern SEXP _GA_gareal_sigmaSelection_Rcpp(SEXP);
extern SEXP _GA_gareal_waCrossover_Rcpp(SEXP, SEXP);
extern SEXP _GA_intersect_asR(SEXP, SEXP);
extern SEXP _GA_optimProbsel_Rcpp(SEXP, SEXP);
extern SEXP _GA_rank_asR(SEXP, SEXP);
extern SEXP _GA_round_double(SEXP, SEXP);
extern SEXP _GA_setdiff_asR(SEXP, SEXP);
extern SEXP _GA_which_asR(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_GA_c_double",                     (DL_FUNC) &_GA_c_double,                     2},
    {"_GA_c_int",                        (DL_FUNC) &_GA_c_int,                        2},
    {"_GA_ga_lrSelection_Rcpp",          (DL_FUNC) &_GA_ga_lrSelection_Rcpp,          3},
    {"_GA_ga_nlrSelection_Rcpp",         (DL_FUNC) &_GA_ga_nlrSelection_Rcpp,         2},
    {"_GA_ga_pmutation_Rcpp",            (DL_FUNC) &_GA_ga_pmutation_Rcpp,            4},
    {"_GA_ga_rwSelection_Rcpp",          (DL_FUNC) &_GA_ga_rwSelection_Rcpp,          1},
    {"_GA_ga_spCrossover_Rcpp",          (DL_FUNC) &_GA_ga_spCrossover_Rcpp,          2},
    {"_GA_ga_tourSelection_Rcpp",        (DL_FUNC) &_GA_ga_tourSelection_Rcpp,        2},
    {"_GA_gabin_Population_Rcpp",        (DL_FUNC) &_GA_gabin_Population_Rcpp,        1},
    {"_GA_gabin_raMutation_Rcpp",        (DL_FUNC) &_GA_gabin_raMutation_Rcpp,        2},
    {"_GA_gabin_uCrossover_Rcpp",        (DL_FUNC) &_GA_gabin_uCrossover_Rcpp,        2},
    {"_GA_gaperm_cxCrossover_Rcpp",      (DL_FUNC) &_GA_gaperm_cxCrossover_Rcpp,      2},
    {"_GA_gaperm_dmMutation_Rcpp",       (DL_FUNC) &_GA_gaperm_dmMutation_Rcpp,       2},
    {"_GA_gaperm_ismMutation_Rcpp",      (DL_FUNC) &_GA_gaperm_ismMutation_Rcpp,      2},
    {"_GA_gaperm_oxCrossover_Rcpp",      (DL_FUNC) &_GA_gaperm_oxCrossover_Rcpp,      2},
    {"_GA_gaperm_pbxCrossover_Rcpp",     (DL_FUNC) &_GA_gaperm_pbxCrossover_Rcpp,     2},
    {"_GA_gaperm_pmxCrossover_Rcpp",     (DL_FUNC) &_GA_gaperm_pmxCrossover_Rcpp,     2},
    {"_GA_gaperm_Population_Rcpp",       (DL_FUNC) &_GA_gaperm_Population_Rcpp,       1},
    {"_GA_gaperm_scrMutation_Rcpp",      (DL_FUNC) &_GA_gaperm_scrMutation_Rcpp,      2},
    {"_GA_gaperm_simMutation_Rcpp",      (DL_FUNC) &_GA_gaperm_simMutation_Rcpp,      2},
    {"_GA_gaperm_swMutation_Rcpp",       (DL_FUNC) &_GA_gaperm_swMutation_Rcpp,       2},
    {"_GA_gareal_blxCrossover_Rcpp",     (DL_FUNC) &_GA_gareal_blxCrossover_Rcpp,     2},
    {"_GA_gareal_laCrossover_Rcpp",      (DL_FUNC) &_GA_gareal_laCrossover_Rcpp,      2},
    {"_GA_gareal_laplaceCrossover_Rcpp", (DL_FUNC) &_GA_gareal_laplaceCrossover_Rcpp, 4},
    {"_GA_gareal_lsSelection_Rcpp",      (DL_FUNC) &_GA_gareal_lsSelection_Rcpp,      1},
    {"_GA_gareal_nraMutation_Rcpp",      (DL_FUNC) &_GA_gareal_nraMutation_Rcpp,      2},
    {"_GA_gareal_Population_Rcpp",       (DL_FUNC) &_GA_gareal_Population_Rcpp,       1},
    {"_GA_gareal_powMutation_Rcpp",      (DL_FUNC) &_GA_gareal_powMutation_Rcpp,      3},
    {"_GA_gareal_raMutation_Rcpp",       (DL_FUNC) &_GA_gareal_raMutation_Rcpp,       2},
    {"_GA_gareal_rsMutation_Rcpp",       (DL_FUNC) &_GA_gareal_rsMutation_Rcpp,       2},
    {"_GA_gareal_sigmaSelection_Rcpp",   (DL_FUNC) &_GA_gareal_sigmaSelection_Rcpp,   1},
    {"_GA_gareal_waCrossover_Rcpp",      (DL_FUNC) &_GA_gareal_waCrossover_Rcpp,      2},
    {"_GA_intersect_asR",                (DL_FUNC) &_GA_intersect_asR,                2},
    {"_GA_optimProbsel_Rcpp",            (DL_FUNC) &_GA_optimProbsel_Rcpp,            2},
    {"_GA_rank_asR",                     (DL_FUNC) &_GA_rank_asR,                     2},
    {"_GA_round_double",                 (DL_FUNC) &_GA_round_double,                 2},
    {"_GA_setdiff_asR",                  (DL_FUNC) &_GA_setdiff_asR,                  2},
    {"_GA_which_asR",                    (DL_FUNC) &_GA_which_asR,                    1},
    {NULL, NULL, 0}
};

void R_init_GA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}