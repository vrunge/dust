// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// madEstimator
double madEstimator(NumericVector& data);
RcppExport SEXP _dust_madEstimator(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(madEstimator(data));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_DUSTMODULE1D();

static const R_CallMethodDef CallEntries[] = {
    {"_dust_madEstimator", (DL_FUNC) &_dust_madEstimator, 1},
    {"_rcpp_module_boot_DUSTMODULE1D", (DL_FUNC) &_rcpp_module_boot_DUSTMODULE1D, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_dust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
