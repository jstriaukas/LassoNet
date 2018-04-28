// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// betanew_lasso_cpp
List betanew_lasso_cpp(NumericMatrix xx, NumericVector xy, NumericVector beta, NumericMatrix M, NumericVector y, double Lambda1, double Lambda2, double iter, double tol);
RcppExport SEXP LassoNet_betanew_lasso_cpp(SEXP xxSEXP, SEXP xySEXP, SEXP betaSEXP, SEXP MSEXP, SEXP ySEXP, SEXP Lambda1SEXP, SEXP Lambda2SEXP, SEXP iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xy(xySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type Lambda1(Lambda1SEXP);
    Rcpp::traits::input_parameter< double >::type Lambda2(Lambda2SEXP);
    Rcpp::traits::input_parameter< double >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    __result = Rcpp::wrap(betanew_lasso_cpp(xx, xy, beta, M, y, Lambda1, Lambda2, iter, tol));
    return __result;
END_RCPP
}
