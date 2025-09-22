#ifndef RCPP_UNIF_H
#define RCPP_UNIF_H

#include <Rcpp.h>

//' Uniform Distribution
//'
//' Moment-generating function (MGF) for the \eqn{Uniform(a, b)} distribution.
//'
//' @param s Vector; argument of MGF.
//' @param lo Lower bound \eqn{a} of distribution.
//' @param hi Lower bound \eqn{b} of distribution.
//' @param log logical; if `TRUE`, return value on log-scale. Otherwise return
//' value on original scale.
//'
//' @return Value of MGF.
//'
//' @examples
//' mgf_unif(0, lo = 0, hi = 1, log = TRUE)
//' mgf_unif(-0.25, lo = 0, hi = 1, log = TRUE)
//' mgf_unif(-0.25, lo = 0, hi = 1, log = FALSE)
//' mgf_unif(+0.25, lo = 0, hi = 1, log = TRUE)
//'
//' @name Uniform
//' @export
// [[Rcpp::export(name = "mgf_unif")]]
Rcpp::NumericVector mgf_unif_rcpp(const Rcpp::NumericVector& s, double lo,
	double hi, bool log = false);

#endif
