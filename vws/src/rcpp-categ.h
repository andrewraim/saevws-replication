#ifndef RCPP_CATEG_H
#define RCPP_CATEG_H

#include <Rcpp.h>

//' Categorical Distribution
//'
//' Draw variates from a categorical distribution.
//'
//' @param n Number of desired draws.
//' @param p Vector of \eqn{k} probabilities for distribution.
//' @param log logical; if `TRUE`, interpret `p` as being specified on the
//' log-scale as `log(p)`. Otherwise, interpret `p` as being specified on the
//' original probability scale.
//' @param one_based logical; if `TRUE`, assume a categorical distribution
//' with support \eqn{\{ 1, \ldots, k \}}. Otherwise, assume support
//' \eqn{\{ 0, \ldots, k - 1 \}}.
//'
//' @return A vector of \eqn{n} draws.
//'
//' @examples
//' p = c(0.1, 0.2, 0.3, 0.4)
//' lp = log(p)
//'
//' set.seed(1234)
//' r_categ(50, p, log = FALSE, one_based = FALSE)
//' r_categ(50, p, log = FALSE, one_based = TRUE)
//'
//' set.seed(1234)
//' r_categ(50, lp, log = TRUE, one_based = FALSE)
//' r_categ(50, lp, log = TRUE, one_based = TRUE)
//'
//' @name Categorical
//' @export
// [[Rcpp::export(name = "r_categ")]]
Rcpp::IntegerVector r_categ_rcpp(unsigned int n, const Rcpp::NumericVector& p,
	bool log = false, bool one_based = false);

#endif

