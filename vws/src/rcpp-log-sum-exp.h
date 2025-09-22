#ifndef RCPP_LOG_SUM_EXP_H
#define RCPP_LOG_SUM_EXP_H

#include <Rcpp.h>

//' Log-Sum-Exp
//'
//' Compute `log(sum(exp(x)))` in a more stable way.
//'
//' @param x A numeric vector
//' @param y A numeric vector
//'
//' @details Computed using the method described by user Ben in StackExchange
//' thread \url{https://stats.stackexchange.com/questions/381936/vectorised-computation-of-logsumexp}.
//'
//' The function `log_sub2_exp` expects that each element of `x` is
//' larger than or equal to its corresponding element in `y`. Otherwise,
//' `NaN` will be returned with a warning.
//'
//' @examples
//' pi = 1:6 / sum(1:6)
//' x = log(2*pi)
//' log(sum(exp(x)))
//' log_sum_exp(x)
//'
//' # Result should be 0
//' x = c(-Inf -Inf, 0)
//' log_sum_exp(x)
//'
//' # Result should be -Inf
//' x = c(-Inf -Inf, -Inf)
//' log_sum_exp(x)
//'
//' # Result should be Inf
//' x = c(-Inf -Inf, Inf)
//' log_sum_exp(x)
//'
//' # Result should be 5 on the original scale
//' out = log_add2_exp(log(3), log(2))
//' exp(out)
//'
//' # Result should be 7 on the original scale
//' out = log_sub2_exp(log(12), log(5))
//' exp(out)
//'
//' @name log_sum_exp
//' @export
// [[Rcpp::export(name = "log_sum_exp")]]
double log_sum_exp_rcpp(const Rcpp::NumericVector& x);

//' @name log_sum_exp
//' @export
// [[Rcpp::export(name = "log_sum_exp_mat")]]
Rcpp::NumericVector log_sum_exp_mat_rcpp(const Rcpp::NumericMatrix& x);

//' @name log_sum_exp
//' @export
// [[Rcpp::export(name = "log_add2_exp")]]
Rcpp::NumericVector log_add2_exp_rcpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& y);

//' @name log_sum_exp
//' @export
// [[Rcpp::export(name = "log_sub2_exp")]]
Rcpp::NumericVector log_sub2_exp_rcpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& y);

#endif

