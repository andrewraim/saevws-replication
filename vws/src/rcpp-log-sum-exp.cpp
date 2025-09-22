#include "rcpp-log-sum-exp.h"
#include "vws.h"


double log_sum_exp_rcpp(const Rcpp::NumericVector& x)
{
	return vws::log_sum_exp(x);
}

Rcpp::NumericVector log_sum_exp_mat_rcpp(const Rcpp::NumericMatrix& x)
{
	return vws::log_sum_exp_mat(x);
}

Rcpp::NumericVector log_add2_exp_rcpp(const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& y)
{
	return vws::log_add2_exp(x, y);
}

Rcpp::NumericVector log_sub2_exp_rcpp(const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& y)
{
	return vws::log_sub2_exp(x, y);
}



