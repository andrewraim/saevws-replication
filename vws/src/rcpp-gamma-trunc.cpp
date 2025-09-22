#include "rcpp-gamma-trunc.h"
#include "vws.h"

Rcpp::NumericVector r_gamma_trunc_rcpp(unsigned int n, double shape,
	double rate, double a, double b)
{
	return vws::r_gamma_trunc(n, shape, rate, a, b);
}

Rcpp::NumericVector d_gamma_trunc_rcpp(const Rcpp::NumericVector& x,
	double shape, double rate, double a, double b, bool log)
{
	return vws::d_gamma_trunc(x, shape, rate, a, b, log);
}

Rcpp::NumericVector p_gamma_trunc_rcpp(const Rcpp::NumericVector& q,
	double shape, double rate, double a, double b, bool lower, bool log)
{
	return vws::p_gamma_trunc(q, shape, rate, a, b, lower, log);
}

Rcpp::NumericVector q_gamma_trunc_rcpp(const Rcpp::NumericVector& p,
	double shape, double rate, double a, double b, bool lower, bool log)
{
	return vws::q_gamma_trunc(p, shape, rate,a, b, lower, log);
}
