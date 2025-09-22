#include "rcpp-invgamma.h"
#include "vws.h"

Rcpp::NumericVector r_invgamma_rcpp(unsigned int n, double a, double b)
{
	return vws::r_invgamma(n, a, b);
}

Rcpp::NumericVector d_invgamma_rcpp(const Rcpp::NumericVector& x, double a,
	double b, bool log)
{
	return vws::d_invgamma(x, a, b, log);
}

Rcpp::NumericVector p_invgamma_rcpp(const Rcpp::NumericVector& q, double a,
	double b, bool lower, bool log)
{
	return vws::p_invgamma(q, a, b, lower, log);
}

Rcpp::NumericVector q_invgamma_rcpp(const Rcpp::NumericVector& p, double a,
	double b, bool lower, bool log)
{
	return vws::q_invgamma(p, a, b, lower, log);
}
