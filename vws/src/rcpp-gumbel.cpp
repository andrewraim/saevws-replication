#include "rcpp-gumbel.h"
#include "vws.h"

Rcpp::NumericVector r_gumbel_rcpp(unsigned int n, double mu, double sigma)
{
	return vws::r_gumbel(n, mu, sigma);
}

Rcpp::NumericVector d_gumbel_rcpp(const Rcpp::NumericVector& x, double mu,
	double sigma, bool log)
{
	return vws::d_gumbel(x, mu, sigma, log);
}

Rcpp::NumericVector p_gumbel_rcpp(const Rcpp::NumericVector& q, double mu,
	double sigma, bool lower, bool log)
{
	return vws::p_gumbel(q, mu, sigma, lower, log);
}

Rcpp::NumericVector q_gumbel_rcpp(const Rcpp::NumericVector& p, double mu,
	double sigma, bool lower, bool log)
{
	return vws::q_gumbel(p, mu, sigma, lower, log);
}
