#include "rcpp-texp.h"
#include "vws.h"

double n_texp_rcpp(double kappa, double lo, double hi, bool log)
{
	return vws::n_texp(kappa, lo, hi, log);
}

double integrate_texp_rcpp(double a, double b, double kappa, double lo,
	double hi, bool log)
{
	return vws::integrate_texp(a, b, kappa, lo, hi, log);
}

Rcpp::NumericVector d_texp_rcpp(const Rcpp::NumericVector& x, double kappa,
	double lo, double hi, bool log)
{
	return vws::d_texp(x, kappa, lo, hi, log);
}

Rcpp::NumericVector p_texp_rcpp(const Rcpp::NumericVector& q, double kappa,
	double lo, double hi, bool lower, bool log)
{
	return vws::p_texp(q, kappa, lo, hi, lower, log);
}

Rcpp::NumericVector q_texp_rcpp(const Rcpp::NumericVector& p, double kappa,
	double lo, double hi, bool lower, bool log)
{
	return vws::q_texp(p, kappa, lo, hi, lower, log);
}

Rcpp::NumericVector r_texp_rcpp(unsigned int n, double kappa, double lo, double hi)
{
	return vws::r_texp(n, kappa, lo, hi);
}

Rcpp::NumericVector mgf_texp_rcpp(const Rcpp::NumericVector& s, double kappa,
	double lo, double hi, bool log)
{
	return vws::mgf_texp(s, kappa, lo, hi, log);
}
