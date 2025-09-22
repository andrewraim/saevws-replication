#ifndef VWS_GAMMA_TRUNC_H
#define VWS_GAMMA_TRUNC_H

#include <Rcpp.h>
#include "log-sum-exp.h"

/*
* Functions for the Gamma distribution with shape and rate parameters,
* truncated to the interval $[a, b]$.
*/

namespace vws {

/*
* Density of Truncated Gamma distribution
*
* - `x`: density argument.
* - `shape`: shape parameter.
* - `rate`: rate parameter.
* - `a`: lower limit of support.
* - `b`: upper limit of support.
* - `log`: if `true`, return log of value. Otherwise return value on the
*   original scale.
*
* Returns density value.
*/
inline double d_gamma_trunc(double x, double shape, double rate, double a, double b,
	bool log)
{
	double scale = 1 / rate;
	double log_num = a < x && x <= b ? R::dgamma(x, shape, scale, true) : R_NegInf;

	// Compute two ways: using upper tail and lower tail CDF. One of these may
	// be -Inf if x is far into the tail of N(mu, sigma2); in this case, the
	// other one may be finite.

	double lpa = R::pgamma(a, shape, scale, true, true);
	double lpb = R::pgamma(b, shape, scale, true, true);
	double clpa = R::pgamma(a, shape, scale, false, true);
	double clpb = R::pgamma(b, shape, scale, false, true);
	double log_den1 = log_sub2_exp(lpb, lpa);
	double log_den2 = log_sub2_exp(clpa, clpb);
	double log_den = std::max(log_den1, log_den2);

	double out = log_num - log_den;
	return log ? out : exp(out);
}

/*
* CDF of Truncated Gamma distribution
*
* - `x`: CDF argument.
* - `shape`: shape parameter.
* - `rate`: rate parameter.
* - `a`: lower limit of support.
* - `b`: upper limit of support.
* - `log`: if `true`, return log of value. Otherwise return value on the
*   original scale.
*
* Returns CDF value.
*/
inline double p_gamma_trunc(double x, double shape, double rate,
	double a, double b, bool lower, bool log)
{
	double scale = 1 / rate;

	if (x <= a || x > b) {
		double out = R_NegInf;
		return log ? out : exp(out);
	}

	double lp_a = R::pgamma(a, shape, scale, true, true);
	double lp_b = R::pgamma(b, shape, scale, true, true);
	double lp_x = R::pgamma(x, shape, scale, true, true);
	double lp_den = log_sub2_exp(lp_b, lp_a);

	double clp_a = R::pgamma(a, shape, scale, false, true);
	double clp_b = R::pgamma(b, shape, scale, false, true);
	double clp_x = R::pgamma(x, shape, scale, false, true);
	double clp_den = log_sub2_exp(clp_a, clp_b);

	double out = (lp_den < clp_den) ?
		log_sub2_exp(lp_x, lp_a) - lp_den :
		log_sub2_exp(clp_a, clp_x) - clp_den;

	return log ? out : exp(out);
}

/*
* Quantile function of Truncated Gamma distribution
*
* - `p`: requested quantile.
* - `shape`: shape parameter.
* - `rate`: rate parameter.
* - `a`: lower limit of support.
* - `b`: upper limit of support.
* - `lower`: if `true`, request $p$ quantile. Otherwise request $1-p$ quantile.
* - `log`: if `true`, assume $p$ is specified on the log-scale. Otherwise,
*    assume it is on the original scale.
*
* Returns quantile value.
*/
inline double q_gamma_trunc(double p, double shape, double rate,
	double a, double b, bool lower, bool log)
{
	double scale = 1 / rate;
	double lp = log ? p : std::log(p);

	double lp_a = R::pgamma(a, shape, scale, true, true);
	double lp_b = R::pgamma(b, shape, scale, true, true);
	double lp_diff = log_sub2_exp(lp_b, lp_a);

	double clp_a = R::pgamma(a, shape, scale, false, true);
	double clp_b = R::pgamma(b, shape, scale, false, true);
	double clp_diff = log_sub2_exp(clp_a, clp_b);

	double out;

	if (lp_diff < clp_diff) {
		double lp_adj = log_add2_exp(lp_a, lp + log_sub2_exp(lp_b, lp_a));
		out = R::qgamma(lp_adj, shape, scale, true, true);
	} else {
		double clp_adj = log_sub2_exp(clp_a, lp + log_sub2_exp(clp_a, clp_b));
		out = R::qgamma(clp_adj, shape, scale, false, true);
	}

    return out;
}

/*
* Draw from Truncated Gamma distribution
*
* - `n`: number of desired draws.
* - `shape`: shape parameter.
* - `rate`: rate parameter.
* - `a`: lower limit of support.
* - `b`: upper limit of support.
*
* Returns vector of $n$ iid draws.
*/
inline Rcpp::NumericVector r_gamma_trunc(unsigned int n, double shape, double rate,
	double a, double b)
{
	Rcpp::NumericVector out(n);

	const Rcpp::NumericVector& u = Rcpp::runif(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_gamma_trunc(u(i), shape, rate, a, b, true, false);
	}

	return out;
}

/*
* Density of Truncated Gamma distribution
*
* - `x`: vector of $n$ density arguments.
* - `shape`: shape parameter.
* - `rate`: rate parameter.
* - `a`: lower limit of support.
* - `b`: upper limit of support.
* - `log`: if `true`, return log of value. Otherwise return value on the
*   original scale.
*
* Returns vector of $n$ values.
*/
inline Rcpp::NumericVector d_gamma_trunc(const Rcpp::NumericVector& x, double shape,
	double rate, double a, double b, bool log)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);

	const Rcpp::NumericVector& u = Rcpp::runif(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_gamma_trunc(x(i), shape, rate, a, b, log);
	}

	return out;
}

/*
* CDF of Truncated Gamma distribution
*
* - `x`: vector of $n$ CDF arguments.
* - `shape`: shape parameter.
* - `rate`: rate parameter.
* - `a`: lower limit of support.
* - `b`: upper limit of support.
* - `log`: if `true`, return log of value. Otherwise return value on the
*   original scale.
*
* Returns vector of $n$ values.
*/
inline Rcpp::NumericVector p_gamma_trunc(const Rcpp::NumericVector& q, double shape,
	double rate, double a, double b, bool lower, bool log)
{
	unsigned int n = q.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = p_gamma_trunc(q(i), shape, rate, a, b, lower, log);
	}

	return out;
}

/*
* Quantile function of Truncated Gamma distribution
*
* - `p`: vector of $n$ quantile arguments.
* - `shape`: shape parameter.
* - `rate`: rate parameter.
* - `a`: lower limit of support.
* - `b`: upper limit of support.
* - `lower`: if `true`, request $p$ quantile. Otherwise request $1-p$ quantile.
* - `log`: if `true`, assume $p$ is specified on the log-scale. Otherwise,
*    assume it is on the original scale.
*
* Returns vector of $n$ values.
*/
inline Rcpp::NumericVector q_gamma_trunc(const Rcpp::NumericVector& p, double shape,
	double rate, double a, double b, bool lower, bool log)
{
	unsigned int n = p.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_gamma_trunc(p(i), shape, rate, a, b, lower, log);
	}

	return out;
}

/*
* Draw from Truncated Gamma distribution
*
* - `shape`: shape parameter.
* - `rate`: rate parameter.
* - `a`: lower limit of support.
* - `b`: upper limit of support.
*
* Returns vector of $n$ iid draws.
*/
inline double r_gamma_trunc(double shape, double rate, double a, double b)
{
	double u = R::runif(0, 1);
	return q_gamma_trunc(u, shape, rate, a, b, true, false);
}

}

#endif
