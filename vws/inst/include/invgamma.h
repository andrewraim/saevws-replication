#ifndef VWS_INVGAMMA_H
#define VWS_INVGAMMA_H

#include <Rcpp.h>
#include "log-sum-exp.h"

/*
* Functions for the Inverse Gamma distribution with shape parameter $a$ and
* rate parameter $b$.
*/

namespace vws {

/*
* Draw from Inverse Gamma distribution
*
* - `a`: shape parameter.
* - `b`: rate parameter.
*
* Returns a draw.
*/
inline double r_invgamma(double a, double b)
{
	return 1 / R::rgamma(a, 1 / b);
}

/*
* Density function of Inverse Gamma distribution
*
* - `x`: argument of density.
* - `a`: shape parameter.
* - `b`: rate parameter.
* - `log`: if `true`, return log of density value. Otherwise return value on
*   the original scale.
*
* Returns a single value.
*/
inline double d_invgamma(double x, double a, double b, bool log = false)
{
	double out = R::dgamma(1 / x, a, 1 / b, true) - 2 * std::log(x);
	return log ? out : exp(out);
}

/*
* CDF of Inverse Gamma distribution
*
* - `x`: argument of CDF.
* - `a`: shape parameter.
* - `b`: rate parameter.
* - `lower`: if `true`, compute $P(X \leq q)$. Otherwise compute $P(X > q)$.
* - `log`: if `true`, return value on the log-scale. Otherwise, return it on
*   the original scale.
*
* Returns a single value.
*/
inline double p_invgamma(double q, double a, double b, bool lower = true,
	bool log = false)
{
	return R::pgamma(1 / q, a, 1 / b, !lower, log);
}

/*
* Quantile function for Inverse Gamma distribution
*
* - `p`: probability of desired quantile.
* - `a`: shape parameter.
* - `b`: rate parameter.
* - `lower`: if `true`, request $p$ quantile. Otherwise request $1-p$ quantile.
* - `log`: if `true`, assume $p$ is specified on the log-scale. Otherwise,
*    assume it is on the original scale.
*
* Returns a single value.
*/
inline double q_invgamma(double p, double a, double b, bool lower = true,
	bool log = false)
{
	if (!log) { p = std::log(p); }

	double out;
	if (p > std::log(1/2)) {
		double cp = log_sub2_exp(0, p);
		out = 1 / R::qgamma(cp, a, 1 / b, lower, true);
	} else {
		out = 1 / R::qgamma(p, a, 1 / b, !lower, true);
	}

	return out;
}

/*
* Draw from Inverse Gamma distribution
*
* - `n`: desired number of draws.
* - `a`: shape parameter.
* - `b`: rate parameter.
*
* Returns a vector of $n$ iid draws.
*/
inline Rcpp::NumericVector r_invgamma(unsigned int n, double a, double b)
{
	return 1 / Rcpp::rgamma(n, a, 1 / b);
}

/*
* Density function for Inverse Gamma distribution
*
* - `x`: a vector of $n$ density arguments.
* - `a`: shape parameter.
* - `b`: rate parameter.
* - `log`: if `true`, return density values on the log-scale. Otherwise, return
*    them on the original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector d_invgamma(const Rcpp::NumericVector& x, double a,
	double b, bool log)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = vws::d_invgamma(x(i), a, b, log);
	}

	return out;
}

/*
* CDF for Inverse Gamma distribution
*
* - `q`: a vector of $n$ CDF arguments.
* - `a`: shape parameter.
* - `b`: rate parameter.
* - `lower`: if `true`, compute $P(X \leq q)$. Otherwise compute $P(X > q)$.
* - `log`: if `true`, return values on the log-scale. Otherwise, return them on
*   the original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector p_invgamma(const Rcpp::NumericVector& q, double a,
	double b, bool lower, bool log)
{
	unsigned int n = q.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = vws::p_invgamma(q(i), a, b, lower, log);
	}

	return out;
}

/*
* Quantile function for Gumbel distribution
*
* - `p`: a vector of $n$ quantile arguments.
* - `a`: shape parameter.
* - `b`: rate parameter.
* - `lower`: if `true`, request $p$ quantile. Otherwise request $1-p$ quantile.
* - `log`: if `true`, assume $p$ is specified on the log-scale. Otherwise,
*    assume it is on the original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector q_invgamma(const Rcpp::NumericVector& p, double a,
	double b, bool lower, bool log)
{
	unsigned int n = p.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = vws::q_invgamma(p(i), a, b, lower, log);
	}

	return out;
}

}

#endif

