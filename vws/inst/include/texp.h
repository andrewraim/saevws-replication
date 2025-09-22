#ifndef VWS_TEXP_H
#define VWS_TEXP_H

#include <Rcpp.h>
#include "log-sum-exp.h"

/*
* Functions for the truncated exponential distribution, with lower limit `lo`
* upper limit `hi`, and rate parameter `kappa`.
*/

namespace vws {

/*
* Normalizing constant of Truncated Exponential distribution
*
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `log`: if `true`, return value on the log-scale. Otherwise return on the
*   original scale
*
* Returns value of normalizing constant.
*/
inline double n_texp(double kappa, double lo, double hi, bool log = false)
{
	double out;

	if (kappa > 0) {
		out = log_sub2_exp(kappa*hi, lo*kappa) - std::log(kappa);
	} else {
		out = log_sub2_exp(lo*kappa, hi*kappa) - std::log(-kappa);
	}

	return log ? out : std::exp(out);
}

/*
* Density function of Truncated Exponential distribution
*
* - `x`: density argument.
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `log`: if `true`, return value on the log-scale. Otherwise return on the
*   original scale
*
* Returns value of density.
*/
inline double d_texp(double x, double kappa, double lo, double hi,
	bool log = false)
{
	double out = kappa*x - n_texp(kappa, lo, hi, true) + std::log(lo <= x && x <= hi);
	return log ? out : std::exp(out);
}

/*
* CDF of Truncated Exponential distribution
*
* - `q`: CDF argument.
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `lower`: if `true`, compute $P(X \leq q)$. Otherwise compute $P(X > q)$.
* - `log`: if `true`, return value on the log-scale. Otherwise, return it on
*   the original scale.
*
* Returns value of CDF.
*/
inline double p_texp(double q, double kappa, double lo, double hi,
	bool lower = true, bool log = false)
{
	double out;

	if (q < lo) {
		out = R_NegInf;
	} else if (q > hi) {
		out = 0;
	} else if (kappa < 0) {
		out = vws::log_sub2_exp(lo*kappa, kappa*q) - vws::log_sub2_exp(lo*kappa, hi*kappa);
	} else {
		out = vws::log_sub2_exp(kappa*q, lo*kappa) - vws::log_sub2_exp(hi*kappa, lo*kappa);
	}

	if (!lower) {
		out = vws::log_sub2_exp(0, out);
	}

	return log ? out : std::exp(out);
}

/*
* Quantile function of Truncated Exponential distribution
*
* - `p`: quantile argument.
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `lower`: if `true`, request $p$ quantile. Otherwise request $1-p$ quantile.
* - `log`: if `true`, assume $p$ is specified on the log-scale. Otherwise,
*   assume it is on the original scale.
*
* Returns value of quantile.
*/
inline double q_texp(double p, double kappa, double lo, double hi,
	bool lower = true, bool log = false)
{
	double lp = log ? p : std::log(p);
	lp = lower ? lp : vws::log_sub2_exp(0, lp);

	double out;

	if (lp > 0 || std::isnan(p)) {
		out = R_NaN;
	} else if (kappa < 0) {
		out = (1 / kappa) * vws::log_sub2_exp(lo*kappa, lp + vws::log_sub2_exp(lo*kappa, hi*kappa));
	} else {
		out = (1 / kappa) * vws::log_add2_exp(lo*kappa, lp + vws::log_sub2_exp(hi*kappa, lo*kappa));
	}

	return out;
}

/*
* Draw from Truncated Exponential distribution
*
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
*
* Returns a draw.
*/
inline double r_texp(double kappa, double lo, double hi)
{
	double u = R::runif(0, 1);
	return q_texp(u, kappa, lo, hi);
}

/*
* Moment-generating function of Truncated Exponential distribution
*
* - `s`: MGF argument.
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `log`: if `true`, return value on the log-scale. Otherwise, return it on
*   the original scale.
*
* Returns value of MGF.
*/
inline double mgf_texp(double s, double kappa, double lo, double hi, bool log)
{
	double out;
	double sk = s + kappa;

	if (sk < 0) {
		out = std::log(kappa) - std::log(-sk) +
			vws::log_sub2_exp(lo*sk, hi*sk) -
			vws::log_sub2_exp(hi*kappa, lo*kappa);
	} else {
		out = std::log(kappa) - std::log(sk) +
			vws::log_sub2_exp(hi*sk, lo*sk) -
			vws::log_sub2_exp(hi*kappa, lo*kappa);
	}

	return log ? out : std::exp(out);
}

/*
* Integrate density of Truncated Exponential distribution
*
* - `a`: lower limit of integral.
* - `b`: upper limit of integral.
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `log`: if `true`, return value on the log-scale. Otherwise, return it on
*   the original scale.
*
* Returns integral, from $a$ to $b$, under the TExp distribution. This
* avoids loss of precision that can happen if we take the difference of two
* `p_texp` calls. (Even if those are computed on the log-scale, there can be
* a loss of precision which may be avoided here).
*/
inline double integrate_texp(double a, double b, double kappa, double lo,
	double hi, bool log = false)
{
	double out;

	if (a < lo && hi < b) {
		out = 0;
	} else if (a < lo) {
		out = p_texp(b, kappa, lo, hi, true, true);
	} else if (hi < b) {
		out = p_texp(a, kappa, lo, hi, false, true);
	} else if (kappa < 0) {
		double lnum = vws::log_sub2_exp(kappa*a, kappa*b);
		double lden = vws::log_sub2_exp(lo*kappa, hi*kappa);
		out = lnum - lden;
	} else {
		double lnum = vws::log_sub2_exp(kappa*b, kappa*a);
		double lden = vws::log_sub2_exp(hi*kappa, lo*kappa);
		out = lnum - lden;
	}

	return log ? out : std::exp(out);
}

/*
* Density function for Truncated Exponential distribution
*
* - `x`: a vector of $n$ density arguments.
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `log`: if `true`, return values on the log-scale. Otherwise, return them
*   on the original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector d_texp(const Rcpp::NumericVector& x, double kappa,
	double lo, double hi, bool log = false)
{
	unsigned int n = x.length();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_texp(x(i), kappa, lo, hi, log);
	}

	return out;
}

/*
* CDF for Truncated Exponential distribution
*
* - `q`: a vector of $n$ CDF arguments.
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `lower`: if `true`, compute $P(X \leq q)$. Otherwise compute $P(X > q)$.
* - `log`: if `true`, return values on the log-scale. Otherwise, return them on
*   the original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector p_texp(const Rcpp::NumericVector& q, double kappa,
	double lo, double hi, bool lower = true, bool log = false)
{
	unsigned int n = q.length();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = p_texp(q(i), kappa, lo, hi, lower, log);
	}

	return out;
}

/*
* Quantile function for Truncated Exponential distribution
*
* - `p`: a vector of $n$ quantile arguments.
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `lower`: if `true`, request $p$ quantile. Otherwise request $1-p$ quantile.
* - `log`: if `true`, assume $p$ is specified on the log-scale. Otherwise,
*    assume it is on the original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector q_texp(const Rcpp::NumericVector& p, double kappa,
	double lo, double hi, bool lower = true, bool log = false)
{
	unsigned int n = p.length();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_texp(p(i), kappa, lo, hi, lower, log);
	}

	return out;
}

/*
* Draw from Truncated Exponential distribution
*
* - `n`: desired number of draws.
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
*
* Returns a vector of $n$ iid draws.
*/
inline Rcpp::NumericVector r_texp(unsigned int n, double kappa, double lo,
	double hi)
{
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	return q_texp(u, kappa, lo, hi);
}

/*
* Moment-generating function of Truncated Exponential distribution
*
* - `s`: a vector of $n$ arguments to MGF.
* - `kappa`: rate parameter.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `log`: if `true`, return values on the log-scale. Otherwise, return them
*   on the original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector mgf_texp(const Rcpp::NumericVector& s, double kappa,
	double lo, double hi, bool log = false)
{
	unsigned int n = s.length();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = mgf_texp(s(i), kappa, lo, hi, log);
	}

	return out;
}

}

#endif
