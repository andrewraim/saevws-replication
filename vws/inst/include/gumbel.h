#ifndef VWS_GUMBEL_H
#define VWS_GUMBEL_H

#include <Rcpp.h>

/*
* Functions for the Gumbel distribution with location parameter $\mu$ and
* scale parameter $\sigma$.
*/

namespace vws {

/*
* Quantile function for Gumbel distribution
*
* - `p`: probability of desired quantile.
* - `mu`: location parameter.
* - `sigma`: scale parameter.
* - `lower`: if `true`, request $p$ quantile. Otherwise request $1-p$ quantile.
* - `log`: if `true`, assume $p$ is specified on the log-scale. Otherwise,
*   assume it is on the original scale.
*
* Returns a single value.
*/
inline double q_gumbel(double p, double mu = 0, double sigma = 1,
	bool lower = true, bool log = false)
{
	double lp0 = log ? p : std::log(p);
	double lp = lower ? lp0 : std::log1p(-std::exp(lp0));
	return mu - sigma * std::log(-lp);
}

/*
* Draw from Gumbel distribution
*
* - `mu`: location parameter.
* - `sigma`: scale parameter.
*
* Returns a draw.
*/
inline double r_gumbel(double mu = 0, double sigma = 1)
{
	double u = R::runif(0, 1);
	return q_gumbel(u, mu, sigma);
}

/*
* Density function for Gumbel distribution
*
* - `x`: argument of density.
* - `mu`: location parameter.
* - `sigma`: scale parameter.
* - `log`: if `true`, return density on the log-scale. Otherwise, return it on
*    the original scale.
*
* Returns a single value.
*/
inline double d_gumbel(double x, double mu = 0, double sigma = 1,
	bool log = false)
{
	double z = (x - mu) / sigma;
	double out = -std::log(sigma) - (z + std::exp(-z));
	return log ? out : std::exp(out);
}

/*
* CDF for Gumbel distribution
*
* - `q`: argument of CDF.
* - `mu`: location parameter.
* - `sigma`: scale parameter.
* - `lower`: if `true`, compute $P(X \leq q)$. Otherwise compute $P(X > q)$.
* - `log`: if `true`, return value on the log-scale. Otherwise, return it on
*   the original scale.
*
* Returns a single value.
*/
inline double p_gumbel(double q, double mu = 0, double sigma = 1,
	bool lower = true, bool log = false)
{
	double z = (q - mu) / sigma;
	double out0 = -std::exp(-z);
	double out = lower ? out0 : std::log1p(-std::exp(out0));
	return log ? out : std::exp(out);
}

/*
* Draw from Gumbel distribution
*
* - `n`: desired number of draws.
* - `mu`: location parameter.
* - `sigma`: scale parameter.
*
* Returns a vector of $n$ iid draws.
*/
inline Rcpp::NumericVector r_gumbel(unsigned int n, double mu = 0,
	double sigma = 1)
{
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_gumbel(u(i), mu, sigma);
	}

	return out;
}

/*
* Density function for Gumbel distribution
*
* - `x`: a vector of $n$ density arguments.
* - `mu`: location parameter.
* - `sigma`: scale parameter.
* - `log`: if `true`, return density values on the log-scale. Otherwise, return
*    them on the original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector d_gumbel(const Rcpp::NumericVector& x, double mu = 0,
	double sigma = 1, bool log = false)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = d_gumbel(x(i), mu, sigma, log);
	}

	return out;
}

/*
* CDF for Gumbel distribution
*
* - `q`: a vector of $n$ CDF arguments.
* - `mu`: location parameter.
* - `sigma`: scale parameter.
* - `lower`: if `true`, compute $P(X \leq q)$. Otherwise compute $P(X > q)$.
* - `log`: if `true`, return values on the log-scale. Otherwise, return them on
*   the original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector p_gumbel(const Rcpp::NumericVector& q, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false)
{
	unsigned int n = q.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = p_gumbel(q(i), mu, sigma, lower, log);
	}

	return out;
}

/*
* Quantile function for Gumbel distribution
*
* - `p`: a vector of $n$ quantile arguments.
* - `mu`: location parameter.
* - `scale`: scale parameter.
* - `lower`: if `true`, request $p$ quantile. Otherwise request $1-p$ quantile.
* - `log`: if `true`, assume $p$ is specified on the log-scale. Otherwise,
*    assume it is on the original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector q_gumbel(const Rcpp::NumericVector& p, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false)
{
	unsigned int n = p.size();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_gumbel(p(i), mu, sigma, lower, log);
	}

	return out;
}

}

#endif

