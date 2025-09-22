#ifndef VWS_CATEG_H
#define VWS_CATEG_H

#include <Rcpp.h>
#include "gumbel.h"

/*
* Functions for the categorical distribution on $\{0, \ldots, k-1\}$ with
* probabilities $\{p_1, \ldots, p_k\}$.
*/

namespace vws {

/*
* Draw from categorical distribution
*
* - `p`: a vector of $k$ probabilities.
* - `log`: if `true`, assume probabilities are given on the log-scale.
*   Otherwise, assume they are on the original scale.
* - `one_based`: if `true`, assume the support of the distribution is
*   $\{1, \ldots, k\}$. Otherwise assume it is $\{0, \ldots, k-1\}$.
*
* Returns a draw.
*/
inline unsigned int r_categ(const Rcpp::NumericVector& p, bool log = false,
	bool one_based = false)
{
	unsigned int k = p.length();
	Rcpp::NumericVector lp;
	if (log) {
		lp = p;
	} else {
		lp = Rcpp::log(p);
	}

	const Rcpp::NumericVector& z = r_gumbel(k);
	return Rcpp::which_max(z + lp) + one_based;
}

/*
* Draw from categorical distribution
*
* - `n`: number of desired draws.
* - `p`: a vector of $k$ probabilities.
* - `log`: if `true`, assume probabilities are given on the log-scale.
*    Otherwise, assume they are on the original scale.
* - `one_based`: if `true`, assume the support of the distribution is
*   $\{1, \ldots, k\}$. Otherwise assume it is $\{0, \ldots, k-1\}$.
*
* Returns a vector of $n$ iid draws.
*/
inline Rcpp::IntegerVector r_categ(unsigned int n,
	const Rcpp::NumericVector& p, bool log = false, bool one_based = false)
{
	Rcpp::IntegerVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = r_categ(p, log, one_based);
	}

	return out;
}

}

#endif

