#ifndef VWS_UNIVARIATE_HELPER_H
#define VWS_UNIVARIATE_HELPER_H

#include <Rcpp.h>

namespace vws {

/*
* This abstract class encapsulates several functions for univariate
* distributions. It is used with the UnivariateConstRegion class.
*
*
*
*/
template <class T>
class UnivariateHelper
{
public:

	/*
	* Density function
	* - `x`: argument
	* - `log`: if `true`, return density values on the log-scale. Otherwise,
	*   return them on the original scale.
	*/
	virtual double d(T x, bool log = false) const = 0;

	/*
	* CDF
	* - `q`: argument
	* - `lower`: if `true`, compute $P(X \leq q)$. Otherwise compute $P(X > q)$.
	* - `log`: if `true`, return values on the log-scale. Otherwise, return
	*   them on the original scale.
	*/
	virtual double p(T q, bool lower = true, bool log = false) const = 0;

	/*
	* Quantile function
	* - `p`: argument
	* - `lower`: if `true`, request $p$ quantile. Otherwise request $1-p$
	*   quantile.
	* - `log`: if `true`, assume $p$ is specified on the log-scale. Otherwise,
	*   assume it is on the original scale.
	*/
	virtual double q(T p, bool lower = true, bool log = false) const = 0;

	/*
	* Support indicator function
	* - `x`: argument
	* Returns `true` if $x$ is in the support of the distribution; `false`
	* otherwise.
	*/
	virtual bool s(T x) const = 0;
};

}

#endif
