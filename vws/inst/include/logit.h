#ifndef VWS_LOGIT_H
#define VWS_LOGIT_H

#include <Rcpp.h>

/*
* Functions for the logit transformation and its inverse.
*/

namespace vws {

/*
* Logit transformation of $p in [0,1]$.
*/
inline double logit(double p)
{
	return R::qlogis(p, 0, 1, true, false);
}

/*
* Inverse logit transformation of $x in R$.
*/
inline double inv_logit(double x)
{
	return R::plogis(x, 0, 1, true, false);
}

/*
* Elementwise logit transformation from $p \in [0,1]^n$.
*/
inline Rcpp::NumericVector logit(const Rcpp::NumericVector& p)
{
	return Rcpp::qlogis(p);
}

/*
* Elementwise inverse logit transformation from $x \in R^n$.
*/
inline Rcpp::NumericVector inv_logit(const Rcpp::NumericVector& x)
{
	return Rcpp::plogis(x);
}

}

#endif
