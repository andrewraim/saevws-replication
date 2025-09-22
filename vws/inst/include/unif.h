#ifndef VWS_UNIF_H
#define VWS_UNIF_H

#include <Rcpp.h>
#include "log-sum-exp.h"

/*
* Functions for the Uniform distribution with lower limit `lo` and upper limit
* `hi`.
*/

namespace vws {

/*
* Moment-generating function for Uniform distribution
*
* - `s`: argument to MGF.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `log`: if `true`, return value on the log-scale. Otherwise, return on the
*   original scale.
*
* Returns a single MGF value.
*/
inline double mgf_unif(double s, double lo, double hi, bool log)
{
	double out;

	if (s < 0) {
		out = log_sub2_exp(lo*s, hi*s) - std::log(hi - lo) - std::log(-s);
	} else if (s > 0) {
		out = log_sub2_exp(hi*s, lo*s) - std::log(hi - lo) - std::log(s);
	} else {
		out = R_NegInf;
	}

	return log ? out : std::exp(out);
}

/*
* Moment-generating function for Uniform distribution
*
* - `s`: a vector of $n$ arguments to MGF.
* - `lo`: lower limit parameter.
* - `hi`: upper limit parameter.
* - `log`: if `true`, return values on the log-scale. Otherwise, return on the
*   original scale.
*
* Returns a vector of $n$ values.
*/
inline Rcpp::NumericVector mgf_unif(const Rcpp::NumericVector& s, double lo,
	double hi, bool log = false)
{
	unsigned int n = s.length();
	Rcpp::NumericVector out(n);

	for (unsigned int i = 0; i < n; i++) {
		out(i) = mgf_unif(s(i), lo, hi, log);
	}

	return out;
}

}

#endif
