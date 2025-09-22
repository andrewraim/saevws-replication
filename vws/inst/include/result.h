#ifndef VWS_RESULT_H
#define VWS_RESULT_H

#include <RcppCommon.h>
#include <vector>

namespace vws {

/*
* Structure for results of `optimize_hybrid` function.
*
* - `par`: value of the optimization variable.
* - `report`: value of the objective function.
* - `method`: description of the method used to find the result.
* - `status`: Corresponds to a code from BFGS if this is used as `method`;
*   otherwise zero.
*
* The `method` field can take on the following values:
*
* - "Brent": Brent optimization method was used.
* - "BFGS": BFGS method was used.
* - "Lower Limit Inf": For a maximization problem, the lower limit was taken
*   as `par` because it had value `inf`.
* - "Upper Limit Inf": For a maximization problem, the upper limit was taken
*   as `par` because it had value `inf`.
* - "Lower Limit NegInf": For a minimization problem, the lower limit was
*   taken as `par` because it had value `-inf`.
* - "Upper Limit NegInf": For a minimization problem, the upper limit was
*   taken as `par` because it had value `-inf`.
* - "Max at Lower Limit": numerical maximization was used, but a larger value
*   was found at the lower limit.
* - "Max at Upper Limit": numerical maximization was used, but a larger value
*   was found at the upper limit.
* - "Min at Lower Limit": numerical minimization was used, but a smaller value
*   was found at the lower limit.
* - "Min at Upper Limit": numerical minimization was used, but a smaller value
*   was found at the upper limit.
*/
struct optimize_hybrid_result {
	double par;
	double value;
	std::string method;
	int status;

	operator SEXP() const;
};

/*
* Structure for results of `rejection` function.
*
*  - `draws`: vector of draws.
*  - `rejects`: vector of rejection counts. The $i$th element contains the
*    count of rejections to obtain the $i$th element of `draws`.
*/
template <typename T>
struct rejection_result
{
	std::vector<T> draws;
	std::vector<unsigned int> rejects;

	operator SEXP() const;
};

}

/*
* SEXP operators below allow each structure above to be converted to an
* `Rcpp::List` using `Rcpp::wrap(...)`.
*/

#include <Rcpp.h>

namespace vws {

inline optimize_hybrid_result::operator SEXP() const
{
	return Rcpp::List::create(
		Rcpp::Named("par") = par,
		Rcpp::Named("value") = value,
		Rcpp::Named("method") = method,
		Rcpp::Named("status") = status
	);
}

template <typename T>
inline rejection_result<T>::operator SEXP() const
{
	return Rcpp::List::create(
		Rcpp::Named("draws") = draws,
		Rcpp::Named("rejects") = rejects
	);
}

}

#endif
