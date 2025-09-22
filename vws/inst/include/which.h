#ifndef VWS_WHICH_H
#define VWS_WHICH_H

#include <Rcpp.h>

namespace vws {

/*
* A C++ implementation of the R `which` function.
*
* - `x`: a vector of $n$ logical values.
* - `one_based`: if `true`, returned indices will be in $\{ 1, \ldots, n \}$.
*   Otherwise, returned indices will be in $\{ 0, \ldots, n-1 \}$.
*
* Returns a vector of indices whose length is between $0$ and $n$. Indices
* correspond to elements of `x` which are `true`.
*/
inline Rcpp::IntegerVector which(const Rcpp::LogicalVector& x, bool one_based = false)
{
	std::vector<unsigned int> idx;

	for (unsigned int i = 0; i < x.size(); i++) {
		if (x(i)) {
			idx.push_back(i + one_based);
		}
	}

	return Rcpp::IntegerVector(idx.begin(), idx.end());
}

}

#endif
