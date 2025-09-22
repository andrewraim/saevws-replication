#include "rcpp-categ.h"
#include "vws.h"

Rcpp::IntegerVector r_categ_rcpp(unsigned int n, const Rcpp::NumericVector& p,
	bool log, bool one_based)
{
	return vws::r_categ(n, p, log, one_based);
}



