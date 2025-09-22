#include "rcpp-unif.h"
#include "vws.h"

Rcpp::NumericVector mgf_unif_rcpp(const Rcpp::NumericVector& s, double lo,
	double hi, bool log)
{
	return vws::mgf_unif(s, lo, hi, log);
}

