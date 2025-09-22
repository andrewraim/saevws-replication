#include "rcpp-rect.h"
#include "vws.h"

Rcpp::NumericVector rect_rcpp(const Rcpp::NumericVector& z,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b)
{
	return vws::rect(z, a, b);
}

Rcpp::NumericVector inv_rect_rcpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b)
{
	return vws::inv_rect(x, a, b);
}

