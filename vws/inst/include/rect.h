#ifndef VWS_RECT_H
#define VWS_RECT_H

#include <Rcpp.h>
#include "logit.h"

/*
* Functions for a "rectangular" transformation and its inverse.
*/


namespace vws {

/*
* Rectangular transformation
*
* - `z`: scalar in the rectangle $[a, b]$.
* - `a`: lower limit which may be `-inf`.
* - `b`: upper limit which may be `+inf`.
*
* Returns a scalar on the real line.
*/
inline double rect(double z, double a, double b)
{
	double x;

	if (std::isinf(a) && std::isinf(b) && a < 0 && b > 0) {
		x = z;
	} else if (std::isinf(a) && a < 0) {
		x = -log(b - z);
	} else if (std::isinf(b) && b > 0) {
		x = log(z - a);
	} else if (std::isfinite(a) && std::isfinite(b)) {
		x = logit((z - a) / (b - a));
	} else {
		Rcpp::stop("Invalid endpoints");
	}

	return x;
}

/*
* Inverse rectangular transformation
*
* - `z`: scalar in the real line.
* - `a`: lower limit which may be `-inf`.
* - `b`: upper limit which may be `+inf`.
*
* Returns a scalar in the rectangle $[a, b]$.
*/
inline double inv_rect(double x, double a, double b)
{
	double z;

	if (std::isinf(a) && std::isinf(b) && a < 0 && b > 0) {
		z = x;
	} else if (std::isinf(a) && a < 0) {
		z = b - exp(-x);
	} else if (std::isinf(b) && b > 0) {
		z = a + exp(x);
	} else if (std::isfinite(a) && std::isfinite(b)) {
		z = (b - a) * inv_logit(x) + a;
	} else {
		Rcpp::stop("Invalid endpoints");
	}

	return z;
}

/*
* Rectangular transformation
*
* - `z`: a vector in the rectangle $[a_1, b_1] \times \cdots \times [a_d, b_d]$.
* - `a`: a vector of $d$ lower limits which each may be `-inf`.
* - `b`: a vector of $d$ upper limits which each may be `+inf`.
*
* Returns a vector of $d$ elements on the real line.
*/
inline Rcpp::NumericVector rect(const Rcpp::NumericVector& z,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b)
{
	unsigned int n = z.size();
	if (n != a.size()) {
		Rcpp::stop("Dimension mismatch");
	}
	if (n != b.size()) {
		Rcpp::stop("Dimension mismatch");
	}
	if (Rcpp::is_true(Rcpp::any(a > b))) {
		Rcpp::stop("any(a > b)");
	}

	Rcpp::NumericVector x(n);

	for (unsigned int i = 0; i < n; i++) {
		x(i) = rect(z(i), a(i), b(i));
	}

	return x;
}

/*
* Inverse rectangular transformation
*
* - `z`: a vector of $d$ elements on the real line.
* - `a`: a vector of $d$ lower limits which each may be `-inf`.
* - `b`: a vector of $d$ upper limits which each may be `+inf`.
*
* Returns a vector in the rectangle $[a_1, b_1] \times \cdots \times [a_d, b_d]$.
*/
inline Rcpp::NumericVector inv_rect(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b)
{
	unsigned int n = x.size();
	if (n != a.size()) {
		Rcpp::stop("Dimension mismatch");
	}
	if (n != b.size()) {
		Rcpp::stop("Dimension mismatch");
	}
	if (is_true(Rcpp::any(a > b))) {
		Rcpp::stop("any(a > b)");
	}

	Rcpp::NumericVector z(n);

	for (unsigned int i = 0; i < n; i++) {
		z(i) = inv_rect(x(i), a(i), b(i));
	}

	return z;
}

}

#endif
