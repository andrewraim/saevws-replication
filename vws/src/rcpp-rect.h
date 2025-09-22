#ifndef RCPP_RECT_H
#define RCPP_RECT_H

#include <Rcpp.h>

//' Rectangular transformation
//'
//' A transformation from Euclidean to a rectangle in \eqn{\mathbb{R}^n}, and
//' its inverse transformation.
//'
//' @param x A point in \eqn{\mathbb{R}^{d}}.
//' @param z A point in the rectangle \eqn{[a_1,b_1] \times \cdots \times [a_d,b_d]}.
//' @param a A vector \eqn{(a_1, \ldots, a_d)}, Elements may be `-Inf`.
//' @param b A vector \eqn{(b_1, \ldots, b_d)}, Elements may be `+Inf`.
//'
//' @examples
//' n = 20
//' x = seq(-5, 5, length.out = n)
//'
//' # Transform x to the interval [-1, 1]
//' a = rep(-1, n)
//' b = rep(+1, n)
//' z = inv_rect(x, a, b)
//' print(z)
//' xx = rect(z, a, b)
//' stopifnot(all(abs(x - xx) < 1e-8))
//'
//' # Transform x to the interval [-Inf, 1]
//' a = rep(-Inf, n)
//' b = rep(+1, n)
//' z = inv_rect(x, a, b)
//' print(z)
//' xx = rect(z, a, b)
//' stopifnot(all(abs(x - xx) < 1e-8))
//'
//' # Transform x to the interval [-1, Inf]
//' a = rep(-1, n)
//' b = rep(+Inf, n)
//' z = inv_rect(x, a, b)
//' print(z)
//' xx = rect(z, a, b)
//' stopifnot(all(abs(x - xx) < 1e-8))
//'
//' # Transform x to the interval [-Inf, Inf]
//' a = rep(-Inf, n)
//' b = rep(+Inf, n)
//' z = inv_rect(x, a, b)
//' print(z)
//' xx = rect(z, a, b)
//' stopifnot(all(abs(x - xx) < 1e-8))
//'
//' @name rect
//' @export
// [[Rcpp::export(name = "rect")]]
Rcpp::NumericVector rect_rcpp(const Rcpp::NumericVector& z,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b);

//' @name rect
//' @export
// [[Rcpp::export(name = "inv_rect")]]
Rcpp::NumericVector inv_rect_rcpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b);

#endif

