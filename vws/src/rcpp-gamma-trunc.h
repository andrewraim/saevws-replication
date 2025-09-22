#ifndef RCPP_GAMMA_TRUNC_H
#define RCPP_GAMMA_TRUNC_H

#include <Rcpp.h>

//' Truncated Gamma Distribution
//'
//' Functions for the Truncated Gamma distribution parameterized by shape and
//' rate, whose density is
//' \deqn{
//' f(x \mid a, b) =
//' \frac{
//' b^a
//' }{
//' \Gamma(a)
//' }
//' x^{-a-1} e^{-b / x} \textrm{I}(x \geq 0)
//' }
//'
//' @param n Number of draws.
//' @param x Vector; argument of density.
//' @param q Vector; argument of cumulative distribution function.
//' @param p Vector; argument of quantile function.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @param a lower limit of support.
//' @param b upper limit of support.
//' @param lower logical; if `TRUE` (default), probabilities are
//' \eqn{P(X \leq x)}; otherwise, \eqn{P(X > x)}.
//' @param log If `TRUE`, return densities and probabilities on the log-scale.
//'
//' @return
//' `d_gamma_trunc` computes the density, `r_gamma_trunc` generates random
//' deviates, `p_gamma_trunc` computes the CDF, and `q_gamma_trunc` computes
//' quantiles.
//'
//' @examples
//' shape = 10
//' rate = 5
//' a = 2
//' b = 5
//' x = r_gamma_trunc(100000, shape, rate, a, b)
//' xx = seq(a, b, length.out = 100)
//'
//' plot(density(x))
//' lines(xx, d_gamma_trunc(xx, shape, rate, a, b), lty = 2, col = "blue", lwd = 2)
//'
//' plot(ecdf(x))
//' lines(xx, p_gamma_trunc(xx, shape, rate, a, b), lty = 2, col = "blue", lwd = 2)
//'
//' pp = seq(0, 1, length.out = 102) |> head(-1) |> tail(-1)
//' qq = quantile(x, probs = pp)
//' plot(pp, qq)
//' lines(pp, q_gamma_trunc(pp, shape, rate, a, b), lty = 2, col = "blue", lwd = 2)
//'
//' @name GammaTrunc
//' @export
// [[Rcpp::export(name = "r_gamma_trunc")]]
Rcpp::NumericVector r_gamma_trunc_rcpp(unsigned int n, double shape,
	double rate, double a, double b);

//' @name GammaTrunc
//' @export
// [[Rcpp::export(name = "d_gamma_trunc")]]
Rcpp::NumericVector d_gamma_trunc_rcpp(const Rcpp::NumericVector& x,
	double shape, double rate, double a, double b, bool log = false);

//' @name GammaTrunc
//' @export
// [[Rcpp::export(name = "p_gamma_trunc")]]
Rcpp::NumericVector p_gamma_trunc_rcpp(const Rcpp::NumericVector& q,
	double shape, double rate, double a, double b, bool lower = true,
	bool log = false);

//' @name GammaTrunc
//' @export
// [[Rcpp::export(name = "q_gamma_trunc")]]
Rcpp::NumericVector q_gamma_trunc_rcpp(const Rcpp::NumericVector& p,
	double shape, double rate, double a,double b, bool lower = true,
	bool log = false);

#endif
