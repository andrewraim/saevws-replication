#ifndef RCPP_INVGAMMA_H
#define RCPP_INVGAMMA_H

#include <Rcpp.h>

//' Inverse Gamma Distribution
//'
//' Functions for the Inverse Gamma distribution parameterized by shape and
//' rate, whose density is
//' \deqn{
//' f(x \mid a, b) =
//' \frac{
//' b^a
//' }{
//' \Gamma(a)
//' }
//' x^{a-1} e^{-b x} \textrm{I}(x \geq 0)
//' }
//'
//' @param n Number of draws.
//' @param x Vector; argument of density.
//' @param q Vector; argument of cumulative distribution function.
//' @param p Vector; argument of quantile function.
//' @param a Shape parameter.
//' @param b Rate parameter.
//' @param lower logical; if `TRUE` (default), probabilities are
//' \eqn{P(X \leq x)}; otherwise, \eqn{P(X > x)}.
//' @param log If `TRUE`, return densities and probabilities on the log-scale.
//'
//' @return
//' `d_invgamma` computes the density, `r_invgamma` generates random deviates,
//' `p_invgamma` computes the CDF, and `q_invgamma` computes quantiles.
//'
//' @examples
//' a = 10
//' b = 5
//' x = r_invgamma(100000, a, b)
//' xx = seq(0, 3, length.out = 100)
//'
//' plot(density(x))
//' lines(xx, d_invgamma(xx, a, b), lty = 2, col = "blue", lwd = 2)
//'
//' plot(ecdf(x))
//' lines(xx, p_invgamma(xx, a, b), lty = 2, col = "blue", lwd = 2)
//'
//' pp = seq(0, 1, length.out = 102) |> head(-1) |> tail(-1)
//' qq = quantile(x, probs = pp)
//' plot(pp, qq)
//' lines(pp, q_invgamma(pp, a, b), lty = 2, col = "blue", lwd = 2)
//'
//' @name InverseGamma
//' @export
// [[Rcpp::export(name = "r_invgamma")]]
Rcpp::NumericVector r_invgamma_rcpp(unsigned int n, double a, double b);

//' @name InverseGamma
//' @export
// [[Rcpp::export(name = "d_invgamma")]]
Rcpp::NumericVector d_invgamma_rcpp(const Rcpp::NumericVector& x, double a,
	double b, bool log = false);

//' @name InverseGamma
//' @export
// [[Rcpp::export(name = "p_invgamma")]]
Rcpp::NumericVector p_invgamma_rcpp(const Rcpp::NumericVector& q, double a,
	double b, bool lower = true, bool log = false);

//' @name InverseGamma
//' @export
// [[Rcpp::export(name = "q_invgamma")]]
Rcpp::NumericVector q_invgamma_rcpp(const Rcpp::NumericVector& p, double a,
	double b, bool lower = true, bool log = false);

#endif
