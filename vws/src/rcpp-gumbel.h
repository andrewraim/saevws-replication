#ifndef RCPP_GUMBEL_H
#define RCPP_GUMBEL_H

#include <Rcpp.h>

//' Gumbel Distribution
//'
//' Functions for the Gumbel distribution with density
//' \deqn{
//' f(x \mid \mu, \sigma) =
//' \frac{1}{\sigma}
//' \exp\{ -\{ (x - \mu) / \sigma + e^{-(x - \mu) / \sigma} \} \}
//' }
//'
//' @param n Number of draws.
//' @param x Vector; argument of density.
//' @param p Vector; argument of cumulative distribution function.
//' @param q Vector; argument of quantile function.
//' @param mu Location parameter.
//' @param sigma Scale parameter.
//' @param lower Logical; if `TRUE` (default), probabilities are
//' \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
//' @param log Logical; if `TRUE`, probabilities p are given as \eqn{log(p)}
//'
//' @return
//' `d_gumbel` computes the density, `r_gumbel` generates random deviates,
//' `p_gumbel` computes the CDF, and `q_gumbel` computes quantiles.
//'
//' @examples
//' mu = 1
//' sigma = 2
//' x = r_gumbel(100000, mu, sigma)
//' xx = seq(min(x), max(x), length.out = 100)
//'
//' plot(density(x))
//' lines(xx, d_gumbel(xx, mu, sigma), lty = 2, col = "blue", lwd = 2)
//'
//' plot(ecdf(x))
//' lines(xx, p_gumbel(xx, mu, sigma), lty = 2, col = "blue", lwd = 2)
//'
//' pp = seq(0, 1, length.out = 102) |> head(-1) |> tail(-1)
//' qq = quantile(x, probs = pp)
//' plot(pp, qq)
//' lines(pp, q_gumbel(pp, mu, sigma), lty = 2, col = "blue", lwd = 2)
//'
//' @name Gumbel
//' @export
// [[Rcpp::export(name = "r_gumbel")]]
 Rcpp::NumericVector r_gumbel_rcpp(unsigned int n, double mu = 0,
 	double sigma = 1);

//' @name Gumbel
//' @export
// [[Rcpp::export(name = "d_gumbel")]]
Rcpp::NumericVector d_gumbel_rcpp(const Rcpp::NumericVector& x, double mu = 0,
	double sigma = 1, bool log = false);

//' @name Gumbel
//' @export
// [[Rcpp::export(name = "p_gumbel")]]
Rcpp::NumericVector p_gumbel_rcpp(const Rcpp::NumericVector& q, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false);

//' @name Gumbel
//' @export
// [[Rcpp::export(name = "q_gumbel")]]
Rcpp::NumericVector q_gumbel_rcpp(const Rcpp::NumericVector& p, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false);

#endif
