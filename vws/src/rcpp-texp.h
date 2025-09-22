#ifndef RCPP_TEXP_H
#define RCPP_TEXP_H

#include <Rcpp.h>

//' Truncated Exponential
//'
//' Functions for the truncated exponential distribution with density
//' \deqn{
//' f(x \mid \kappa, a, b) =
//' \frac{
//' \kappa e^{\kappa x}
//' }{
//' e^{\kappa b} - e^{\kappa a}
//' } \textrm{I}(a \leq x \leq b)
//' }
//'
//' @param n Number of desired draws.
//' @param x Vector; argument of density.
//' @param p Vector; argument of cumulative distribution function.
//' @param q Vector; argument of quantile function.
//' @param s Vector; argument of moment-generating function.
//' @param kappa Parameter of distribution.
//' @param a Lower limit of integral.
//' @param b Upper limit of integral.
//' @param lo Lower limit of support.
//' @param hi Upper limit of support.
//' @param log logical; if `TRUE`, probabilities are interpreted on the
//' log-scale. Otherwise, probabilities are interpreted on the original scale.
//' @param lower logical; if `TRUE`, probabilities are \eqn{P(X \leq x)};
//' otherwise \eqn{P(X > x)}.
//'
//' @return
//' \item{n_texp}{Normalizing constant.}
//' \item{integrate_texp}{Compute integral of density over an interval.}
//' \item{d_texp}{Compute density.}
//' \item{p_texp}{Compute cumulative probabilities.}
//' \item{q_texp}{Compute quantiles.}
//' \item{r_texp}{Generate variates.}
//' \item{mgf_texp}{Compute moment generating function.}
//'
//' @examples
//' set.seed(1234)
//' kappa = 2
//' lo = -0.5
//' hi = +0.5
//'
//' x = r_texp(100000, kappa, lo, hi)
//' xx = seq(min(x), max(x), length.out = 100)
//'
//' hist(x, probability = TRUE)
//' lines(xx, d_texp(xx, kappa, lo, hi), lty = 2, col = "blue", lwd = 2)
//'
//' plot(ecdf(x))
//' lines(xx, p_texp(xx, kappa, lo, hi), lty = 2, col = "blue", lwd = 2)
//'
//' pp = seq(0, 1, length.out = 102) |> head(-1) |> tail(-1)
//' qq = quantile(x, probs = pp)
//' plot(pp, qq)
//' lines(pp, q_texp(pp, kappa, lo, hi), lty = 2, col = "blue", lwd = 2)
//'
//' n_texp(kappa, lo, hi)
//' n_texp(kappa, lo, hi, log = TRUE)
//'
//' integrate_texp(0, hi, kappa, lo, hi)
//' integrate_texp(0, hi, kappa, lo, hi, log = TRUE)
//'
//' mgf_texp(0, kappa, lo, hi)
//' mgf_texp(0, kappa, lo, hi, log = TRUE)
//' mgf_texp(0.5, kappa, lo, hi)
//' mgf_texp(0.5, kappa, lo, hi, log = TRUE)
//'
//' @name TExp
//' @export
// [[Rcpp::export(name = "n_texp")]]
double n_texp_rcpp(double kappa, double lo, double hi, bool log = false);

//' @name TExp
//' @export
// [[Rcpp::export(name = "integrate_texp")]]
double integrate_texp_rcpp(double a, double b, double kappa, double lo,
	double hi, bool log = false);

//' @name TExp
//' @export
// [[Rcpp::export(name = "d_texp")]]
Rcpp::NumericVector d_texp_rcpp(const Rcpp::NumericVector& x, double kappa,
	double lo, double hi, bool log = false);

//' @name TExp
//' @export
// [[Rcpp::export(name = "p_texp")]]
Rcpp::NumericVector p_texp_rcpp(const Rcpp::NumericVector& q, double kappa,
	double lo, double hi, bool lower = true, bool log = false);

//' @name TExp
//' @export
// [[Rcpp::export(name = "q_texp")]]
Rcpp::NumericVector q_texp_rcpp(const Rcpp::NumericVector& p, double kappa,
	double lo, double hi, bool lower = true, bool log = false);

//' @name TExp
//' @export
// [[Rcpp::export(name = "r_texp")]]
Rcpp::NumericVector r_texp_rcpp(unsigned int n, double kappa, double lo,
	double hi);

//' @name TExp
//' @export
// [[Rcpp::export(name = "mgf_texp")]]
Rcpp::NumericVector mgf_texp_rcpp(const Rcpp::NumericVector& s, double kappa,
	double lo, double hi, bool log = false);

#endif
