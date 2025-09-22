#ifndef RCPP_OPTIMIZE_HYBRID_H
#define RCPP_OPTIMIZE_HYBRID_H

#include <Rcpp.h>

//' Hybrid Univariate Optimization
//'
//' Use Brent's method if a bounded search interval is specified. Otherwise use
//' BFGS method.
//'
//' @param f Objective function. Should take a scalar as an argument.
//' @param init Initial value for optimization variable.
//' @param lower Lower bound for search; may be \eqn{-\infty}.
//' @param upper Upper bound for search; may be \eqn{+\infty}.
//' @param maximize logical; if `TRUE`, optimization will be a maximization.
//' Otherwise, it will be a minimization.
//' @param maxiter Maximum number of iterations.
//'
//' @returns
//' \item{par}{Value of optimization variable.}
//' \item{value}{Value of optimization function.}
//' \item{method}{Description of result.}
//' \item{status}{Status code from BFGS or `0` otherwise.}
//'
//' @examples
//' f = function(x) { x^2 }
//' optimize_hybrid(f, init = 0, lower = -1, upper = 2, maximize = FALSE)
//' optimize_hybrid(f, init = 0, lower = -1, upper = Inf, maximize = FALSE)
//' optimize_hybrid(f, init = 0, lower = -Inf, upper = 1, maximize = FALSE)
//' optimize_hybrid(f, init = 0, lower = 0, upper = Inf, maximize = FALSE)
//' optimize_hybrid(f, init = 0, lower = -Inf, upper = 0, maximize = FALSE)
//'
//' f = function(x) { 1 - x^2 }
//' optimize_hybrid(f, init = 0, lower = -1, upper = 1, maximize = TRUE)
//' optimize_hybrid(f, init = 0, lower = -1, upper = 0, maximize = TRUE)
//' optimize_hybrid(f, init = 0, lower = 0, upper = 1, maximize = TRUE)
//'
//' @export
// [[Rcpp::export(name = "optimize_hybrid")]]
Rcpp::List optimize_hybrid_rcpp(const Rcpp::Function& f, double init,
	double lower, double upper, bool maximize = false,
	unsigned maxiter = 10000);

#endif

