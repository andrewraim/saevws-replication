#ifndef LOCAL_UTIL_H
#define LOCAL_UTIL_H

#include <RcppArmadillo.h>
#include "vws.h"

inline void stopifnot(bool cond, const char* fmt, ...)
{
	if (cond) { return; }

	// Insert placeholders into formatted string; see
	// <https://stackoverflow.com/q/1056411>
	char msg[256];
	va_list args;
	va_start(args, fmt);
	vsnprintf(msg, 255, fmt, args);
	va_end(args);

	Rcpp::stop(std::string(msg) + " is not TRUE");
}

/*
inline void logger(const char* fmt, ...)
{
	const Rcpp::Datetime& dt = Rcpp::Datetime(time(NULL));

	// Insert placeholders into formatted string; see
	// <https://stackoverflow.com/q/1056411>
	char msg[256];
	va_list args;
	va_start(args, fmt);
	vsnprintf(msg, 255, fmt, args);
	va_end(args);

	Rprintf("%04d:%02d:%02d %02d:%02d:%02d - %s", dt.getYear(), dt.getMonth(),
		dt.getDay(), dt.getHours(), dt.getMinutes(), dt.getSeconds(), msg);
}
*/

inline arma::mat crossprod(const arma::mat& X)
{
	return arma::trans(X) * X;
}

inline arma::mat crossprod(const arma::mat& X, const arma::mat& Y)
{
	return arma::trans(X) * Y;
}

inline arma::vec rnorm(const arma::vec& mu, const arma::vec& sigma)
{
	stopifnot(mu.n_elem == sigma.n_elem, "mu.n_elem == sigma.n_elem");
	unsigned int n = mu.n_elem;
	return sigma % arma::randn(n) + mu;
}

inline double d_invgamma(double x, double shape, double scale, bool log)
{
	double out = R_NegInf;
	if (x > 0) {
		double lnc = std::lgamma(shape) - shape * std::log(scale);
		out = -(shape + 1)*std::log(x) - scale / x - lnc;
	}
	return log ? out : exp(out);
}

inline arma::vec dlnorm(const arma::vec& x, const arma::vec& mean, const arma::vec& sd, bool log)
{
	unsigned int n = x.n_elem;
	if (n != mean.n_elem){ Rcpp::stop("n != mean.n_elem"); }
	if (n != sd.n_elem){ Rcpp::stop("n != sd.n_elem"); }

	arma::vec out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = R::dlnorm(x(i), mean(i), sd(i), log);
	}

	return out;
}

inline double r_invgamma(double a, double b)
{
	return 1 / R::rgamma(a, 1 / b);
}

inline double dot(const arma::vec& x)
{
	return arma::dot(x, x);
}

inline arma::mat r_mvnorm_prec(const arma::vec& mu, const arma::mat& Omega)
{
	stopifnot(mu.n_elem == Omega.n_rows, "mu.n_elem == Omega.n_rows");
	stopifnot(Omega.n_rows == Omega.n_cols, "Omega.n_rows == Omega.n_cols");
	unsigned int k = mu.n_elem;
	const arma::vec& z = arma::randn(k);
	const arma::mat& A = arma::chol(Omega);
	return arma::solve(A, z) + mu;
}

inline double r_lnorm_trunc(double meanlog, double sdlog, double lo, double hi,
	bool lower = true, bool log = false)
{
	double p = R::runif(0, 1);
	double lpp = log ? p : std::log(p);
	lpp = lower ? lpp : vws::log_sub2_exp(0, lpp);

	double lpa = R::plnorm(lo, meanlog, sdlog, true, true);
	double lpb = R::plnorm(hi, meanlog, sdlog, true, true);
	double lp = vws::log_sub2_exp(lpb, lpa);

	double clpa = R::plnorm(lo, meanlog, sdlog, false, true);
	double clpb = R::plnorm(hi, meanlog, sdlog, false, true);
	double clp = vws::log_sub2_exp(clpa, clpb);

	double lpr = std::max(lp, clp);

	double lq;
	if (std::isinf(lpp) || std::isinf(lpr)) {
		lq = lpa;
	} else {
		lq = vws::log_add2_exp(lpa, lpp + lpr);
	}

	// Protect against log-probabilities greater than zero, which can happen
	// numerically (assuming there are no mistakes).
	lq = std::min(lq, 0.0);

	double out = R::qlnorm(lq, meanlog, sdlog, true, true);

	// Protect against quantiles outside of the support, which can happen
	// numerically (assuming there are no mistakes).
	return std::max(std::min(out, hi), lo);
}

inline Rcpp::NumericVector log_cumsum_exp(const Rcpp::NumericVector& x)
{
	unsigned int k = x.length();
	Rcpp::NumericVector out(k);
	double lsum = R_NegInf;

	for (unsigned int i = 0; i < k; i++) {
		lsum = vws::log_add2_exp(x(i), lsum);
		out(i) = lsum;
	}

	return out;
}

inline unsigned int q_categ(double p, const Rcpp::NumericVector& cprobs)
{
	unsigned int k = cprobs.length();
	unsigned int lo = 0;
	unsigned int hi = k-1;

	while (lo + 1 < hi) {
		unsigned int i = std::ceil((lo + hi) / 2.0);
		bool ind = (p <= cprobs(i));
		lo = ind*lo + (1-ind)*i;
		hi = ind*i + (1-ind)*hi;
	}

	return p <= cprobs(lo) ? lo : hi;
}

#endif
