#include "rcpp-optimize-hybrid.h"
#include "vws.h"
#include "fntl.h"

Rcpp::List optimize_hybrid_rcpp(const Rcpp::Function& f, double init,
	double lower, double upper, bool maximize, unsigned maxiter)
{
	const fntl::dfd& ff = [&](double x) -> double {
		const Rcpp::NumericVector& out = f(x);
		return out(0);
	};

	const auto& out = vws::optimize_hybrid(ff, init, lower, upper, maximize, maxiter);

	return Rcpp::List::create(
		Rcpp::Named("par") = out.par,
		Rcpp::Named("value") = out.value,
		Rcpp::Named("method") = out.method
	);
}

