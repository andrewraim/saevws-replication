#ifndef VWS_OPTIMIZE_HYBRID_H
#define VWS_OPTIMIZE_HYBRID_H

#include <Rcpp.h>
#include "fntl.h"
#include "logit.h"
#include "result.h"
#include "rect.h"

namespace vws {

/*
* A hybrid optimization for univariate problems. Use Brent if both bounds are
* finite; otherwise, use BFGS. In the latter case, if one of the bounds is
* finite, it is enforced via a transformation.
*
* - `f`: objective function.
* - `init`: initial value used with BFGS.
* - `lower`: lower bound.
* - `upper`: upper bound.
* - `maximize`: if `true`, optimization will be a maximization. Otherwise it
*   is a minimization.
* - `maxiter`: maximum number of iterations.
*
* Returns a `optimize_hybrid_result` structure.
*/
inline optimize_hybrid_result optimize_hybrid(const fntl::dfd& f, double init,
	double lower, double upper, bool maximize, unsigned maxiter = 100000)
{
	if (lower > upper) {
		Rcpp::stop("lower > upper");
	}

	double f_lower = f(lower);
	double f_upper = f(upper);
	bool f_lower_pos_inf = (f_lower > 0) && std::isinf(f_lower);
	bool f_lower_neg_inf = (f_lower < 0) && std::isinf(f_lower);
	bool f_upper_pos_inf = (f_upper > 0) && std::isinf(f_upper);
	bool f_upper_neg_inf = (f_upper < 0) && std::isinf(f_upper);

	optimize_hybrid_result out;
	out.status = 0;

	if (maximize && f_lower_pos_inf) {
		out.par = lower;
		out.value = R_PosInf;
		out.method = "Lower Limit Inf";
		return out;
	}

	if (maximize && f_upper_pos_inf) {
		out.par = upper;
		out.value = R_PosInf;
		out.method = "Upper Limit Inf";
		return out;
	}

	if (!maximize && f_lower_neg_inf) {
		out.par = lower;
		out.value = R_NegInf;
		out.method = "Lower Limit NegInf";
		return out;
	}

	if (!maximize && f_upper_neg_inf) {
		out.par = upper;
		out.value = R_NegInf;
		out.method = "Upper Limit NegInf";
		return out;
	}

	if (std::isfinite(lower) && std::isfinite(upper)) {
		// If both endpoints are finite, use Brent's method
	    fntl::optimize_args args;
	    args.fnscale = maximize ? -1 : 1;
	    args.maxiter = maxiter;

		const auto& opt_out = fntl::optimize_brent(f, lower, upper, args);
		out.par = opt_out.par;
		out.value = opt_out.value;
		out.method = "Brent";
	} else {
		// If one or both endpoints are infinite, use BFGS method
		fntl::bfgs_args args;
		args.fnscale = maximize ? -1 : 1;
		args.maxit = maxiter;

		// Transform to the interval (a,b]
		const fntl::dfd& tx = [&](double x) {
			return inv_rect(x, lower, upper);
		};

		// Compose function with transformation
	    const fntl::dfv& ff = [&](const Rcpp::NumericVector& x) {
			return f(tx(x(0)));
		};

		const auto& vinit = Rcpp::NumericVector::create(init);
		const fntl::bfgs_result& opt_out = fntl::bfgs(vinit, ff, args);

		out.par = tx(opt_out.par[0]);
		out.value = opt_out.value;
		out.method = "BFGS";
		out.status = fntl::to_underlying(opt_out.status);
	}

	// In case the function is strictly increasing or decreasing, check the
	// objective value at the endpoints.
	if (maximize) {
		if (f_lower > out.value) {
			out.par = lower;
			out.value = f_lower;
			out.method = "Max at Lower Limit";
		}
		if (f_upper > out.value) {
			out.par = upper;
			out.value = f_upper;
			out.method = "Max at Upper Limit";
		}
	} else {
		if (f_lower < out.value) {
			out.par = lower;
			out.value = f_lower;
			out.method = "Min at Lower Limit";
		}
		if (f_upper < out.value) {
			out.par = upper;
			out.value = f_upper;
			out.method = "Min at Upper Limit";
		}
	}

	return out;
}

}

#endif
