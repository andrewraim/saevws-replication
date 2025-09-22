#ifndef VWS_REJECTION_H
#define VWS_REJECTION_H

#include <Rcpp.h>
#include "logger.h"
#include "Region.h"
#include "FMMProposal.h"
#include "result.h"
#include "typedefs.h"

namespace vws {

/*
* Structure of optional arguments for rejection sampler.
*
*  - `max_rejects`: Maximum number of rejections to tolerate before bailing out.
*
*  - `report_period`: specifies the period in which progress should be reported
*    (printed to the screen as a log message).
*
*  - `max_rejects_action`: what should happen if `max_rejects` rejections have
*    been obtained during sampling. The default action `STOP` results in an
*    exception being thrown; here, any successful draws that may have been
*    obtained are not returned.
*
*  - `log_ratio_ub`: it is possible numerically for log-ratio
*    $\log[f_0(x) / h_0(x)]$ to be greater than zero. This condition should
*    not occur otherwise, and usually indicates a mistake in user code. This
*    argument is the maximum value allowed where an exception will not be
*    thrown.
*/
struct rejection_args
{
	unsigned int max_rejects = std::numeric_limits<unsigned int>::max();
	unsigned int report_period = std::numeric_limits<unsigned int>::max();
	error_action max_rejects_action = error_action::STOP;
	double log_ratio_ub = 1e-5;
};

/*
*  Accept-reject algorithm using a VWS proposal.
*
*  - `h`: a VWS proposal.
*  - `n`: number of desired draws.
*  - `args`: additional arguments
*
*  Returns a structure with saved draws and rejection counts.
*/
template <typename T, typename R>
inline rejection_result<T>
rejection(const FMMProposal<T,R>& h, unsigned int n, const rejection_args& args)
{
	std::vector<T> draws;
	std::vector<unsigned int> rejects(n, 0L);

	unsigned int N_rejects = 0;
	bool accept = false;

	unsigned int max_rejects = args.max_rejects;
	unsigned int report_period = args.report_period;
	error_action max_rejects_action = args.max_rejects_action;
	double log_ratio_ub = args.log_ratio_ub;

	// The constant M in the acceptance ratio is always M = 1.
	double log_M = 0;

	for (unsigned int i = 0; i < n; i++) {
		accept = false;
		while (!accept && N_rejects < max_rejects) {
			double v = ::R::runif(0, 1);
			const std::vector<T>& draws_new = h.r(1);
			const T& x = draws_new[0];
			double log_fx = h.d_target_unnorm(x);
			double log_hx = h.d(x, false, true);
			double log_ratio = log_fx - log_hx - log_M;

			if (log_ratio > log_ratio_ub) {
				Rcpp::stop("log_ratio = %g > %g  x = %g  log f(x) = %g  log h(x) = %g",
					log_ratio, log_ratio_ub, x, log_fx, log_hx);
			} else if (log(v) < log_ratio) {
				// Accept x as a draw from f(x)
				draws.push_back(x);
				accept = true;
			} else {
				// Reject x
				N_rejects++;
				rejects[i]++;
			}

			// Report progress after `report` candidates
			unsigned int N_accepts = i + accept;
			if ((N_rejects + N_accepts) % report_period == 0) {
				logger("After %d candidates, %d accepts and %d rejects\n",
					N_accepts + N_rejects, N_accepts, N_rejects);
			}
		}
	}

	if (N_rejects >= max_rejects) {
		switch(max_rejects_action) {
			case error_action::STOP:
				Rcpp::stop("Reached maximum number of rejects: %d\n", max_rejects);
				break;
			case error_action::WARNING:
				Rcpp::warning("Reached maximum number of rejects: %d\n", max_rejects);
				break;
			case error_action::MESSAGE:
				Rprintf("Reached maximum number of rejects: %d\n", max_rejects);
				break;
			default:
				break;
		}
	}

	rejection_result<T> out;
	out.draws = draws;
	out.rejects = rejects;
	return out;
}

/*
*  Accept-reject algorithm using VWS proposal.
*
*  - `h`: a VWS proposal.
*  - `n`: number of desired draws.
*
*  Returns a structure with saved draws and rejection counts.
*/
template <typename T, typename R>
inline rejection_result<T> rejection(const FMMProposal<T,R>& h, unsigned int n)
{
	rejection_args ctrl;
	return rejection(h, n, ctrl);
}

}

#endif
