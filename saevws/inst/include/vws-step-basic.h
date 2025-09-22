#ifndef VWS_STEP_BASIC_H
#define VWS_STEP_BASIC_H

#include <RcppArmadillo.h>
#include "ConstSAEMajorizer.h"
#include "VWSStepOutput.h"
#include "LognormalHelper.h"
#include "local-util.h"


/*
* VWS step without tuning.
*/
inline VWSStepOutput vws_step_basic(const arma::vec& mu, double tau,
	const arma::vec& kappa, const arma::vec& lambda, unsigned int N, double tol,
	unsigned int max_rejects)
{
	unsigned int m = mu.n_elem;

	VWSStepOutput out;
	out.sigma2 = arma::vec(m);
	out.rejects = arma::zeros<arma::uvec>(m);
	out.log_bound = arma::vec(m);

	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.max_rejects_action = vws::error_action::STOP;

	for (unsigned int i = 0; i < m; i++) {
		if (i % 100 == 0) {
			Rcpp::checkUserInterrupt();
		}

		const vws::uv_weight_function& w =
		[&](double x, bool log = true) {
			// return d_invgamma(x, kappa(i) - 1, lambda(i), log);
			return d_invgamma(x, kappa(i), lambda(i), log);
		};

		// Restrict range to something smaller than (0, Inf] to avoid numerical
		// issues in matrix computation of likelihood.
		LognormalHelper helper(mu(i), tau);
		vws::UnivariateConstRegion supp(0, R_PosInf, w, helper);
		vws::FMMProposal<double, vws::UnivariateConstRegion> h({ supp });

		const Rcpp::NumericVector& refine_out = h.refine(N - 1, tol);
		const vws::rejection_result<double>& vws_out = vws::rejection(h, 1, args);

		out.sigma2(i) = vws_out.draws[0];
		out.rejects(i) = vws_out.rejects[0];
		out.log_bound(i) = refine_out[refine_out.length() - 1];
	}

	return out;
}

#endif
