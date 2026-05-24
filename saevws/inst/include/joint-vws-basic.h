#ifndef SAEVWS_JOINT_VWS_BASIC_H
#define SAEVWS_JOINT_VWS_BASIC_H

#include <RcppArmadillo.h>
#include "joint-sae-majorizer.h"
#include "joint-vws-output.h"
#include "lognormal-helper.h"
#include "local-util.h"

/*
* VWS step without tuning.
*/
inline joint_vws_output joint_vws_basic(const arma::vec& mu, double tau,
	const arma::vec& kappa, const arma::vec& lambda, unsigned int N, double tol,
	unsigned int max_rejects)
{
	unsigned int m = mu.n_elem;

	joint_vws_output out;
	out.sigma2 = arma::vec(m);
	out.rejects = arma::zeros<arma::uvec>(m);
	out.log_bound = arma::vec(m);

	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.action = fntl::error_action::STOP;

	for (unsigned int i = 0; i < m; i++) {
		if (i % 100 == 0) {
			Rcpp::checkUserInterrupt();
		}

		const vws::dfdb& w =
		[&](double x, bool log = true) {
			return d_invgamma(x, kappa(i), lambda(i), log);
		};

		fntl::density df = [&](double x, bool log = false) {
			return R::dlnorm(x, mu(i), tau, log);
		};
		fntl::cdf pf = [&](double q, bool lower = true, bool log = false) {
			return R::plnorm(q, mu(i), tau, lower, log);
		};
		fntl::quantile qf = [&](double p, bool lower = true, bool log = false) {
			return R::qlnorm(p, mu(i), tau, lower, log);
		};
		vws::univariate_helper helper(df, pf, qf);

		// Restrict range to something smaller than (0, Inf] to avoid numerical
		// issues in matrix computation of likelihood.
		vws::real_const_region supp(0, R_PosInf, w, helper);
		vws::fmm_proposal<double, vws::real_const_region> h({ supp });

		const Rcpp::NumericVector& refine_out = h.refine(N - 1, tol);
		const vws::rejection_result<double>& vws_out = vws::rejection(h, 1, args);

		out.sigma2(i) = vws_out.draws[0];
		out.rejects(i) = vws_out.rejects[0];
		out.log_bound(i) = refine_out[refine_out.length() - 1];
	}

	return out;
}

#endif

