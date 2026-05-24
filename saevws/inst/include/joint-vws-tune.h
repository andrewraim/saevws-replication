#ifndef SAEVWS_JOINT_VWS_TUNE_H
#define SAEVWS_JOINT_VWS_TUNE_H

#include <RcppArmadillo.h>
#include "joint-sae-majorizer.h"
#include "joint-vws-output.h"
#include "local-util.h"

/*
* Self-tuned VWS step. Before removing knots, check that the overall rejection
* rate does not increase above `tol_suff`. This check involves some
* computational overhead.
*
* This function has at least one argument that cannot be invoked via Rcpp.
*/
inline joint_vws_output joint_vws_tune(std::vector<joint_sae_majorizer>& proposals,
	const arma::vec& mu, double tau, const arma::vec& kappa,
	const arma::vec& lambda, unsigned int max_rejects, double tol_suff,
	double tol_merge)
{
	unsigned int m = mu.n_elem;

	joint_vws_output out;
	out.sigma2 = arma::vec(m);
	out.rejects = arma::zeros<arma::uvec>(m);
	out.updates = arma::zeros<arma::uvec>(m);
	out.log_bound = arma::vec(m);

	for (unsigned int i = 0; i < m; i++)
	{
		joint_sae_majorizer& maj = proposals[i];
		joint_sae_majorizer_output maj_out = maj.get_output(mu(i), tau, kappa(i), lambda(i));

		if (i % 100 == 0) {
			Rcpp::checkUserInterrupt();
		}

		bool accept = false;
		while (!accept && out.rejects(i) < max_rejects)
		{
			const vws::dfdb& w =
			[&](double x, bool log = true) {
				return d_invgamma(x, kappa(i), lambda(i), log);
			};

			double v = R::runif(0, 1);

			/*
			* Take one draw from the VWS proposal:
			* 1. Draw finite mixture component idx from Categ(exp(lxu))
			* 2. Draw from truncated lognormal with index idx.
			*/
			unsigned int idx = vws::r_categ(1, maj_out.lxu, true, false)(0);
			double x = r_lnorm_trunc(mu(i), tau, maj_out.lower(idx),
				maj_out.upper(idx));
			double log_fx = w(x, true);
			double log_hx = maj.w_major(idx, kappa(i), lambda(i), true);
			double log_ratio = log_fx - log_hx;

			if (log(v) < log_ratio) {
				// Accept x as a draw from f(x)
				out.sigma2(i) = x;
				accept = true;
			} else {
				/* Reject x and update proposal according to thresholds */
				out.rejects(i)++;

				if (maj_out.log_bound < log(tol_suff)) {
					/*
					 * Drop regions that contribute very little. But only if
					 * overall bound is small enough, and dropping the region
					 * does not put us back over tol_suff threshold.
					*/
					for (unsigned int j = 0; j < maj_out.upper.length() - 1; j++) {
						if (maj_out.log_bound_regions(j) >= log(tol_merge)) {
							continue;
						}

						// Make a copy and drop the knot in the copy.
						joint_sae_majorizer maj0 = maj;
						maj0.delete_knots({ maj_out.upper(j) });

						// If the bound of the copy has not increased beyond
						// the threshold, replace the original with the copy.
						// Also make sure to save the updated output object.
						joint_sae_majorizer_output maj0_out = maj0.get_output(mu(i), tau, kappa(i), lambda(i));
						if (maj0_out.log_bound < log(tol_suff)) {
							maj = maj0;
							maj_out = maj0_out;
							out.updates(i)++;
						}
					}
				} else {
					/*
					 * Add the rejected draw as a knot and save the updated
					 * output object.
					*/
					maj.add_knot(x);
					maj_out = maj.get_output(mu(i), tau, kappa(i), lambda(i));
					out.updates(i)++;
				}
			}
		}

		if (out.rejects(i) >= max_rejects) {
			Rcpp::stop("Too many rejects in vws step");
		}

		out.log_bound(i) = maj_out.log_bound;
	}

	return out;
}

#endif
