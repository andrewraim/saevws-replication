#include <RcppArmadillo.h>
#include <chrono>
#include "saevws.h"
#include "armspp"

const double SEC_PER_MICROSEC = 1e-6;

// [[Rcpp::export]]
Rcpp::List gibbs_unmatch_cpp(const arma::vec& y, const arma::vec& sigma,
	const arma::mat& X, const Rcpp::List& init, const Rcpp::List& control,
	const Rcpp::List& fixed)
{
	unsigned int m = y.n_elem;
	unsigned int d = X.n_cols;

	stopifnot(X.n_rows == m, "X.n_rows == m");
	stopifnot(sigma.n_elem == m, "sigma2.n_elems == m");

	const arma::mat& XtX = crossprod(X);

	stopifnot(control.inherits("control_unmatch"), "control inherits from control_unmatch");
	unsigned int R = control["R"];
	unsigned int burn = control["burn"];
	unsigned int thin = control["thin"];
	unsigned int report = control["report"];
	const arma::uvec& save_latent = control["save_latent"];
	const Rcpp::List& inner_ctrl = control["inner"];
	const Rcpp::String& inner_method = inner_ctrl["method"];
	unsigned int max_rejects = inner_ctrl["max_rejects"];
	double tol_suff = inner_ctrl["tol_suff"];
	double tol_merge = inner_ctrl["tol_merge"];
	unsigned int N = inner_ctrl["N"];
	unsigned int tune = inner_ctrl["tune"];
	double amh_varprop_init = inner_ctrl["amh_varprop_init"];
	double amh_varprop_c = inner_ctrl["amh_varprop_c"];
	double amh_varprop_eps = inner_ctrl["amh_varprop_eps"];

	unsigned int rep_keep = 0;
	unsigned int R_keep = std::ceil((R - burn) / double(thin));

	// Set up histories
	arma::mat beta_hist(R_keep, d);
	arma::vec tau2_hist(R_keep);
	arma::mat mu_hist(R_keep, save_latent.size());

	arma::uvec mu_rejects_hist(R);
	arma::uvec mu_comps_hist(R);
	arma::uvec mu_tunes_hist(R);
	arma::uvec mu_tuned_hist(R);
	arma::uvec mu_rejects_areas(m);
	arma::vec mu_mem_hist(R);
	mu_tunes_hist.fill(0);
	mu_tuned_hist.fill(0);
	mu_rejects_areas.fill(0);
	double avg_mu_comps = 0;

	// This is used if vws_method == "vws-tune" or vws_method == "vws-basic"
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = 1e6;
	args.tol_suff = tol_suff;
	args.tol_merge = tol_merge;

	// Set up initial values
	stopifnot(init.inherits("init_unmatch"), "init inherits from init_unmatch");
	arma::vec beta = init["beta"];
	arma::vec mu = init["mu"];
	double tau2 = init["tau2"];

	arma::vec Xbeta = X * beta;

	// This is only used if vws_method == "vws-tune"
	std::vector<unmatch_sae_proposal> proposals;
	for (unsigned int i = 0; i < m; i++) {
	 	unmatch_sae_proposal x(y(i), sigma(i), Xbeta(i), std::sqrt(tau2));
	 	proposals.push_back(x);
	}

	// This is only used if inner_method == "arms"
	std::mt19937_64 rng(static_cast<uint_fast64_t>(UINT_FAST64_MAX * R::unif_rand()));
	arma::mat arms_quantiles(m, 3);
	arms_quantiles.col(0).fill(0.1);
	arms_quantiles.col(1).fill(1);
	arms_quantiles.col(2).fill(5);

	// This block is only used if inner_method == "amh"
	arma::vec mu_mean(m);
	arma::vec mu_g(m);
	arma::vec mu_varprop(m);
	mu_mean.fill(0);
	mu_g.fill(0);
	mu_varprop.fill(amh_varprop_init);

	// Set up fixed parameters
	stopifnot(fixed.inherits("fixed_unmatch"), "fixed inherits from fixed_unmatch");

	// Note: flat prior is assumed for each parameter in this model

	// Set up timers
	double elapsed_beta = 0;
	double elapsed_tau2 = 0;
	double elapsed_mu = 0;

	for (unsigned int rep = 0; rep < R; rep++)
	{
		// Draw [mu | rest]
		if (!fixed["mu"]) {
			auto st = std::chrono::system_clock::now();

			if (strcmp(inner_method.get_cstring(), "imh") == 0) {
				/*
				* Independent Metropolis sampling step from You & Rao (2002)
				*/

				for (unsigned int i = 0; i < m; i++) {
					double mu_prop = R::rlnorm(Xbeta(i), std::sqrt(tau2));

					double u = R::runif(0, 1);
					double log_num = R::dnorm(mu_prop, y(i), sigma(i), true);
					double log_den = R::dnorm(mu(i), y(i), sigma(i), true);
					const double log_ratio = std::min(log_num - log_den, 0.0);
					if (std::log(u) < log_ratio) {
						mu(i) = mu_prop;
					} else {
						mu_rejects_hist(rep)++;
						mu_rejects_areas(i)++;
					}
				}
			} else if (strcmp(inner_method.get_cstring(), "amh") == 0) {
				/*
				* Adaptive Metropolis-Hastings (AMH) sampling from Haario,
				* Saksman, & Tamminen (2005)
				*/

				for (unsigned int i = 0; i < m; i++) {
					// Draw a candidate and decide whether to accept it
					double mu_prev = mu(i);
					double mu_prop = R::rnorm(mu_prev, std::sqrt(mu_varprop(i)));

					double u = R::runif(0, 1);
					double log_num = R::dlnorm(mu_prop, Xbeta(i), std::sqrt(tau2), true) +
						R::dnorm(mu_prop, y(i), sigma(i), true);
					double log_den = R::dlnorm(mu_prev, Xbeta(i), std::sqrt(tau2), true) +
						R::dnorm(mu_prev, y(i), sigma(i), true);
					double log_ratio = std::min(log_num - log_den, 0.0);
					if (std::log(u) < log_ratio) {
						mu(i) = mu_prop;
					} else {
						mu_rejects_hist(rep)++;
						mu_rejects_areas(i)++;
					}

					// Adapt the proposal distribution
					double t = rep;
					double mu_mean_prev = mu_mean(i);
					mu_mean(i) = (t * mu_mean(i) + mu(i)) / (t + 1);

					if (t == 0) {
						mu_g(i) = std::pow(mu(i), 2) / (t + 1);
					} else  {
						mu_g(i) = (t - 1) / t * mu_g(i) +
							std::pow(mu_mean_prev, 2) +
							std::pow(mu(i), 2) / t -
							(t + 1) / t * std::pow(mu_mean(i), 2);
					}

					if (rep < 10) {
						mu_varprop(i) = amh_varprop_init;
					} else if (rep < tune) {
						mu_varprop(i) = std::pow(amh_varprop_c, 2)  * (mu_g(i) + amh_varprop_eps);
					}
				}
			} else if (strcmp(inner_method.get_cstring(), "arms") == 0) {
				/*
				 * After sampling, grab a few quantiles from the proposal
				 * to use in the next round of the Gibbs sampler. This is the
				 * suggestion in Gilks et al (1992).
				*/
				for (unsigned int i = 0; i < m; i++) {
					const Rcpp::NumericVector& points = Rcpp::wrap(arms_quantiles.row(i));
					arms_unmatch_functor armsfun(y(i), sigma(i), Xbeta(i), std::sqrt(tau2));
					armspp::ARMS<double, arms_unmatch_functor, Rcpp::NumericVector::const_iterator>
					mu_dist(
						armsfun,  // log-density functor
						0,        // lower
						1000,     // upper
						0,        // convex adjustment
						points.begin(),
						points.size(),
						100,      // max_points
						true,     // use metropolis or not
						mu(i)     // previous value
					);
					double mu_i = mu_dist(rng);
					bool is_reject = std::fabs(mu_i - mu(i)) < 1e-8;
					mu(i) = mu_i;
					mu_rejects_hist(rep) += is_reject;
					mu_rejects_areas(i) += is_reject;
					arms_quantiles(i,0) = mu_dist.envelopeQuantile(0.05);
					arms_quantiles(i,1) = mu_dist.envelopeQuantile(0.50);
					arms_quantiles(i,2) = mu_dist.envelopeQuantile(0.95);
				}
			} else if (strcmp(inner_method.get_cstring(), "vws-tune") == 0) {
				/*
				* Self-tuned VWS using vws package
				*
				* Do self-tuning if iteration is before `tune` argument.
				* Otherwise, update target parameter values and try to proceed
				* with partitions from past tuning.
				*/
				for (unsigned int i = 0; i < m; i++) {
					proposals[i].update(Xbeta(i), std::sqrt(tau2));

					if (rep < tune) {
						const auto& vws_out = vws::rejection_tune(proposals[i], 1, args);
						mu(i) = vws_out.draws[0];
						mu_rejects_areas(i) += vws_out.rejects[0];
						mu_rejects_hist(rep) += vws_out.rejects[0];
						mu_tunes_hist(rep) += vws_out.tunes[0];
						mu_tuned_hist(rep) += (vws_out.tunes[0] > 0);
					} else {
						const auto& vws_out = vws::rejection(proposals[i], 1, args);
						mu(i) = vws_out.draws[0];
						mu_rejects_areas(i) += vws_out.rejects[0];
						mu_rejects_hist(rep) += vws_out.rejects[0];
					}
				}
			} else if (strcmp(inner_method.get_cstring(), "vws-basic") == 0) {
				// VWS without tuning using vws package
				for (unsigned int i = 0; i < m; i++) {
					unmatch_sae_proposal h(y(i), sigma(i), Xbeta(i), std::sqrt(tau2));
			    	h.refine(N - 1, tol_suff);
			    	const auto& vws_out = vws::rejection(h, 1, args);
					mu(i) = vws_out.draws[0];
					mu_rejects_areas(i) += vws_out.rejects[0];
					mu_rejects_hist(rep) += vws_out.rejects[0];
				}
			} else {
				Rcpp::stop("Unrecognized method in inner_ctrl");
			}

			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_mu += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [beta | rest]
		if (!fixed["beta"]) {
			auto st = std::chrono::system_clock::now();
			const arma::vec& mm = arma::solve(XtX, crossprod(X, arma::log(mu)));
			const arma::mat& Omega = (1 / tau2) * XtX;
			beta = r_mvnorm_prec(mm, Omega);
			Xbeta = X * beta;
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_beta += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [tau2 | rest]
		if (!fixed["tau2"]) {
			auto st = std::chrono::system_clock::now();
			double aa = m / 2.0;
			double bb = 1 / 2.0 * dot(log(mu) - Xbeta);
			tau2 = r_invgamma(aa, bb);
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_tau2 += td.count() * SEC_PER_MICROSEC;
		}

		// Save total number of mixture components at this point
		mu_comps_hist(rep) = 0;
		for (unsigned int i = 0; i < m; i++) {
		 	mu_comps_hist(rep) += proposals[i].size();
		}
		avg_mu_comps = mu_comps_hist(rep) / double(m);

		// Save total memory usage for VWS proposals
		mu_mem_hist(rep) = 0;
		for (unsigned int i = 0; i < m; i++) {
		 	mu_mem_hist(rep) += mem(proposals[i]);
		}

		if (rep >= burn && rep % thin == 0) {
			beta_hist.row(rep_keep) = beta.t();
			tau2_hist[rep_keep] = tau2;

			for (unsigned int l = 0; l < save_latent.size(); l++) {
				unsigned int i = save_latent(l);
				mu_hist(rep_keep, l) = mu(i);
			}

			rep_keep++;
		}

		if ((rep + 1) % report == 0) {
			unsigned int s = (rep >= report) ? rep - report : 0;
			unsigned int rejects = arma::sum(mu_rejects_hist(arma::span(s, rep)));

			if (strcmp(inner_method.get_cstring(), "vws-tune") == 0)
			{
				unsigned int tunes = arma::sum(mu_tunes_hist(arma::span(s, rep)));
				logger("[%d] avg-N: %0.4f  tunes: %d  rejects: %d\n", rep + 1,
					avg_mu_comps, tunes, rejects);
			} else {
				logger("[%d] rejects: %d\n", rep + 1, rejects);
			}
		}

		Rcpp::checkUserInterrupt();

	}

	Rcpp::List elapsed = Rcpp::List::create(
		Rcpp::Named("beta") = elapsed_beta,
		Rcpp::Named("tau2") = elapsed_tau2,
		Rcpp::Named("mu") = elapsed_mu
	);

	return Rcpp::List::create(
		Rcpp::Named("beta_hist") = beta_hist,
		Rcpp::Named("tau2_hist") = tau2_hist,
		Rcpp::Named("mu_hist") = mu_hist,
		Rcpp::Named("R_keep") = R_keep,
		Rcpp::Named("elapsed") = elapsed,
		Rcpp::Named("R") = R,
		Rcpp::Named("burn") = burn,
		Rcpp::Named("thin") = thin,
		Rcpp::Named("inner_method") = inner_method,
		Rcpp::Named("mu_rejects_hist") = mu_rejects_hist,
		Rcpp::Named("mu_rejects_areas") = mu_rejects_areas,
		Rcpp::Named("mu_comps_hist") = mu_comps_hist,
		Rcpp::Named("mu_tunes_hist") = mu_tunes_hist,
		Rcpp::Named("mu_tuned_hist") = mu_tuned_hist,
		Rcpp::Named("mu_mem_hist") = mu_mem_hist,
		Rcpp::Named("m") = m
	);
}

