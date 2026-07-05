#include <RcppArmadillo.h>
#include <chrono>
#include "saevws.h"
#include "armspp"

const double SEC_PER_MICROSEC = 1e-6;

// [[Rcpp::export]]
Rcpp::List gibbs_joint_cpp(const arma::vec& y, const arma::vec& s2,
	const arma::mat& X, const arma::mat& Z, const arma::vec& df,
	const Rcpp::List& init, const Rcpp::List& control, const Rcpp::List& fixed)
{
	unsigned long mem_start = rss_kb();

	unsigned int m = y.n_elem;
	unsigned int d1 = X.n_cols;
	unsigned int d2 = Z.n_cols;

	stopifnot(s2.n_elem == m, "s2.length() == m");
	stopifnot(df.n_elem == m, "df.length() == m");
	stopifnot(X.n_rows == m, "X.nrow() == m");
	stopifnot(Z.n_rows == m, "Z.nrow() == m");
	stopifnot(arma::all(df > 1), "all(df > 1)");

	const arma::mat& XtX = crossprod(X);
	const arma::mat& ZtZ = crossprod(Z);
	const arma::mat& XtX_inv = arma::inv(XtX);
	const arma::mat& ZtZ_inv = arma::inv(ZtZ);

	stopifnot(control.inherits("control_joint"), "control inherits from control_joint");
	unsigned int R = control["R"];
	unsigned int burn = control["burn"];
	unsigned int thin = control["thin"];
	unsigned int report = control["report"];
	const arma::uvec& save_latent = control["save_latent"];
	const Rcpp::List& inner_ctrl = control["inner"];
	bool record_mem = control["record_mem"];

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
	arma::mat beta_hist(R_keep, d1);
	arma::mat gamma_hist(R_keep, d2);
	arma::vec phi2_hist(R_keep);
	arma::vec tau2_hist(R_keep);
	arma::mat sigma2_hist(R_keep, save_latent.size());
	arma::mat theta_hist(R_keep, save_latent.size());

	arma::vec mem_hist(R);
	arma::uvec sigma2_rejects_hist(R);
	arma::uvec sigma2_comps_hist(R);
	arma::uvec sigma2_tunes_hist(R);
	arma::uvec sigma2_tuned_hist(R);
	arma::uvec sigma2_rejects_areas(m);
	mem_hist.fill(arma::datum::nan);
	sigma2_tunes_hist.fill(0);
	sigma2_tuned_hist.fill(0);
	sigma2_rejects_areas.fill(0);
	double avg_sigma2_comps = 0;

	// This is used if vws_method == "vws-tune" or vws_method == "vws-basic"
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = 1e6;
	args.tol_suff = tol_suff;
	args.tol_merge = tol_merge;

	// Set up initial values
	stopifnot(init.inherits("init_joint"), "init inherits from init_joint");
	arma::vec beta = init["beta"];
	arma::vec gamma = init["gamma"];
	arma::vec sigma2 = init["sigma2"];
	arma::vec theta = init["theta"];
	double tau2 = init["tau2"];
	double phi2 = init["phi2"];

	arma::vec Xbeta = X * beta;
	arma::vec Zgamma = Z * gamma;

	const arma::vec& kappa = (df - 1) / 2.0;
	arma::vec lambda = arma::pow(y - theta, 2) / 2.0 + df % s2 / 2.0;

	// This block is only used if vws_method == "vws-tune"
	std::vector<joint_sae_proposal> proposals;
		for (unsigned int i = 0; i < m; i++) {
		joint_sae_proposal x(Zgamma(i), std::sqrt(tau2), kappa(i), lambda(i));
		proposals.push_back(x);
	}

	// This block is only used if inner_method == "arms"
	std::mt19937_64 rng(static_cast<uint_fast64_t>(double(UINT_FAST64_MAX) * R::unif_rand()));
	arma::mat arms_quantiles(m, 3);
	arms_quantiles.col(0).fill(0.1);
	arms_quantiles.col(1).fill(1);
	arms_quantiles.col(2).fill(5);

	// This block is only used if inner_method == "amh"
	arma::vec zeta_mean(m);
	arma::vec zeta_g(m);
	arma::vec zeta_var_prop(m);
	zeta_mean.fill(0);
	zeta_g.fill(0);
	zeta_var_prop.fill(amh_varprop_init);

	// Set up fixed parameters
	stopifnot(fixed.inherits("fixed_joint"), "fixed inherits from fixed_joint");

	// Note: flat prior is assumed for each parameter in this model

	// Set up timers
	double elapsed_beta = 0;
	double elapsed_gamma = 0;
	double elapsed_phi2 = 0;
	double elapsed_tau2 = 0;
	double elapsed_sigma2 = 0;
	double elapsed_theta = 0;

	for (unsigned int rep = 0; rep < R; rep++)
	{
		// Draw [theta | rest]
		if (!fixed["theta"]) {
			auto st = std::chrono::system_clock::now();
			const arma::vec& ww = phi2 / (phi2 + sigma2);
			const arma::vec& mm = ww % y + (1 - ww) % Xbeta;
			const arma::vec& ss = arma::sqrt(ww % sigma2);
			theta = rnorm(mm, ss);
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_theta += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [beta | rest]
		if (!fixed["beta"]) {
			auto st = std::chrono::system_clock::now();
			const arma::vec& mm = arma::solve(XtX, crossprod(X, theta));
			const arma::mat& Omega = (1 / phi2) * XtX;
			beta = r_mvnorm_prec(mm, Omega);
			Xbeta = X * beta;
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_beta += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [gamma | rest]
		if (!fixed["gamma"]) {
			auto st = std::chrono::system_clock::now();
			const arma::vec& mm = arma::solve(ZtZ, crossprod(Z, log(sigma2)));
			const arma::mat& Omega = (1 / tau2) * ZtZ;
			gamma = r_mvnorm_prec(mm, Omega);
			Zgamma = Z * gamma;
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_gamma += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [phi2 | rest]
		if (!fixed["phi2"]) {
			auto st = std::chrono::system_clock::now();
			double aa = m/2.0 - 1;
			double bb = 1/2.0 * dot(theta - Xbeta);
			phi2 = r_invgamma(aa, bb);
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_phi2 += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [tau2 | rest]
		if (!fixed["tau2"]) {
			auto st = std::chrono::system_clock::now();
			double aa = m/2.0 - 1;
			double bb = 1/2.0 * dot(log(sigma2) - Zgamma);
			tau2 = r_invgamma(aa, bb);
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_tau2 += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [sigma2 | rest]
		if (!fixed["sigma2"]) {
			auto st = std::chrono::system_clock::now();
			lambda = arma::pow(y - theta, 2) / 2.0 + df % s2 / 2.0;

			if (strcmp(inner_method.get_cstring(), "imh") == 0) {
				/*
				* Independent Metropolis sampling step from You (2021)
				*/

				for (unsigned int i = 0; i < m; i++) {
					double sigma2_prop = r_invgamma(kappa(i), lambda(i));
					double u = R::runif(0, 1);
					double log_num = R::dlnorm(sigma2_prop, Zgamma(i), std::sqrt(tau2), true);
					double log_den = R::dlnorm(sigma2(i), Zgamma(i), std::sqrt(tau2), true);
					double log_ratio = std::min(log_num - log_den, 0.0);
					if (std::log(u) < log_ratio) {
						sigma2(i) = sigma2_prop;
					} else {
						sigma2_rejects_hist(rep)++;
						sigma2_rejects_areas(i)++;
					}
				}
			} else if (strcmp(inner_method.get_cstring(), "amh") == 0) {
				/*
				* Adaptive Metropolis-Hastings (AMH) sampling from Haario,
				* Saksman, & Tamminen (2005).
				*
				* Adjust the numerator and denominator with the Jacobian of the
				* transformation, since we exponentiate the candidate (which is
				* real-valued).
				*/

				for (unsigned int i = 0; i < m; i++) {
					double sigma2_prev = sigma2(i);
					double u = R::runif(0, 1);
					double phi = std::log(sigma2_prev);

					double phi_prop = R::rnorm(phi, std::sqrt(zeta_var_prop(i)));

					double sigma2_prop = std::exp(phi_prop);
					double log_num = R::dlnorm(sigma2_prop, Zgamma(i), std::sqrt(tau2), true) +
						d_invgamma(sigma2_prop, kappa(i), lambda(i), true) +
						phi_prop;
					double log_den = R::dlnorm(sigma2_prev, Zgamma(i), std::sqrt(tau2), true) +
						d_invgamma(sigma2_prev, kappa(i), lambda(i), true) +
						phi;
					double log_ratio = std::min(log_num - log_den, 0.0);
					if (std::log(u) < log_ratio) {
						sigma2(i) = sigma2_prop;
					} else {
						sigma2_rejects_hist(rep)++;
						sigma2_rejects_areas(i)++;
					}

					// Adapt the proposal distribution
					double t = rep;
					double zeta_cur = std::log(sigma2(i));
					double zeta_mean_prev = zeta_mean(i);
					zeta_mean(i) = (t * zeta_mean(i) + zeta_cur) / (t + 1);

					if (t == 0) {
						zeta_g(i) = std::pow(zeta_cur, 2) / (t + 1);
					} else  {
						zeta_g(i) = (t - 1) / t * zeta_g(i) +
							std::pow(zeta_mean_prev, 2) +
							std::pow(zeta_cur, 2) / t -
							(t + 1) / t * std::pow(zeta_mean(i), 2);
					}

					if (rep < 10) {
						zeta_var_prop(i) = amh_varprop_init;
					} else if (rep < tune) {
						zeta_var_prop(i) = std::pow(amh_varprop_c, 2)  * (zeta_g(i) + amh_varprop_eps);
					}
				}
			} else if (strcmp(inner_method.get_cstring(), "arms") == 0) {
				/*
				 * Adaptive Rejection Metropolis Sampling (ARMS)
				 * After sampling, grab a few quantiles from the proposal
				 * to use in the next round of the Gibbs sampler. This is the
				 * suggestion in Gilks et al (1992).
				*/
				for (unsigned int i = 0; i < m; i++) {
					const Rcpp::NumericVector& points = Rcpp::wrap(arms_quantiles.row(i));
					arms_joint_functor armsfun(Zgamma(i), std::sqrt(tau2), kappa(i), lambda(i));
					armspp::ARMS<double, arms_joint_functor, Rcpp::NumericVector::const_iterator>
					sigma2_dist(
						armsfun,  // log-density functor
						0,        // lower
						1000,     // upper
						0,        // convex adjustment
						points.begin(),
						points.size(),
						100,      // max_points
						true,     // use metropolis or not
						sigma2(i) // previous value
					);
					double sigma2_i = sigma2_dist(rng);
					bool is_reject = std::fabs(sigma2_i - sigma2(i)) < 1e-8;
					sigma2(i) = sigma2_i;
					sigma2_rejects_hist(rep) += is_reject;
					sigma2_rejects_areas(i) += is_reject;
					arms_quantiles(i,0) = sigma2_dist.envelopeQuantile(0.05);
					arms_quantiles(i,1) = sigma2_dist.envelopeQuantile(0.50);
					arms_quantiles(i,2) = sigma2_dist.envelopeQuantile(0.95);
				}
			} else if (strcmp(inner_method.get_cstring(), "vws-tune") == 0) {
				/*
				* Self-tuned VWS using vws package
				*
				* Do self-tuning if iteration is before `tune` argument.
				* Otherwise, update target parameter values and try to proceed
				* with partitions from past tuning.
				*
				* We have to call `update` either way if we want to ensure
				* the proposal will work as an envelope for rejection sampling.
				*/

				for (unsigned int i = 0; i < m; i++) {
					proposals[i].update(Zgamma(i), std::sqrt(tau2), kappa(i), lambda(i));

					if (rep < tune) {
						const auto& vws_out = vws::rejection_tune(proposals[i], 1, args);
						sigma2(i) = vws_out.draws[0];
						sigma2_rejects_hist(rep) += vws_out.rejects[0];
						sigma2_rejects_areas(i) += vws_out.rejects[0];
						sigma2_tunes_hist(rep) += vws_out.tunes[0];
						sigma2_tuned_hist(rep) += (vws_out.tunes[0] > 0);
					} else {
						const auto& vws_out = vws::rejection(proposals[i], 1, args);
						sigma2(i) = vws_out.draws[0];
						sigma2_rejects_hist(rep) += vws_out.rejects[0];
						sigma2_rejects_areas(i) += vws_out.rejects[0];
					}
				}
			} else if (strcmp(inner_method.get_cstring(), "vws-basic") == 0) {
				// VWS without tuning using vws package
				for (unsigned int i = 0; i < m; i++) {
					joint_sae_proposal h(Zgamma(i), std::sqrt(tau2), kappa(i), lambda(i));
					h.refine(N - 1, tol_suff);
					const auto& vws_out = vws::rejection(h, 1, args);
					sigma2(i) = vws_out.draws[0];
					sigma2_rejects_hist(rep) += vws_out.rejects[0];
					sigma2_rejects_areas(i) += vws_out.rejects[0];
				}
			} else {
				Rcpp::stop("Unrecognized method in inner_ctrl");
			}

			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_sigma2 += td.count() * SEC_PER_MICROSEC;
		}

		// Save total number of regions at this point
		sigma2_comps_hist(rep) = 0;
		for (unsigned int i = 0; i < m; i++) {
			sigma2_comps_hist(rep) += proposals[i].size();
		}
		avg_sigma2_comps = sigma2_comps_hist(rep) / double(m);

		// Save total memory usage in MB
		if (record_mem) {
			mem_hist(rep) = rss_kb() / 1024.0 - mem_start / 1024.0;
		}

		if (rep >= burn && rep % thin == 0) {
			beta_hist.row(rep_keep) = beta.t();
			gamma_hist.row(rep_keep) = gamma.t();
			phi2_hist[rep_keep] = phi2;
			tau2_hist[rep_keep] = tau2;

			for (unsigned int l = 0; l < save_latent.size(); l++) {
				unsigned int i = save_latent(l);
				sigma2_hist(rep_keep, l) = sigma2(i);
				theta_hist(rep_keep, l) = theta(i);
			}

			rep_keep++;
		}

		if ((rep + 1) % report == 0) {
			unsigned int s = (rep >= report) ? rep - report : 0;
			unsigned int rejects = arma::sum(sigma2_rejects_hist(arma::span(s, rep)));

			if (strcmp(inner_method.get_cstring(), "vws-tune") == 0)
			{
				unsigned int tunes = arma::sum(sigma2_tunes_hist(arma::span(s, rep)));
				logger("[%d] avg-N: %0.4f  tunes: %d  rejects: %d\n", rep + 1,
					avg_sigma2_comps, tunes, rejects);
			} else {
				logger("[%d] rejects: %d\n", rep + 1, rejects);
			}
		}

		Rcpp::checkUserInterrupt();
	}

	Rcpp::List elapsed = Rcpp::List::create(
		Rcpp::Named("beta") = elapsed_beta,
		Rcpp::Named("gamma") = elapsed_gamma,
		Rcpp::Named("phi2") = elapsed_phi2,
		Rcpp::Named("tau2") = elapsed_tau2,
		Rcpp::Named("sigma2") = elapsed_sigma2,
		Rcpp::Named("theta") = elapsed_theta
	);

	return Rcpp::List::create(
		Rcpp::Named("beta") = beta_hist,
		Rcpp::Named("gamma") = gamma_hist,
		Rcpp::Named("phi2") = phi2_hist,
		Rcpp::Named("tau2") = tau2_hist,
		Rcpp::Named("sigma2") = sigma2_hist,
		Rcpp::Named("theta") = theta_hist,
		Rcpp::Named("R_keep") = R_keep,
		Rcpp::Named("elapsed") = elapsed,
		Rcpp::Named("R") = R,
		Rcpp::Named("burn") = burn,
		Rcpp::Named("thin") = thin,
		Rcpp::Named("inner_method") = inner_method,
		Rcpp::Named("sigma2_rejects") = sigma2_rejects_hist,
		Rcpp::Named("sigma2_rejects_areas") = sigma2_rejects_areas,
		Rcpp::Named("sigma2_comps") = sigma2_comps_hist,
		Rcpp::Named("sigma2_tunes") = sigma2_tunes_hist,
		Rcpp::Named("sigma2_tuned") = sigma2_tuned_hist,
		Rcpp::Named("mem") = mem_hist,
		Rcpp::Named("m") = m,
		Rcpp::Named("inner_method") = inner_method
	);
}

