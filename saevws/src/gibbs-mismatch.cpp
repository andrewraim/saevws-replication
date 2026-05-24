#include <RcppArmadillo.h>
#include <chrono>
#include "saevws.h"
#include "armspp"

const double SEC_PER_MICROSEC = 1e-6;

// [[Rcpp::export]]
Rcpp::List gibbs_mismatch_cpp(const arma::vec& y, const arma::vec& sigma,
	const arma::mat& X, const Rcpp::List& init, const Rcpp::List& control,
	const Rcpp::List& fixed)
{
	unsigned int m = y.n_elem;
	unsigned int d = X.n_cols;

	stopifnot(X.n_rows == m, "X.n_rows == m");
	stopifnot(sigma.n_elem == m, "sigma2.n_elems == m");

	const arma::mat& XtX = crossprod(X);

	stopifnot(control.inherits("control_mismatch"), "control inherits from control_mismatch");
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

	unsigned int rep_keep = 0;
	unsigned int R_keep = std::ceil((R - burn) / double(thin));

	// Set up histories
	arma::mat beta_hist(R_keep, d);
	arma::vec tau2_hist(R_keep);
	arma::mat mu_hist(R_keep, save_latent.size());

	arma::uvec mu_rejections_hist(R);
	arma::uvec mu_knots_hist(R);
	arma::uvec mu_tunes_hist(R);
	arma::uvec mu_rejections_areas(m);
	mu_rejections_areas.fill(0);
	double avg_mu_knots = 0;

	// This is used if vws_method == "vws-tune" or vws_method == "vws-basic"
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = 1e6;
	args.tol_suff = tol_suff;
	args.tol_merge = tol_merge;

	// Set up initial values
	stopifnot(init.inherits("init_mismatch"), "init inherits from init_mismatch");
	arma::vec beta = init["beta"];
	arma::vec mu = init["mu"];
	double tau2 = init["tau2"];

	arma::vec Xbeta = X * beta;

	// This is only used if vws_method == "vws-tune"
	std::vector<mismatch_sae_proposal> proposals;
	for (unsigned int i = 0; i < m; i++) {
	 	mismatch_sae_proposal x(y(i), sigma(i), Xbeta(i), std::sqrt(tau2));
	 	proposals.push_back(x);
	}

	// This is only used if inner_method == "arms"
	std::mt19937_64 rng(static_cast<uint_fast64_t>(UINT_FAST64_MAX * R::unif_rand()));
	arma::mat arms_quantiles(m, 3);
	arms_quantiles.col(0).fill(0.1);
	arms_quantiles.col(1).fill(1);
	arms_quantiles.col(2).fill(5);

	// Set up fixed parameters
	stopifnot(fixed.inherits("fixed_mismatch"), "fixed inherits from fixed_mismatch");

	// Note: flat prior is assumed for each parameter in this model

	// Set up timers
	double elapsed_beta = 0;
	double elapsed_tau2 = 0;
	double elapsed_mu = 0;

	for (unsigned int rep = 0; rep < R; rep++) {
		Rcpp::checkUserInterrupt();

		if ((rep + 1) % report == 0) {
			if (strcmp(inner_method.get_cstring(), "vws-tune") == 0)
			{
	        	logger("Starting rep %d with avg mu knots %g\n",
	        		rep + 1, avg_mu_knots);
			} else {
	        	logger("Starting rep %d\n", rep + 1);
			}
		}

		// Draw [mu | rest]
		if (!fixed["mu"]) {
			auto st = std::chrono::system_clock::now();
			// kappa = (df - 1) / 2.0;
			// lambda = arma::pow(y - theta, 2) / 2.0 + df % s2 / 2.0;

			if (strcmp(inner_method.get_cstring(), "imh") == 0) {
				// Independent Metropolis sampling step from You & Rao (2002)
				const arma::vec& u = arma::randu(m);
				arma::vec mu_prop(m);
				arma::vec log_num(m);
				arma::vec log_den(m);
				for (unsigned int i = 0; i < m; i++) {
					mu_prop(i) = R::rlnorm(Xbeta(i), std::sqrt(tau2));
					log_num(i) = R::dnorm(mu_prop(i), y(i), sigma(i), true);
					log_den(i) = R::dnorm(mu(i), y(i), sigma(i), true);
				}
				const arma::vec& log_ratio = arma::min(log_num - log_den, arma::zeros(m));
				const arma::uvec& idx = arma::find(arma::log(u) < log_ratio);
				mu(idx) = mu_prop.elem(idx);
				mu_rejections_hist(rep) = m - idx.n_elem;
				mu_rejections_areas += (arma::log(u) >= log_ratio);
				mu_tunes_hist(rep) = 0;
			} else if (strcmp(inner_method.get_cstring(), "arms") == 0) {
				/*
				 * After sampling, grab a few quantiles from the proposal
				 * to use in the next round of the Gibbs sampler. This is the
				 * suggestion in Gilks et al (1992).
				*/
				for (unsigned int i = 0; i < m; i++) {
					const Rcpp::NumericVector& points = Rcpp::wrap(arms_quantiles.row(i));
					arms_mismatch_functor armsfun(y(i), sigma(i), Xbeta(i), std::sqrt(tau2));
					armspp::ARMS<double, arms_mismatch_functor, Rcpp::NumericVector::const_iterator>
					mu_dist(
						armsfun,  // log-density functor
						0,        // lower
						1000,     // upper
						0,        // convex adjustment
						points.begin(),
						points.size(),
						100,      // max_points
						true,     // use metropolis or not
						mu(i) // previous value
					);
					double mu_i = mu_dist(rng);
					bool is_reject = std::fabs(mu_i - mu(i)) < 1e-8;
					mu(i) = mu_i;
					mu_rejections_hist(rep) += is_reject;
					mu_rejections_areas(i) += is_reject;
					arms_quantiles(i,0) = mu_dist.envelopeQuantile(0.05);
					arms_quantiles(i,1) = mu_dist.envelopeQuantile(0.50);
					arms_quantiles(i,2) = mu_dist.envelopeQuantile(0.95);
				}
			} else if (strcmp(inner_method.get_cstring(), "vws-tune") == 0) {
				// Self-tuned VWS using vws package
				for (unsigned int i = 0; i < m; i++) {
				 	proposals[i].update(Xbeta(i), std::sqrt(tau2));
			     	const auto& vws_out = vws::rejection_tune(proposals[i], 1, args);
				 	mu(i) = vws_out.draws[0];
				 	mu_rejections_areas(i) += vws_out.rejects[0];
				 	mu_rejections_hist(rep) += vws_out.rejects[0];
				 	mu_tunes_hist(rep) += vws_out.tunes[0];
				}
			} else if (strcmp(inner_method.get_cstring(), "vws-basic") == 0) {
				// VWS without tuning using vws package
				for (unsigned int i = 0; i < m; i++) {
					mismatch_sae_proposal h(y(i), sigma(i), Xbeta(i), std::sqrt(tau2));
			    	h.refine(N - 1, tol_suff);
			    	const auto& vws_out = vws::rejection(h, 1, args);
					mu(i) = vws_out.draws[0];
					mu_rejections_areas(i) += vws_out.rejects[0];
					mu_rejections_hist(rep) += vws_out.rejects[0];
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

		// Save total number of knots at this point
		mu_knots_hist(rep) = 0;
		for (unsigned int i = 0; i < m; i++) {
		 	mu_knots_hist(rep) += proposals[i].size();
		}
		avg_mu_knots = mu_knots_hist(rep) / double(m);

		if (rep >= burn && rep % thin == 0) {
			beta_hist.row(rep_keep) = beta.t();
			tau2_hist[rep_keep] = tau2;

			for (unsigned int l = 0; l < save_latent.size(); l++) {
				unsigned int i = save_latent(l);
				mu_hist(rep_keep, l) = mu(i);
			}

			rep_keep++;
		}
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
		Rcpp::Named("mu_rejections_hist") = mu_rejections_hist,
		Rcpp::Named("mu_rejections_areas") = mu_rejections_areas,
		Rcpp::Named("mu_knots_hist") = mu_knots_hist,
		Rcpp::Named("mu_tunes_hist") = mu_tunes_hist,
		Rcpp::Named("m") = m
	);
}

