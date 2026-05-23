#include <RcppArmadillo.h>
#include <chrono>
#include "local-util.h"
#include "joint-sae-proposal.h"
#include "joint-vws-output.h"
#include "joint-vws-basic.h"
#include "joint-vws-tune.h"
#include "arms-joint-functor.h"
#include "armspp"

const double SEC_PER_MICROSEC = 1e-6;

// [[Rcpp::export]]
Rcpp::List gibbs_joint_cpp(const arma::vec& y, const arma::vec& s2,
	const arma::mat& X, const arma::mat& Z, const arma::vec& df,
	const Rcpp::List& init, const Rcpp::List& control, const Rcpp::List& fixed)
{
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
	const Rcpp::String& inner_method = inner_ctrl["method"];
	unsigned int max_rejects = inner_ctrl["max_rejects"];
	double tol_suff = inner_ctrl["tol_suff"];
	double tol_merge = inner_ctrl["tol_merge"];
	unsigned int N = inner_ctrl["N"];

	unsigned int rep_keep = 0;
	unsigned int R_keep = std::ceil((R - burn) / double(thin));

	// Set up histories
	arma::mat beta_hist(R_keep, d1);
	arma::mat gamma_hist(R_keep, d2);
	arma::vec phi2_hist(R_keep);
	arma::vec tau2_hist(R_keep);
	arma::mat sigma2_hist(R_keep, save_latent.size());
	arma::mat theta_hist(R_keep, save_latent.size());

	arma::uvec sigma2_rejections_hist(R);
	arma::uvec sigma2_knots_hist(R);
	arma::uvec sigma2_tunes_hist(R);
	arma::uvec sigma2_rejections_areas(m);
	sigma2_rejections_areas.fill(0);
	double avg_sigma2_knots = 0;

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

	arma::vec kappa = (df - 1) / 2.0;
	arma::vec lambda = arma::pow(y - theta, 2) / 2.0 + df % s2 / 2.0;
	arma::vec tau(m);
	tau.fill(std::sqrt(tau2));

	// This is only used if vws_method == "vws-tune"
	// std::vector<SAEProposal> proposals;
	// for (unsigned int i = 0; i < m; i++) {
	// 	SAEProposal x(Zgamma(i), tau(i), kappa(i), lambda(i));
	// 	proposals.push_back(x);
	// }

	// This is only used if inner_method == "vws-tune"
	std::vector<JointSAEMajorizer> proposals;
	for (unsigned int i = 0; i < m; i++) {
		JointSAEMajorizer maj;
		proposals.push_back(maj);
	}

	// This is only used if inner_method == "arms"
	std::mt19937_64 rng(static_cast<uint_fast64_t>(UINT_FAST64_MAX * R::unif_rand()));
	arma::mat arms_quantiles(m, 3);
	arms_quantiles.col(0).fill(0.1);
	arms_quantiles.col(1).fill(1);
	arms_quantiles.col(2).fill(5);

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

	for (unsigned int rep = 0; rep < R; rep++) {
		Rcpp::checkUserInterrupt();

		if ((rep + 1) % report == 0) {
			if (strcmp(inner_method.get_cstring(), "vws-tune") == 0)
			{
	        	logger("Starting rep %d with avg sigma2 knots %g\n",
	        		rep + 1, avg_sigma2_knots);
			} else {
	        	logger("Starting rep %d\n", rep + 1);
			}
		}

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
			tau.fill(std::sqrt(tau2));
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_tau2 += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [sigma2 | rest]
		if (!fixed["sigma2"]) {
			auto st = std::chrono::system_clock::now();
			kappa = (df - 1) / 2.0;
			lambda = arma::pow(y - theta, 2) / 2.0 + df % s2 / 2.0;

			if (strcmp(inner_method.get_cstring(), "imh") == 0) {
				// Independent Metropolis sampling step from You (2021)
				const arma::vec& u = arma::randu(m);
				arma::vec sigma2_prop(m);
				for (unsigned int i = 0; i < m; i++) {
					sigma2_prop(i) = r_invgamma(kappa(i), lambda(i));
				}
				const arma::vec& log_num = dlnorm(sigma2_prop, Zgamma, tau, true);
				const arma::vec& log_den = dlnorm(sigma2, Zgamma, tau, true);
				const arma::vec& log_ratio = arma::min(log_num - log_den, arma::zeros(m));
				const arma::uvec& idx = arma::find(arma::log(u) < log_ratio);
				sigma2(idx) = sigma2_prop.elem(idx);
				sigma2_rejections_hist(rep) = m - idx.n_elem;
				sigma2_rejections_areas += (arma::log(u) >= log_ratio);
				sigma2_tunes_hist(rep) = 0;
			} else if (strcmp(inner_method.get_cstring(), "arms") == 0) {
				/*
				 * After sampling, grab a few quantiles from the proposal
				 * to use in the next round of the Gibbs sampler. This is the
				 * suggestion in Gilks et al (1992).
				*/
				for (unsigned int i = 0; i < m; i++) {
					const Rcpp::NumericVector& points = Rcpp::wrap(arms_quantiles.row(i));
					ARMSJointFunctor armsfun(Zgamma(i), tau(i), kappa(i), lambda(i));
					armspp::ARMS<double, ARMSJointFunctor, Rcpp::NumericVector::const_iterator>
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
					sigma2_rejections_hist(rep) += is_reject;
					sigma2_rejections_areas(i) += is_reject;
					arms_quantiles(i,0) = sigma2_dist.envelopeQuantile(0.05);
					arms_quantiles(i,1) = sigma2_dist.envelopeQuantile(0.50);
					arms_quantiles(i,2) = sigma2_dist.envelopeQuantile(0.95);
				}
			} else if (strcmp(inner_method.get_cstring(), "vws-tune") == 0) {
				// Self-tuned VWS using vws package
				// for (unsigned int i = 0; i < m; i++) {
				// 	proposals[i].update(Zgamma(i), tau(i), kappa(i), lambda(i));
			    // 	const auto& vws_out = vws::rejection_tune(proposals[i], 1, args);
				// 	sigma2 = vws_out.draws[0];
				// 	sigma2_rejections_hist(rep) += vws_out.rejects[0];
				// 	sigma2_tunes_hist(rep) += vws_out.tunes[0];
				// }

				// Self-tuned VWS using customized implementation
				const JointVWSOutput& vws_out = joint_vws_tune(proposals,
				 	Zgamma, std::sqrt(tau2), kappa, lambda, max_rejects,
				 	tol_suff, tol_merge);
				sigma2 = vws_out.sigma2;
				sigma2_rejections_hist(rep) = arma::sum(vws_out.rejects);
				sigma2_tunes_hist(rep) = arma::sum(vws_out.updates);
			} else if (strcmp(inner_method.get_cstring(), "vws-basic") == 0) {
				// VWS without tuning using vws package
				// for (unsigned int i = 0; i < m; i++) {
				// 	proposals[i].update(Zgamma(i), tau(i), kappa(i), lambda(i));
			    // 	const auto& vws_out = vws::rejection(proposals[i], 1, args);
				// 	sigma2 = vws_out.draws[0];
				// 	sigma2_rejections_hist(rep) += vws_out.rejects[0];
				// 	sigma2_tunes_hist(rep) += vws_out.tunes[0];
				// }

				// VWS without tuning using customized implementation
				const JointVWSOutput& vws_out = joint_vws_basic(Zgamma,
					std::sqrt(tau2), kappa, lambda, N, tol_suff, max_rejects);
				sigma2 = vws_out.sigma2;
				sigma2_rejections_hist(rep) = arma::sum(vws_out.rejects);
				sigma2_tunes_hist(rep) = arma::sum(vws_out.updates);
			} else {
				Rcpp::stop("Unrecognized method in inner_ctrl");
			}

			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_sigma2 += td.count() * SEC_PER_MICROSEC;
		}

		// Save total number of knots at this point
		sigma2_knots_hist(rep) = 0;
		// for (unsigned int i = 0; i < m; i++) {
		// 	sigma2_knots_hist(rep) += proposals[i].size();
		// }
		for (unsigned int i = 0; i < proposals.size(); i++) {
			sigma2_knots_hist(rep) += proposals[i].get_knots().length();
		}
		avg_sigma2_knots = sigma2_knots_hist(rep) / double(m);

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
		Rcpp::Named("beta_hist") = beta_hist,
		Rcpp::Named("gamma_hist") = gamma_hist,
		Rcpp::Named("phi2_hist") = phi2_hist,
		Rcpp::Named("tau2_hist") = tau2_hist,
		Rcpp::Named("sigma2_hist") = sigma2_hist,
		Rcpp::Named("theta_hist") = theta_hist,
		Rcpp::Named("R_keep") = R_keep,
		Rcpp::Named("elapsed") = elapsed,
		Rcpp::Named("R") = R,
		Rcpp::Named("burn") = burn,
		Rcpp::Named("thin") = thin,
		Rcpp::Named("sigma2_rejections_hist") = sigma2_rejections_hist,
		Rcpp::Named("sigma2_rejections_areas") = sigma2_rejections_areas,
		Rcpp::Named("sigma2_knots_hist") = sigma2_knots_hist,
		Rcpp::Named("sigma2_tunes_hist") = sigma2_tunes_hist,
		Rcpp::Named("m") = m
	);
}

