#include <RcppArmadillo.h>
#include <chrono>
#include "local-util.h"
#include "VWSStepOutput.h"
#include "vws-step-basic.h"
#include "vws-step-tune.h"

const double SEC_PER_MICROSEC = 1e-6;

// [[Rcpp::export]]
Rcpp::List gibbs_cpp(const arma::vec& y, const arma::vec& s2,
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

	stopifnot(control.inherits("gibbs_control"), "control inherits from gibbs_control");
	unsigned int R = control["R"];
	unsigned int burn = control["burn"];
	unsigned int thin = control["thin"];
	unsigned int report = control["report"];
	const arma::uvec& save_latent = control["save_latent"];
	const Rcpp::List& vws_ctrl = control["vws"];
	const Rcpp::String& vws_method = vws_ctrl["method"];
	unsigned int max_rejects = vws_ctrl["max_rejects"];
	double tol1 = vws_ctrl["tol1"];
	double tol2 = vws_ctrl["tol2"];
	unsigned int N = vws_ctrl["N"];

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
	arma::uvec sigma2_knot_updates_hist(R);
	double avg_sigma2_knots = 0;

	// This is only used if vws_method == "vws-tune"
	std::vector<ConstSAEMajorizer> proposals;
	for (unsigned int i = 0; i < m; i++) {
		ConstSAEMajorizer maj;
		proposals.push_back(maj);
	}

	// Set up initial values
	stopifnot(init.inherits("gibbs_init"), "init inherits from gibbs_init");
	arma::vec beta = init["beta"];
	arma::vec gamma = init["gamma"];
	arma::vec sigma2 = init["sigma2"];
	arma::vec theta = init["theta"];
	double tau2 = init["tau2"];
	double phi2 = init["phi2"];

	double tau = std::sqrt(tau2);
	arma::vec Xbeta = X * beta;
	arma::vec Zgamma = Z * gamma;

	// Set up fixed parameters
	stopifnot(fixed.inherits("gibbs_fixed"), "fixed inherits from gibbs_fixed");

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
			if (strcmp(vws_method.get_cstring(), "vws-tune") == 0)
			{
	        	vws::logger("Starting rep %d with avg sigma2 knots %g\n",
	        		rep + 1, avg_sigma2_knots);
			} else {
	        	vws::logger("Starting rep %d\n", rep + 1);
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
			tau = std::sqrt(tau2);
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_tau2 += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [sigma2 | rest]
		if (!fixed["sigma2"]) {
			auto st = std::chrono::system_clock::now();
			const arma::vec& kappa = (df - 1) / 2.0;
			const arma::vec& lambda = arma::pow(y - theta, 2) / 2.0 + df % s2 / 2.0;
			arma::vec tau_vec(m);
			tau_vec.fill(tau);

			if (strcmp(vws_method.get_cstring(), "imh") == 0) {
				// Independent Metropolis sampling step from You (2021)
				const arma::vec& u = arma::randu(m);
				arma::vec sigma2_prop(m);
				for (unsigned int i = 0; i < m; i++) {
					sigma2_prop(i) = r_invgamma(kappa(i), lambda(i));
				}
				const arma::vec& log_num = dlnorm(sigma2_prop, Zgamma, tau_vec, true);
				const arma::vec& log_den = dlnorm(sigma2, Zgamma, tau_vec, true);
				const arma::vec& log_ratio = arma::min(log_num - log_den, arma::zeros(m));
				const arma::uvec& idx = arma::find(arma::log(u) < log_ratio);
				sigma2(idx) = sigma2_prop.elem(idx);
				sigma2_rejections_hist(rep) = m - idx.n_elem;
				sigma2_knot_updates_hist(rep) = 0;
			} else if (strcmp(vws_method.get_cstring(), "vws-tune") == 0) {
				// Self-tuned VWS
				const VWSStepOutput& vws_out = vws_step_tune(proposals,
					Zgamma, tau, kappa, lambda, max_rejects, tol1, tol2);
				sigma2 = vws_out.sigma2;
				sigma2_rejections_hist(rep) = arma::sum(vws_out.rejects);
				sigma2_knot_updates_hist(rep) = arma::sum(vws_out.updates);
			} else if (strcmp(vws_method.get_cstring(), "vws-basic") == 0) {
				// VWS without tuning
				const VWSStepOutput& vws_out = vws_step_basic(Zgamma, tau,
					kappa, lambda, N, tol1, max_rejects);
				sigma2 = vws_out.sigma2;
				sigma2_rejections_hist(rep) = arma::sum(vws_out.rejects);
				sigma2_knot_updates_hist(rep) = arma::sum(vws_out.updates);
			} else {
				Rcpp::stop("Unrecognized method in vws_ctrl");
			}

			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_sigma2 += td.count() * SEC_PER_MICROSEC;
		}

		// Save total number of knots at this point
		sigma2_knots_hist(rep) = 0;
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
		Rcpp::Named("sigma2_knots_hist") = sigma2_knots_hist,
		Rcpp::Named("sigma2_knot_updates_hist") = sigma2_knot_updates_hist,
		Rcpp::Named("m") = m
	);
}
