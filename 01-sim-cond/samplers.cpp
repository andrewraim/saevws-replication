// [[Rcpp::depends(saevws, vws, fntl, RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>
#include "saevws.h"

const double SEC_PER_MICROSEC = 1e-6;

// Metropolis sampling step from You (2021)
// [[Rcpp::export]]
Rcpp::List r_metro(unsigned int n, double init, double mu, double tau,
	double kappa, double lambda)
{
	Rcpp::NumericVector draws(n+1);
	double rejects = 0;
	double x = init;
	draws(0) = x;

	for (unsigned int i = 0; i < n; i++) {
		double u = R::runif(0, 1);
		double z = r_invgamma(kappa, lambda);
		double log_num = R::dlnorm(z, mu, tau, true);
		double log_den = R::dlnorm(x, mu, tau, true);
		double log_ratio = std::min(log_num - log_den, 0.0);
		bool accept = log(u) < log_ratio;
		x = accept ? z : x;
		draws(i+1) = x;
		rejects += !accept;
	}

	return Rcpp::List::create(
		Rcpp::Named("draws") = draws,
		Rcpp::Named("rejects") = rejects
	);
}

// [[Rcpp::export]]
Rcpp::List r_target_old(unsigned int n, double mu, double tau, double kappa,
	double lambda, double tol_suff, double tol_merge, unsigned int max_rejects)
{
	auto st = std::chrono::system_clock::now();

	std::vector<ConstSAEMajorizer> proposals = { ConstSAEMajorizer() };

	unsigned int one = 1L;
	arma::vec mu_vec(one);
	arma::vec lambda_vec(one);
	arma::vec kappa_vec(one);

	mu_vec.fill(mu);
	kappa_vec.fill(kappa);
	lambda_vec.fill(lambda);

	arma::vec draws(n);
	arma::vec log_bounds(n);
	arma::uvec rejections(n);
	arma::uvec knots(n);

	for (unsigned int i = 0; i < n; i++) {
		const VWSStepOutput& vws_out = vws_step_tune(proposals, mu_vec,
			tau, kappa_vec, lambda_vec, max_rejects, tol_suff, tol_merge);
		draws(i) = vws_out.sigma2[0];
		rejections(i) = vws_out.rejects[0];
		log_bounds(i) = vws_out.log_bound[0];
		knots(i) = proposals[0].get_knots().length();
	}

	auto et = std::chrono::system_clock::now();
	auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
	auto elapsed = td.count() * SEC_PER_MICROSEC;

	return Rcpp::List::create(
		Rcpp::Named("draws") = draws,
		Rcpp::Named("log_bounds") = log_bounds,
		Rcpp::Named("rejections") = rejections,
		Rcpp::Named("regions") = knots,
		Rcpp::Named("elapsed") = elapsed
	);
}


// [[Rcpp::export]]
Rcpp::List r_target(unsigned int n, double mu, double tau, double kappa,
	double lambda, double tol_suff, double tol_merge, unsigned int max_rejects,
	unsigned int report = 1e8)
{
	auto st = std::chrono::system_clock::now();

	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.report = report;
	args.action = fntl::error_action::STOP;
	args.tol_suff = tol_suff;
	args.tol_merge = tol_merge;

	const vws::dfdb& w = [&](double x, bool log = true) -> double {
		double out = (x > 0) ? d_invgamma(x, kappa, lambda, true) : R_NegInf;
		return log ? out : std::exp(out);
	};

	fntl::density df = [&](double x, bool log = false) {
		return R::dlnorm(x, mu, tau, log);
	};

	fntl::cdf pf = [&](double q, bool lower = true, bool log = false) {
		return R::plnorm(q, mu, tau, lower, log);
	};

	fntl::quantile qf = [&](double p, bool lower = true, bool log = false) {
		return R::qlnorm(p, mu, tau, lower, log);
	};

	double mode = lambda / (kappa + 1);

	// We have a simple closed-form max and min that we can use here.

	const vws::optimizer& maxopt =
	[=](const vws::dfdb& w, double lo, double hi, bool log)
	{
		double x = (mode <= lo) ? lo :
			(mode > hi) ? hi :
			mode;
		return w(x, true);
	};

	const vws::optimizer& minopt =
	[=](const vws::dfdb& w, double lo, double hi, bool log)
	{
		double hi_out = w(hi, true);
		double lo_out = w(lo, true);
		double out = (mode <= lo) ? hi_out :
			(mode > hi) ? lo_out :
			std::min(lo_out, hi_out);
		return log ? out : exp(out);
	};

	vws::UnivariateHelper helper(df, pf, qf);
	vws::RealConstRegion supp(0, R_PosInf, w, helper, maxopt, minopt);
	vws::FMMProposal<double, vws::RealConstRegion> h(supp);

	// Start with just one region and tune from there.
	const vws::rejection_result<double>& out = vws::rejection_tune(h, n, args);

	auto et = std::chrono::system_clock::now();
	auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
	auto elapsed = td.count() * SEC_PER_MICROSEC;

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("log_bounds") = out.log_bounds,
		Rcpp::Named("rejections") = out.rejects,
		Rcpp::Named("regions") = out.regions,
		Rcpp::Named("elapsed") = elapsed
	);
}
