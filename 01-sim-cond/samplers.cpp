// [[Rcpp::depends(saevws, vws, fntl, RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>
#include "local-util.h"
#include "vws-step-tune-v1.h"

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
		double z = vws::r_invgamma(kappa, lambda);
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
Rcpp::List r_target(unsigned int n, double mu, double tau, double kappa,
	double lambda, double tol1, double tol2, unsigned int max_rejects)
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
		const VWSStepOutput& vws_out = vws_step_tune_v1(proposals, mu_vec,
			tau, kappa_vec, lambda_vec, max_rejects, tol1, tol2);
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
		Rcpp::Named("knots") = knots,
		Rcpp::Named("elapsed") = elapsed
	);
}

