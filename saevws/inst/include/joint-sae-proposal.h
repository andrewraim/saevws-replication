#ifndef SAEVWS_JOINT_SAE_PROPOSAL_H
#define SAEVWS_JOINT_SAE_PROPOSAL_H

// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

/*
* Define a subclass of fmm_proposal that we can expose to R via Modules
*/
class joint_sae_proposal : public vws::fmm_proposal<double, vws::real_const_region>
{
public:
	joint_sae_proposal(
		double mu,
		double tau,
		double kappa,
		double lambda)
	: vws::fmm_proposal<double, vws::real_const_region>(supp(mu, tau, kappa, lambda))
	{
	}

	void update(double mu, double tau, double kappa, double lambda);

private:
	vws::real_const_region supp(double mu, double tau, double kappa, double lambda);
};

/*
* Implementation of member functions is below
*/

inline void joint_sae_proposal::update(double mu, double tau, double kappa, double lambda)
{
	const vws::dfdb& w = [=](double x, bool log = true) -> double {
		double out = (x > 0) ? d_invgamma(x, kappa, lambda, true) : R_NegInf;
		return log ? out : std::exp(out);
	};

	const fntl::density& df = [=](double x, bool log = false) {
		return R::dlnorm(x, mu, tau, log);
	};

	const fntl::cdf& pf = [=](double q, bool lower = true, bool log = false) {
		return R::plnorm(q, mu, tau, lower, log);
	};

	const fntl::quantile& qf = [=](double p, bool lower = true, bool log = false) {
		return R::qlnorm(p, mu, tau, lower, log);
	};

	vws::univariate_helper helper(df, pf, qf);

		double mode = lambda / (kappa + 1);

	// We have a simple closed-form max and min that we can use here.

	const vws::optimizer& maxopt =
	[=](const vws::dfdb& w, double lo, double hi, bool log)
	{
		double x = (mode <= lo) ? lo :
			(mode > hi) ? hi :
			mode;
		return w(x, log);
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

	// Update weight function in each region in the proposal
	std::set<vws::real_const_region>::iterator itr = _regions.begin();
	for (; itr != _regions.end(); ++itr) {
		vws::real_const_region& reg = const_cast<vws::real_const_region&>(*itr);
		reg.set_w(w);
		reg.set_helper(helper);
		reg.set_maxopt(maxopt);
		reg.set_minopt(minopt);
		reg.init();
	}
}

inline vws::real_const_region joint_sae_proposal::supp(
	double mu,
	double tau,
	double kappa,
	double lambda)
{
	const vws::dfdb& w = [=](double x, bool log = true) -> double {
		double out = (x > 0) ? d_invgamma(x, kappa, lambda, true) : R_NegInf;
		return log ? out : std::exp(out);
	};

	fntl::density df = [=](double x, bool log = false) {
		return R::dlnorm(x, mu, tau, log);
	};

	fntl::cdf pf = [=](double q, bool lower = true, bool log = false) {
		return R::plnorm(q, mu, tau, lower, log);
	};

	fntl::quantile qf = [=](double p, bool lower = true, bool log = false) {
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
		return w(x, log);
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

	vws::univariate_helper helper(df, pf, qf);
	vws::real_const_region out(0, R_PosInf, w, helper, maxopt, minopt);
	return out;
}

#endif
