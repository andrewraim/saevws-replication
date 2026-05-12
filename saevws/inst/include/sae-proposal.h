#ifndef SAE_PROPOSAL_H
#define SAE_PROPOSAL_H

// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

/*
* Define a subclass of FMMProposal that we can expose to R via Modules
*/
class SAEProposal : public vws::FMMProposal<double, vws::RealConstRegion>
{
public:
	SAEProposal(
		double mu,
		double tau,
		double kappa,
		double lambda)
	: vws::FMMProposal<double, vws::RealConstRegion>(supp(mu, tau, kappa, lambda))
	{
	}

	void update(double mu, double tau, double kappa, double lambda);

private:
	vws::RealConstRegion supp(double mu, double tau, double kappa, double lambda);
};

/*
* Implementation of member functions is below
*/

inline void SAEProposal::update(double mu, double tau, double kappa, double lambda)
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

	vws::UnivariateHelper helper(df, pf, qf);

	// std::set<vws::RealConstRegion>::iterator itr = _regions.begin();
	// for (; itr != _regions.end(); ++itr) {
	//	itr->set(w);
	//	itr->set(helper);
	// }

	for (unsigned int i = 0; i < _regions_vec.size(); i++) {
		_regions_vec[i].set(w);
		_regions_vec[i].set(helper);
	}
}

inline vws::RealConstRegion SAEProposal::supp(
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
		double out = d_invgamma(x, kappa, lambda, log);
		// Rprintf("maxopt: x=%g lo=%g hi=%g mode=%g out=%g\n",
		// 	x, lo, hi, mode, out);
		return out;
	};

	const vws::optimizer& minopt =
	[=](const vws::dfdb& w, double lo, double hi, bool log)
	{
		double hi_out = d_invgamma(hi, kappa, lambda, true);
		double lo_out = d_invgamma(lo, kappa, lambda, true);

		double out = (mode <= lo) ? hi_out :
		           (mode > hi) ? lo_out :
		           std::min(lo_out, hi_out);

		// Rprintf("minopt: lo=%g hi=%g mode=%g lo_out=%g hi_out=%g\n",
		// 	lo, hi, mode, lo_out, hi_out);
		return log ? out : exp(out);
	};

	vws::UnivariateHelper helper(df, pf, qf);
	vws::RealConstRegion out(0, R_PosInf, w, helper, maxopt, minopt);
	return out;
}

#endif
