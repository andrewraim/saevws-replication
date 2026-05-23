#ifndef MISMATCH_SAE_PROPOSAL_H
#define MISMATCH_SAE_PROPOSAL_H

// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

/*
* Define a subclass of FMMProposal that we can expose to R via Modules
*/
class MismatchSAEProposal : public vws::FMMProposal<double, vws::RealConstRegion>
{
public:
	MismatchSAEProposal(
		double y,
		double sigma,
		double xbeta,
		double tau)
	: vws::FMMProposal<double, vws::RealConstRegion>(supp(y, sigma, xbeta, tau))
	{
	}

	void update(double xbeta, double tau);

private:
	vws::RealConstRegion supp(double y, double sigma, double xbeta, double tau);
};

/*
* Implementation of member functions is below
*/

inline void MismatchSAEProposal::update(double xbeta, double tau)
{
	const vws::dfdb& w = [=](double mu, bool log = true) -> double {
		double out = (mu > 0) ? R::dlnorm(mu, xbeta, tau, true) : R_NegInf;
		return log ? out : std::exp(out);
	};

	// We have a simple closed-form max and min that we can use here.
	double tau2 = tau * tau;
	double mode = exp(xbeta - tau2);

	const vws::optimizer& maxopt =
	[=](const vws::dfdb& w, double lo, double hi, bool log) {
		double y = (mode > hi) ? hi :
			(mode < lo) ? lo :
			mode;
		double out = w(y, true);
		return log ? out : exp(out);
	};

	const vws::optimizer& minopt =
	[=](const vws::dfdb& w, double lo, double hi, bool log) {
		double lwa = w(lo, true);
		double lwb = w(hi, true);
		double out = std::min(lwa, lwb);
		return log ? out : exp(out);
	};

	// Update weight function in each region in the proposal
	std::set<vws::RealConstRegion>::iterator itr = _regions.begin();
	for (; itr != _regions.end(); ++itr) {
		vws::RealConstRegion& reg = const_cast<vws::RealConstRegion&>(*itr);
		reg.set_w(w);
		reg.set_maxopt(maxopt);
		reg.set_minopt(minopt);
		reg.init();
	}
}

inline vws::RealConstRegion MismatchSAEProposal::supp(
	double y,
	double sigma,
	double xbeta,
	double tau)
{
	const vws::dfdb& w = [=](double mu, bool log = true) -> double {
		double out = (mu > 0) ? R::dlnorm(mu, xbeta, tau, true) : R_NegInf;
		return log ? out : std::exp(out);
	};

	const fntl::density& df = [=](double mu, bool log = false) {
		return R::dnorm(mu, y, sigma, log);
	};

	const fntl::cdf& pf = [=](double q, bool lower = true, bool log = false) {
		return R::pnorm(q, y, sigma, lower, log);
	};

	const fntl::quantile& qf = [=](double p, bool lower = true, bool log = false) {
		return R::qnorm(p, y, sigma, lower, log);
	};

	// We have a simple closed-form max and min that we can use here.
	double tau2 = tau * tau;
	double mode = exp(xbeta - tau2);

	const vws::optimizer& maxopt =
	[=](const vws::dfdb& w, double lo, double hi, bool log) {
		double y = (mode > hi) ? hi :
			(mode < lo) ? lo :
			mode;
		double out = w(y, true);
		return log ? out : exp(out);
	};

	const vws::optimizer& minopt =
	[=](const vws::dfdb& w, double lo, double hi, bool log) {
		double lwa = w(lo, true);
		double lwb = w(hi, true);
		double out = std::min(lwa, lwb);
		return log ? out : exp(out);
	};

	vws::UnivariateHelper helper(df, pf, qf);
	vws::RealConstRegion out(0, R_PosInf, w, helper, maxopt, minopt);
	return out;
}

#endif
