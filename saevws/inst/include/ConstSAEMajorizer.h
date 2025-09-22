#ifndef CONST_SAE_MAJORIZER_H
#define CONST_SAE_MAJORIZER_H

#include <Rcpp.h>
#include "vws.h"
#include "local-util.h"

struct ConstSAEMajorizerOutput {
	Rcpp::NumericVector lower;
	Rcpp::NumericVector upper;
	Rcpp::NumericVector lxu;
	Rcpp::NumericVector lxl;
	Rcpp::NumericVector log_bound_regions;
	double log_bound;

	ConstSAEMajorizerOutput operator=(const ConstSAEMajorizerOutput& x) {
		this->lower = x.lower;
		this->upper = x.upper;
		this->lxu = x.lxu;
		this->lxl = x.lxl;
		this->log_bound_regions = x.log_bound_regions;
		this->log_bound = x.log_bound;
		return *this;
	}
};

class ConstSAEMajorizer {
public:
	ConstSAEMajorizer()
	: _knots(), _lower(), _upper()
	{
		_lower.push_back(0);
		_upper.push_back(R_PosInf);
	}

	ConstSAEMajorizer(const Rcpp::NumericVector& knots)
	: _knots(), _lower(), _upper()
	{
		_knots.insert(knots.begin(), knots.end());
		_lower.push_back(0);
		_lower.insert(_lower.end(), knots.begin(), knots.end());
		_upper.insert(_upper.end(), knots.begin(), knots.end());
		_upper.push_back(R_PosInf);
	}

	void add_knot(double x) {
		_knots.insert(x);

		_lower.clear();
		_lower.push_back(0);
		_lower.insert(_lower.end(), _knots.begin(), _knots.end());

		_upper.clear();
		_upper.insert(_upper.end(), _knots.begin(), _knots.end());
		_upper.push_back(R_PosInf);
	}

	void delete_knots(const std::vector<double>& x) {
		for (unsigned int i = 0; i < x.size(); i++) {
			_knots.erase(x[i]);
		}

		_lower.clear();
		_lower.push_back(0);
		_lower.insert(_lower.end(), _knots.begin(), _knots.end());

		_upper.clear();
		_upper.insert(_upper.end(), _knots.begin(), _knots.end());
		_upper.push_back(R_PosInf);
	}

	void update() {

	}

	Rcpp::NumericVector get_knots() const {
		return Rcpp::NumericVector(_knots.begin(), _knots.end());
	}

	Rcpp::NumericVector get_lower() const {
		return Rcpp::NumericVector(_lower.begin(), _lower.end());
	}

	Rcpp::NumericVector get_upper() const {
		return Rcpp::NumericVector(_upper.begin(), _upper.end());
	}

	// This version of the weight function takes the index `idx` of the region
	// to use in the calculation. This avoids a search when we already know
	// which one to use.
	double w_major(unsigned int idx, double kappa, double lambda, bool log) const
	{
		double mode = lambda / (kappa + 1);
		double lo = _lower[idx];
		double hi = _upper[idx];

		double x;
		if (mode <= lo) {
			x = lo;
		} else if (mode > hi) {
			x = hi;
		} else {
			x = mode;
		}
		return d_invgamma(x, kappa, lambda, log);
	}

	double w_minor(unsigned int idx, double kappa, double lambda, bool log) const
	{
		double mode = lambda / (kappa + 1);
		double lo = _lower[idx];
		double hi = _upper[idx];

		double hi_out = d_invgamma(hi, kappa, lambda, log);
		double lo_out = d_invgamma(lo, kappa, lambda, log);

		double out;
		if (mode <= lo) {
			out = hi_out;
		} else if (mode > hi) {
			out = lo_out;
		} else {
			out = std::min(lo_out, hi_out);
		}

		return log ? out : exp(out);
	}

	ConstSAEMajorizerOutput get_output(double mu, double tau, double kappa,
		double lambda) const
	{
		unsigned int N = _lower.size();
		const Rcpp::NumericVector& lower = get_lower();
		const Rcpp::NumericVector& upper = get_upper();

		Rcpp::NumericVector log_w_min(N);
		Rcpp::NumericVector log_w_max(N);

		for (unsigned int j = 0; j < N; j++) {
			log_w_max(j) = w_major(j, kappa, lambda, true);
			log_w_min(j) = w_minor(j, kappa, lambda, true);
		}

		ConstSAEMajorizerOutput out;
		out.lower = lower;
		out.upper = upper;

		const Rcpp::NumericVector& lp_region_lower = Rcpp::plnorm(lower, mu, tau, true, true);
		const Rcpp::NumericVector& lp_region_upper = Rcpp::plnorm(upper, mu, tau, true, true);
		const Rcpp::NumericVector& clp_region_lower = Rcpp::plnorm(lower, mu, tau, false, true);
		const Rcpp::NumericVector& clp_region_upper = Rcpp::plnorm(upper, mu, tau, false, true);
		const Rcpp::NumericVector& lp_region1 = vws::log_sub2_exp(lp_region_upper, lp_region_lower);
		const Rcpp::NumericVector& lp_region2 = vws::log_sub2_exp(clp_region_lower, clp_region_upper);
		const Rcpp::NumericVector& lp_region = Rcpp::pmax(lp_region1, lp_region2);

		out.lxu = log_w_max + lp_region;
		out.lxl = log_w_min + lp_region;
		out.log_bound_regions = vws::log_sub2_exp(out.lxu, out.lxl) - vws::log_sum_exp(out.lxu);
		out.log_bound = vws::log_sum_exp(out.log_bound_regions);

		return out;
	}

	ConstSAEMajorizer operator=(const ConstSAEMajorizer& x) {
		this->_knots = x._knots;
		this->_lower = x._lower;
		this->_upper = x._upper;
		return *this;
	}

private:
	std::set<double> _knots;
	std::vector<double> _lower;
	std::vector<double> _upper;
};

#endif
