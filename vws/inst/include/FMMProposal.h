#ifndef VWS_FMM_PROPOSAL_H
#define VWS_FMM_PROPOSAL_H

#include <Rcpp.h>
#include "typedefs.h"
#include "log-sum-exp.h"
#include "categ.h"
#include "which.h"
#include "logger.h"
#include "Region.h"
#include <iterator>

namespace vws {

/*
* FMMProposal is a class which represents a VWS proposal: a finite mixture that
* implements certain operations. It takes two template arguments:
*
* - `T`: the data type of the support. For example, for a continuous
*   univariate support, we may use `double`.
* - `R`: a subclass of `Region` toused by the proposal. `Region` subclasses
*   implement any problem-specific logic used in the proposal.
*/
template <class T, class R>
class FMMProposal
{
public:
	/*
	* Constructor for FMMProposal based on a vector of one or more Regions.
	*/
	FMMProposal(const std::vector<R>& regions);

	/*
	* Access region characteristics
	*
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*
	* Return vectors which are ordered so that elements correspond to regions
	* $1, \ldots, N$.
	*
 	* - `get_xi_upper`: get $(\overline{\xi}_1, \ldots \overline{\xi}_1)$.
	* - `get_xi_lower`: get $(\underline{\xi}_1, \ldots \underline{\xi}_1)$.
	* - `get_bifurcatable`: vector of indicators of whether regions are
	*   bifurcatable.
	* - `get_pi`: mixing weights $(\pi_1, \ldots, \pi_N)$.
	* - `rejection_bound_regions`: vector of contributions to the rejection
	*   bound. These sum to the value of `rejection_bound`.
	*/
	Rcpp::NumericVector get_xi_upper(bool log = true) const;
	Rcpp::NumericVector get_xi_lower(bool log = true) const;
	Rcpp::LogicalVector get_bifurcatable() const;
	Rcpp::NumericVector get_pi(bool log = false) const;
	Rcpp::NumericVector rejection_bound_regions(bool log = false) const;

	/*
	* Access number of regions.
	*/
	unsigned int get_N() const {
		return _log_xi_upper->size();
	}

	/*
	* The following provide const iterators to internal data structures. These
	* may be more efficient (but less convenient) than the accessors because
	* copies are not made here.
	*/
	typename std::set<R>::const_iterator regions_begin() const {
		return _regions.begin();
	}
	typename std::set<R>::const_iterator regions_end() const {
		return _regions.end();
	}
	Rcpp::NumericVector::const_iterator log_xi_upper_begin() const {
		return _log_xi_upper->begin();
	}
	Rcpp::NumericVector::const_iterator log_xi_upper_end() const {
		return _log_xi_upper->end();
	}
	Rcpp::NumericVector::const_iterator log_xi_lower_begin() const {
		return _log_xi_lower->begin();
	}
	Rcpp::NumericVector::const_iterator log_xi_lower_end() const {
		return _log_xi_lower->end();
	}
	Rcpp::LogicalVector::const_iterator bifurcatable_begin() const {
		return _bifurcatable->begin();
	}
	Rcpp::LogicalVector::const_iterator bifurcatable_end() const {
		return _bifurcatable->end();
	}

	/*
	* Upper bound for rejection probability.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	double rejection_bound(bool log = false) const;

	/*
	* Normalizing constant for the proposal distribution.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	double nc(bool log = false) const;

	/*
	*  Generate draws from proposal distribution.
	*  - `n`: desired number of draws.
	*  Returns a vector of draws.
	*/
	std::vector<T> r(unsigned int n = 1) const;

	/*
	*  Generate draws from proposal distribution.
	*  - `n`: desired number of draws.
	*  Returns a pair. First element is vector of draws. Second element is
	*  vector with indices, each represents the region from which the draw was
	*  taken.
	*/
	std::pair<std::vector<T>, std::vector<unsigned int>> r_ext(unsigned int n = 1) const;

	/*
	* Compute density of proposal distribution.
	* - `x`: argument of density.
	* - `normalize`: whether to compute normalized or unnormalized density.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	double d(const T& x, bool normalize = true, bool log = false) const;

	/*
	* Compute $\overline{w}_(x)$.
	*
	* - `x`: a density argument.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	double w_major(const T& x, bool log = true) const;

	/*
	* Compute $f_0(x) = w(x) g(x)$.
	*
	* - `x`: a density argument.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	double d_target_unnorm(const T& x, bool log = true) const;

	/*
	* Prepare summary table of regions which compose the proposal distribution.
	*/
	Rcpp::DataFrame summary() const;

	/*
	* Print summary table of regions which compose the proposal distribution.
	* Limit the display to `n` regions.
	*/
	void print(unsigned int n = 5) const;

	/*
	* Refine the proposal.
	*
	* - `knots`: vector of $N$ knots at which to bifurcate.
	*
	* This function modifies the proposal object in place.
	*
	* It returns a vector of $N+1$ values that represent rejection bounds (on
	* the log-scale) corresponding to the bifurcations from `knots`: the
	* element with index `j` is the bound after applying the first $j$ `knots`.
	*/
	Rcpp::NumericVector refine(const std::vector<T>& knots);

	/*
	* Refine the proposal.
	*
	* - `N`: number of additional desired bifurcations.
	* - `tol`: tolerance for rejection bound.
	* - `greedy`: whether to use greedy selection to select regions for
	*   bifurcation.
	* - `report`: specifies the period in which progress should be reported
	*   (printed to the screen as a log message).
	*
	* This function modifies the proposal object in place.
	*
	* If `greedy` is `true`, the selected region to bifurcate will be the one
	* with the largest contribution to the rejection bound. Otherwise, the
	* selection is random, with probability proportional to regions'
	* contribution. Note that in either case, non-bifurcatable regions are
	* excluded from the candidates.
	*
	* If tolerance `tol` is achieved for the rejection bound before $N$
	* bifurcation steps are attempted, the operation is considered to be
	* complete.
	*
	* This function returns a vector of up to $N+1$ values that represent
	* rejection bounds (on the log-scale) that correspond to the bifurcation
	* steps: the element with index `j` corresponds to the bound at the $j$th
	* step.
	*/
	Rcpp::NumericVector refine(unsigned int N, double tol = 0,
		bool greedy = false, unsigned int report = uint_max);

private:

	/*
	* Recompute internal auxiliary data structures.
	*
	* The following data structures are based on the set `_regions`:
	* - `_regions_vec`
	* - `_log_xi_upper`
	* - `_log_xi_lower`
	* - `_bifurcatable`
	*
	* Here we construct them from the current state of `_regions`.
	*/
	void recache();

	std::set<R> _regions;
	std::vector<R> _regions_vec;
	std::unique_ptr<Rcpp::NumericVector> _log_xi_upper;
	std::unique_ptr<Rcpp::NumericVector> _log_xi_lower;
	std::unique_ptr<Rcpp::LogicalVector> _bifurcatable;
};

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::refine(const std::vector<T>& knots)
{
	unsigned int N = knots.size() + 1;

	std::vector<double> log_bdd_hist;

	log_bdd_hist.push_back(rejection_bound(true));

	for (unsigned int j = 0; j < knots.size(); j++) {
		unsigned int L = _regions.size();

		const T& x = knots[j];
		const R& reg0 = _regions.begin()->singleton(x);

		typename std::set<R>::const_iterator itr_lower = _regions.upper_bound(reg0);
		--itr_lower;

		if (itr_lower == _regions.end()) {
			Rcpp::warning("Could not find region which is upper bound for knots[%d]", j);
			continue;
		}

		// In this case, x was smaller (according to operator<) than all of the
		// regions in the set. So print a warning and continue.
		if (!itr_lower->s(x)) {
			Rcpp::warning("Could not find a region that contains knots[%d]", j);
			continue;
		}

		if (!itr_lower->is_bifurcatable()) {
			Rcpp::warning("Region that contains knots[%d] is not bifurcatable", j);
			continue;
		}

		const R& r = *itr_lower;

		// Bifurcate the selected region and replace the selected region with
		// the two new ones.
		const std::pair<R,R>& bif_out = r.bifurcate(x);
		typename std::set<R>::iterator itr = _regions.find(r);
		_regions.erase(itr);
		_regions.insert(bif_out.first);
		_regions.insert(bif_out.second);
		recache();

		log_bdd_hist.push_back(rejection_bound(true));
	}

	return Rcpp::NumericVector(log_bdd_hist.begin(), log_bdd_hist.end());
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::refine(unsigned int N, double tol,
	bool greedy, unsigned int report)
{
	if (tol < 0) {
		Rcpp::stop("tol must be nonnegative");
	}

	std::vector<double> log_bdd_hist;
	log_bdd_hist.push_back(rejection_bound(true));

	for (unsigned int j = 0; j < N; j++) {
		// If we can beat the tolerance before we reach N steps, return now
		if (log_bdd_hist[j] <= std::log(tol)) {
			break;
		}

		unsigned int L = _regions.size();

		// Each region's contribution to the rejection rate
		Rcpp::NumericVector log_volume = rejection_bound_regions(true);

		// Identify the regions which are bifurcatable; for the rest, we set
		// their log-volume to -inf so they will not be selected.
		unsigned int n_bif = 0;
		for (unsigned int l = 0; l < L; l++) {
			log_volume(l) = (*_bifurcatable)(l) ? log_volume(l) : R_NegInf;
			n_bif += (*_bifurcatable)(l);
		}

		if (Rcpp::any(Rcpp::is_nan(log_volume))) {
			Rcpp::stop("nan found in log_volume");
		}

		if (n_bif == 0) {
			Rcpp::warning("No regions left to bifurcate");
			break;
		}

		unsigned int jdx;
		if (greedy) {
			jdx = Rcpp::which_max(log_volume);
		} else {
			jdx = r_categ(log_volume, true);
		}

		const R& r = _regions_vec[jdx];

		// Split the target region and make another proposal with it.
		//
		// TBD: we may not need to recache each time through the loop.
		// Perhaps we can insert into the end and zero out the removed entries?
		const std::pair<R,R>& bif_out = r.bifurcate();
		typename std::set<R>::iterator itr = _regions.find(r);
		_regions.erase(itr);
		_regions.insert(bif_out.first);
		_regions.insert(bif_out.second);
		recache();

		log_bdd_hist.push_back(rejection_bound(true));

		if (j % report == 0 && report < uint_max) {
			logger("After %d steps log Pr{rejection} <= %g\n", j, log_bdd_hist[j+1]);
		}
	}

	return Rcpp::NumericVector(log_bdd_hist.begin(), log_bdd_hist.end());
}

template <class T, class R>
FMMProposal<T,R>::FMMProposal(const std::vector<R>& regions)
: _regions(), _regions_vec()
{
	_regions.insert(regions.begin(), regions.end());
	recache();
}

template <class T, class R>
void FMMProposal<T,R>::recache()
{
	/*
	* TBD: This might be a good place to check for no overlaps, and perhaps to
	* make sure there are no gaps?
	*/

	unsigned int N = _regions.size();
	_regions_vec.clear();
	_log_xi_upper = std::make_unique<Rcpp::NumericVector>(N);
	_log_xi_lower = std::make_unique<Rcpp::NumericVector>(N);
	_bifurcatable = std::make_unique<Rcpp::LogicalVector>(N);

	unsigned int j = 0;
	for (auto itr = _regions.begin(); itr != _regions.end(); ++itr) {
		_regions_vec.push_back(*itr);
		(*_log_xi_upper)[j] = itr->get_xi_upper(true);
		(*_log_xi_lower)[j] = itr->get_xi_lower(true);
		(*_bifurcatable)[j] = itr->is_bifurcatable();
		j++;
	}
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::get_xi_upper(bool log) const
{
	Rcpp::NumericVector out(_log_xi_upper->begin(), _log_xi_upper->end());
	if (log) { return out; } else { return Rcpp::exp(out); }
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::get_xi_lower(bool log) const
{
	Rcpp::NumericVector out(_log_xi_lower->begin(), _log_xi_lower->end());
	if (log) { return out; } else { return Rcpp::exp(out); }
}

template <class T, class R>
Rcpp::LogicalVector FMMProposal<T,R>::get_bifurcatable() const
{
	Rcpp::LogicalVector out(_bifurcatable->begin(), _bifurcatable->end());
	return out;
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::get_pi(bool log) const
{
	Rcpp::NumericVector lxu(_log_xi_upper->begin(), _log_xi_upper->end());
	const Rcpp::NumericVector& out = lxu - vws::log_sum_exp(lxu);
	if (log) { return out; } else { return Rcpp::exp(out); }
}


template <class T, class R>
double FMMProposal<T,R>::rejection_bound(bool log) const
{
	// Each region's contribution to the rejection rate bound
	const Rcpp::NumericVector& lxl = *_log_xi_lower;
	const Rcpp::NumericVector& lxu = *_log_xi_upper;
	const Rcpp::NumericVector& log_bound = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu);

	// Overall rejection rate bound
	double out = log_sum_exp(log_bound);
	if (log) { return out; } else { return exp(out); }
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::rejection_bound_regions(bool log) const
{
	// Each region's contribution to the rejection rate bound.
	const Rcpp::NumericVector& lxl = *_log_xi_lower;
	const Rcpp::NumericVector& lxu = *_log_xi_upper;
	const Rcpp::NumericVector& out = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu);
	if (log) { return out; } else { return exp(out); }
}

template <class T, class R>
double FMMProposal<T,R>::nc(bool log) const
{
	const Rcpp::NumericVector& lxu = *_log_xi_upper;
	double out = log_sum_exp(lxu);
	return log ? out : exp(out);
}

template <class T, class R>
typename std::vector<T> FMMProposal<T,R>::r(unsigned int n) const
{
	return r_ext(n).first;
}

template <class T, class R>
typename std::pair<std::vector<T>, std::vector<unsigned int>>
FMMProposal<T,R>::r_ext(unsigned int n) const
{
	// Draw from the mixing weights, which are given on the log scale and not
	// normalized.
	const Rcpp::NumericVector& lxu = *_log_xi_upper;
	const Rcpp::IntegerVector& idx = r_categ(n, lxu, true);

	// Draw the values from the respective mixture components.
	std::vector<T> x;
	for (unsigned int i = 0; i < n; i++) {
		const std::vector<T>& draws = _regions_vec[idx(i)].r(1);
		x.push_back(draws[0]);
	}

	return std::make_pair(x, Rcpp::as<std::vector<unsigned int>>(idx));
}

template <class T, class R>
double FMMProposal<T,R>::d(const T& x, bool normalize, bool log) const
{
	double log_nc = normalize ? nc(true) : 0;

	/*
	* Search for the region containing x using std::set. This should require
	* something like O(log N) time, where N is the number of regions.
	*
	* 1. Create a singleton region `x_singleton` from the point x.
	* 2. Find the region r in  _regions whose upper bound is x. This should be
	*    one position past the region we want.
	* 3. Rewind the iterator one step. This should yield the region with x.
	* 4. Make sure this region contains x; otherwise, throw an error.
	*/
	const R& x_singleton = _regions.begin()->singleton(x);

	typename std::set<R>::const_iterator itr_lower = _regions.upper_bound(x_singleton);
	--itr_lower;

	double out;

	if (itr_lower == _regions.end()) {
		// Could not find region which is upper bound for x. Assume that x is
		// outside of the support
		out = R_NegInf;
	} else if (!itr_lower->s(x)) {
		// Could not find a region that contains x. Assume that x is outside of
		// the support.
		out = R_NegInf;
	} else {
		out = itr_lower->w_major(x, true) + itr_lower->d_base(x, true) - log_nc;
	}

	return log ? out : exp(out);
}

template <class T, class R>
double FMMProposal<T,R>::d_target_unnorm(const T& x, bool log) const
{
	double out = _regions.begin()->w(x, true) + _regions.begin()->d_base(x, true);
	return log ? out : exp(out);
}

template <class T, class R>
double FMMProposal<T,R>::w_major(const T& x, bool log) const
{
	/*
	* Search for the region containing x using std::set, similar to the `d`
	* method.
	*/
	const R& x_singleton = _regions.begin()->singleton(x);

	typename std::set<R>::const_iterator itr_lower = _regions.upper_bound(x_singleton);
	--itr_lower;

	double out;

	if (itr_lower == _regions.end()) {
		// Could not find region which is upper bound for x. Assume that x is
		// outside of the support
		out = R_NegInf;
	} else if (!itr_lower->s(x)) {
		// Could not find a region that contains x. Assume that x is outside of
		// the support.
		out = R_NegInf;
	} else {
		out = itr_lower->w_major(x, true);
	}

	return log ? out : exp(out);
}

template <class T, class R>
Rcpp::DataFrame FMMProposal<T,R>::summary() const
{
	unsigned int N = _regions.size();
	Rcpp::StringVector v1(N);
	Rcpp::NumericVector v2(N);
	Rcpp::NumericVector v3(N);

	unsigned int j = 0;

	for (auto itr = _regions.begin(); itr != _regions.end(); ++itr) {
		v1[j] = itr->description();
		v2[j] = itr->get_xi_upper(true);
		v3[j] = itr->get_xi_lower(true);
		j++;
	}

	return Rcpp::DataFrame::create(
		Rcpp::Named("Region") = v1,
		Rcpp::Named("log_xi_upper") = v2,
		Rcpp::Named("log_xi_lower") = v3
	);
}

template <class T, class R>
void FMMProposal<T,R>::print(unsigned int n) const
{
	const Rcpp::DataFrame& tbl = summary();
	unsigned int N = tbl.nrows();

	const Rcpp::IntegerVector& idx = Rcpp::seq(0, std::min(n,N) - 1);

	const Rcpp::StringVector& v1 = tbl["Region"];
	const Rcpp::NumericVector& v2 = tbl["log_xi_upper"];
	const Rcpp::NumericVector& v3 = tbl["log_xi_lower"];

	const Rcpp::DataFrame& tbl_head = Rcpp::DataFrame::create(
		Rcpp::Named("Region") = v1[idx],
		Rcpp::Named("log_xi_upper") = v2[idx],
		Rcpp::Named("log_xi_lower") = v3[idx]
	);

	Rprintf("FMM Proposal with %d regions\n", N);

	Rcpp::print(tbl_head);

	if (N > n) {
		Rprintf("There are %d more regions not displayed\n", N - n);
	}
}

}

#endif
