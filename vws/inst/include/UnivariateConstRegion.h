#ifndef VWS_UNIVARIATE_CONST_REGION_H
#define VWS_UNIVARIATE_CONST_REGION_H

#include <Rcpp.h>
#include <memory>
#include "fntl.h"
#include "Region.h"
#include "UnivariateHelper.h"
#include "optimize-hybrid.h"
#include "log-sum-exp.h"
#include "typedefs.h"

namespace vws {

/*
* Univariate Region with Constant Majorizer
*
* A subclass of Region based on univariate intervals and a constant majorizer
* for the weight function. This version is for continuous supports.
*/
class UnivariateConstRegion : public Region<double>
{
protected:
	double _a;
	double _b;
	const uv_weight_function* _w;
	const UnivariateHelper<double>* _helper;
	double _log_w_max;
	double _log_w_min;
	double _log_prob;

public:
	/*
	* Construct a singleton region based on interval $(a,a]$.
	* - `a`: Lower and upper limit of interval.
	* - `w`: Weight function for the target distribution.
	* - `helper`: a subclass of `UnivariateHelper`.
	*/
	UnivariateConstRegion(double a, const uv_weight_function& w,
		const UnivariateHelper<double>& helper);

	/*
	* Construct a region based on interval $(a,b]$.
	* - `a` Lower limit of interval.
	* - `b` Upper limit of interval.
	* - `w` Weight function for the target distribution.
	* - `helper`: a subclass of `UnivariateHelper`.
	*/
	UnivariateConstRegion(double a, double b, const uv_weight_function& w,
		const UnivariateHelper<double>& helper);

	/*
	* The following functions override abstract methods in `Region`. See that
	* class' documentation for their interfaces.
	*/
	double d_base(const double& x, bool log = false) const;
	std::vector<double> r(unsigned int n) const;
	double d(const double& x, bool log = false) const;
	bool s(const double& x) const;
	double w(const double& x, bool log = true) const;
	double w_major(const double& x, bool log = true) const;
	bool is_bifurcatable() const;
	double get_lower() const { return _a; }
	double get_upper() const { return _b; }
	double get_xi_upper(bool log = true) const;
	double get_xi_lower(bool log = true) const;
	std::string description() const;
	void print() const;

	/*
	* Maximize or minimize the function $w$ over this region. Optimization
	* is carried out with the `optimize_hybrid` function.
	*
	* - `maximize`: if `true` do maximization; otherwise do minimization.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*
	* Returns the optimized value of $w$.
	*/
	double optimize(bool maximize = true, bool log = true) const;


	/*
	* A midpoint between limits $a$ and $b$ of region. If $a$ and $b$ are both
	* finite, return the standard midpoint. If both are infinite, return zero.
	* If only $a$ is finite, return a larger point in the support. If only $b$
	* is finite, return a smaller point in the support.
	*/
	double midpoint() const;

	/*
	* Return a pair of regions that result from bifurcating this region. The
	* bifurcation point is chosen to be the midpoint of $(a, b]$.
	*/
	std::pair<UnivariateConstRegion,UnivariateConstRegion> bifurcate() const;

	/*
	* Return a pair of regions that result from bifurcating this region at $x$.
	*/
	std::pair<UnivariateConstRegion,UnivariateConstRegion> bifurcate(const double& x) const;

	/*
	* Return a region based on the singleton interval $(x, x]$, using this
	* object's weight function, base distribution, etc.
	*/
	UnivariateConstRegion singleton(const double& x) const;

	/*
	* Region $(a_1, b_1]$ is considered "less than" $(a_2, b_2]$ if $a_1 < a_2$.
	*/
	bool operator<(const UnivariateConstRegion& x) const {
		return _a < x._a;
	}

	/*
	* Region $(a_1, b_1]$ is considered "equal to" $(a_2, b_2]$ if $a_1 = a_2$
	* and $b_1 = b_2$.
	*/
	bool operator==(const UnivariateConstRegion& x) const {
		return _a == x._a && _b == x._b;
	}

	/*
	* Set this Region to be equal to `x`.
	*/
	const UnivariateConstRegion& operator=(const UnivariateConstRegion& x) {
		_a = x._a;
		_b = x._b;
		_w = x._w;
		_helper = x._helper;
		_log_w_max = x._log_w_max;
		_log_w_min = x._log_w_min;
		_log_prob = x._log_prob;
		return *this;
	}
};

inline UnivariateConstRegion::UnivariateConstRegion(double a,
	const uv_weight_function& w, const UnivariateHelper<double>& helper)
: _a(a), _b(a), _w(&w), _helper(&helper)
{
	_log_w_max = (*_w)(a, true);
	_log_w_min = (*_w)(a, true);
	_log_prob = R_NegInf;
}

inline UnivariateConstRegion::UnivariateConstRegion(double a, double b,
	const uv_weight_function& w, const UnivariateHelper<double>& helper)
: _a(a), _b(b), _w(&w), _helper(&helper)
{
	if (a > b) {
		Rcpp::stop("a > b");
	}

	_log_w_max = optimize(true);
	_log_w_min = optimize(false);

	// Compute $P(a < X <= b)$ for $X \sim g$ on the log scale.
	_log_prob = log_sub2_exp(_helper->p(_b, true, true), _helper->p(_a, true, true));
}

inline double UnivariateConstRegion::d_base(const double& x, bool log) const
{
	return _helper->d(x, log);
}

inline std::vector<double> UnivariateConstRegion::r(unsigned int n) const
{
	// Generate a draw from $g_j$; i.e., the density $g$ truncated to this
	// region. Compute g$q((pb - pa) * u + pa) on the log scale.
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	double log_pa = _helper->p(_a, true, true);
	const Rcpp::NumericVector& log_p = log_add2_exp(_log_prob + log(u), Rcpp::rep(log_pa, n));

	std::vector<double> out;
	for (unsigned int i = 0; i < n; i++) {
		out.push_back(_helper->q(log_p(i), true, true));
	}

	return out;
}

inline double UnivariateConstRegion::d(const double& x, bool log) const
{
	double out;
	if (!s(x)) {
		out = R_NegInf;
	} else {
		out = _helper->d(x, true) - log_sub2_exp(_helper->p(_b, true, true), _helper->p(_a, true, true));
	}
	return log ? out : exp(out);
}

inline bool UnivariateConstRegion::s(const double& x) const
{
	return (_a < x && x <= _b) && _helper->s(x);
}

inline double UnivariateConstRegion::w(const double& x, bool log) const
{
	return (*_w)(x, log);
}

inline double UnivariateConstRegion::w_major(const double& x, bool log) const
{
	double out = _helper->s(x) ? _log_w_max : R_NegInf;
	return log ? out : exp(out);
}

inline double UnivariateConstRegion::midpoint() const
{
	double out;

	if (std::isinf(_a) && std::isinf(_b) && _a < 0 && _b > 0) {
		// In this case, we have an interval (-inf, inf). Make a split at zero.
		out = 0;
	} else if (std::isinf(_a) && _a < 0) {
		// Left endpoint is -inf. Split based on right endpoint.
		out = _b - std::fabs(_b) - 1;
	} else if (std::isinf(_b) && _b > 0) {
		// Right endpoint is inf. Split based on left endpoint.
		out = _a + std::fabs(_a) + 1;
	} else {
		out = (_a + _b) / 2;
	}

	return out;
}

inline std::pair<UnivariateConstRegion,UnivariateConstRegion>
UnivariateConstRegion::bifurcate() const
{
	return bifurcate(midpoint());
}

inline std::pair<UnivariateConstRegion,UnivariateConstRegion>
UnivariateConstRegion::bifurcate(const double& x) const
{
	UnivariateConstRegion r1(_a, x, *_w, *_helper);
	UnivariateConstRegion r2(x, _b, *_w, *_helper);
	return std::make_pair(r1, r2);
}

inline UnivariateConstRegion UnivariateConstRegion::singleton(const double& x) const
{
	return UnivariateConstRegion(x, *_w, *_helper);
}

inline bool UnivariateConstRegion::is_bifurcatable() const
{
	return true;
}

inline double UnivariateConstRegion::get_xi_upper(bool log) const
{
	double log_xi_upper = _log_w_max + _log_prob;
	return log ? log_xi_upper : exp(log_xi_upper);
}

inline double UnivariateConstRegion::get_xi_lower(bool log) const
{
	double log_xi_lower = _log_w_min + _log_prob;
	return log ? log_xi_lower : exp(log_xi_lower);
}

inline std::string UnivariateConstRegion::description() const
{
	char buf[32];
	sprintf(buf, "(%g, %g]", _a, _b);
	return buf;
}

inline void UnivariateConstRegion::print() const
{
	Rprintf("Region<double> (%g, %g]\n", _a, _b);
}

inline double UnivariateConstRegion::optimize(bool maximize, bool log) const
{
	// Pass the log-weight function to `optimize_hybrid`.
    const fntl::dfd& f = [&](double x) -> double { return (*_w)(x, true); };
	const auto& out = optimize_hybrid(f, 0, _a, _b, maximize);

	if ( (out.value > 0) && std::isinf(out.value) ) {
		Rcpp::stop("Infinite value found in optimize. Cannot be used with UnivariateConstRegion");
	}

	return log ? out.value : exp(out.value);
}

}

#endif
