#ifndef VWS_REGION_H
#define VWS_REGION_H

#include <Rcpp.h>
#include <memory>

namespace vws {

/*
* A Region contains all of the problem-specific logic for a VWS sampler.
* This is an abstract class that defines the interface.
*/
template <class T>
class Region
{
public:
	/*
	* Density function $g$ of the base distribution.
	* - `x`: argument of density.
	* - `log` if `true`, return value on the log-scale. Otherwise, return it on
	*   the original scale.
	*/
	virtual double d_base(const T& x, bool log = false) const = 0;

	/*
	* Generate a vector of $n$ draws from $g_j$ specific to this region.
	*/
	virtual std::vector<T> r(unsigned int n) const = 0;

	/*
	* Density of $g_j$ specific to this region.
	* - `x`: argument of density.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	virtual double d(const T& x, bool log = false) const = 0;

	/*
	* Indicator of whether $x$ is in the support for $g_j$ specific to this
	* region.
	*/
	virtual bool s(const T& x) const = 0;

	/*
	* The weight function $w$.
	* - `x`: argument.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	virtual double w(const T& x, bool log = true) const = 0;

	/*
	* Majorized weight function $\overline{w}_j$ for this region.
	* - `x`: Argument to weight function.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it on
	*   the original scale.
	*/
	virtual double w_major(const T& x, bool log = true) const = 0;

	/*
	* Indicator of whether this region is bifurcatable.
	*/
	virtual bool is_bifurcatable() const = 0;

	/*
	* The quantity $\overline{\xi}_j$ for this region.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	virtual double get_xi_upper(bool log = true) const = 0;

	/*
	* The quantity $\underline{\xi}_j$ for this region.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	virtual double get_xi_lower(bool log = true) const = 0;

	/*
	* A string that describes this region.
	*/
	virtual std::string description() const = 0;

	/*
	* Print a description of this region.
	*/
	virtual void print() const = 0;
};

}

#endif
