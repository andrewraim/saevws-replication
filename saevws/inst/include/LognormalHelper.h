#ifndef LOGNORMAL_HELPER_H
#define LOGNORMAL_HELPER_H

#include "vws.h"

/*
* Subclass of `LognormalHelper` for the Uniform distribution.
*/

class LognormalHelper : public vws::UnivariateHelper<double>
{
public:
	LognormalHelper(double mean, double sd)
	: _mean(mean), _sd(sd)
	{
	}

	double d(double x, bool log = false) const {
		return R::dlnorm(x, _mean, _sd, log);
	}
	double p(double q, bool lower = true, bool log = false) const {
		return R::plnorm(q, _mean, _sd, lower, log);
	}
	double q(double p, bool lower = true, bool log = false) const {
		return R::qlnorm(p, _mean, _sd, lower, log);
	}
	bool s(double x) const {
		return 0 <= x;
	}
	const LognormalHelper& operator=(const LognormalHelper& x) {
		_mean = x._mean;
		_sd = x._sd;
		return *this;
	}

private:
	double _mean;
	double _sd;
};

#endif

