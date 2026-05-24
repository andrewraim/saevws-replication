#ifndef ARMS_MISMATCH_FUNCTOR
#define ARMS_MISMATCH_FUNCTOR

#include "local-util.h"

class arms_mismatch_functor {
public:
	arms_mismatch_functor(double y, double sigma, double xbeta, double tau)
		: _y(y), _sigma(sigma), _xbeta(xbeta), _tau(tau), _nEvaluations(0)
	{
	}

	double operator()(double x)
	{
		++_nEvaluations;
		double out1 = R::dnorm(x, _y, _sigma, true);
		double out2 = R::dlnorm(x, _xbeta, _tau, true);
		return out1 + out2;
	}

	int nEvaluations() const { return _nEvaluations; }

private:
	double _y;
	double _sigma;
	double _xbeta;
	double _tau;
	int _nEvaluations;
};

#endif
