#ifndef ARMS_FUNCTOR
#define ARMS_FUNCTOR

class ARMSJointFunctor {
public:
	ARMSJointFunctor(double mu, double tau, double kappa, double lambda)
		: _mu(mu), _tau(tau), _kappa(kappa), _lambda(lambda), _nEvaluations(0)
	{
	}

	double operator()(double x)
	{
		++_nEvaluations;
		double out1 = d_invgamma(x, _kappa, _lambda, true);
		double out2 = R::dlnorm(x, _mu, _tau, true);
		return out1 + out2;
	}

	int nEvaluations() const { return _nEvaluations; }

private:
	double _mu;
	double _tau;
	double _kappa;
	double _lambda;
	int _nEvaluations;
};

#endif

