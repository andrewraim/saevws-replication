Rcpp::sourceCpp("samplers.cpp")

d_invgamma = function(x, shape, scale, log = FALSE)
{
	n = length(x)
	out = rep(-Inf, n)

	lnc = base::lgamma(shape) - shape * base::log(scale)

	for (i in 1:n) {
		if (x[i] > 0) {
			out[i] = -(shape + 1)*base::log(x[i]) - scale / x[i] - lnc
		}
	}
	if (log) { return(out) } else { return(exp(out)) }
}

d_target_unnorm = function(x, mu, tau, kappa, lambda, log = FALSE)
{
	n = length(x)
	out = rep(-Inf, n)
	for (i in 1:n) {
		if (x[i] > 0) {
			out[i] = d_invgamma(x[i], kappa, lambda, log = TRUE) +
				dlnorm(x[i], mu, tau, log = TRUE)
		}
	}
	if (log) { return(out) } else { return(exp(out)) }
}

n_target = function(mu, tau, kappa, lambda, log = FALSE)
{
	f = function(x) { d_target_unnorm(x, mu, tau, kappa, lambda) }
	int_out = integrate(f, lower = 0, upper = Inf)
	out = base::log(int_out$value)
	if (log) { return(out) } else { return(exp(out)) }
}

d_target = function(x, mu, tau, kappa, lambda, log = FALSE)
{
	lnc = n_target(mu, tau, kappa, lambda, log = TRUE)
	out = d_target_unnorm(x, mu, tau, kappa, lambda, log = TRUE) - lnc
	if (log) { return(out) } else { return(exp(out)) }
}

w_target = function(x, kappa, lambda, log = TRUE)
{
	n = length(x)
	out = rep(-Inf, n)
	for (i in 1:n) {
		if (x[i] > 0) {
			out[i] = d_invgamma(x[i], kappa, lambda, log = TRUE)
		}
	}
	if (log) { return(out) } else { return(exp(out)) }
}

