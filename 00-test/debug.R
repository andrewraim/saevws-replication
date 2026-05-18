Rcpp::sourceCpp("../01-sim-cond/samplers.cpp")

# ----- Data setup -----
ff = file.path("..", "data", "saipe.csv")
dat_saipe = read_csv(ff) %>%
	mutate(df = 0.36 * sqrt(hu_sampled)) %>%
	mutate(y = log(pov_count)) %>%
	mutate(s2 = pov_se^2 / pov_count^2) %>%
	filter(pov_count > 1 & df > 1) %>%
	filter(!is.na(snap)) # Removing one county with missing snap

m = nrow(dat_saipe)
X = model.matrix(~ log1p(snap) + log1p(pep), data = dat_saipe)
Z = model.matrix(~ log(hu_sampled) , data = dat_saipe)
df = dat_saipe$df
y = dat_saipe$y
s2 = dat_saipe$s2
d1 = ncol(X)
d2 = ncol(Z)

theta = y
kappa = (df - 1) / 2.0
lambda = (y - theta)^2 / 2.0 + df * s2 / 2.0

gamma = c(1.8477, -0.9028)
tau = 0.3062679

Zgamma = Z %*% gamma


# ----- Focus on one obs -----

## This obs seems to have very low df
## This is influcing the invgamma distribution so that it is positive on larger
## supports. We tend to draw large numbers from it which are not accepted in
## the MH step.
idx = 2655

tol_suff = 0.85
tol_merge = 1e-4

out1 = r_metro(n = 2000, init = 0, Zgamma[idx], tau, kappa[idx], lambda[idx])
out2 = r_target(n = 2000, Zgamma[idx], tau, kappa[idx], lambda[idx], tol_suff, tol_merge, max_rejects)

plot(out1$draws[10:2000], type = "l")
plot(out2$draws[10:2000], type = "l")

# Metropolis step
sigma2 = 0.3
#sigma2_prop = raim::r_invgamma(1, kappa[idx], lambda[idx])
sigma2_prop = 0.35
log_num = dlnorm(sigma2_prop, Zgamma[idx], tau, log = TRUE)
log_den = dlnorm(sigma2, Zgamma[idx], tau, log = TRUE)
log_ratio = min(log_num - log_den, 0)

curve(raim::d_invgamma(x, kappa[idx], lambda[idx], log = TRUE), xlim = c(0, 10))
curve(dlnorm(x, Zgamma[idx], tau, log = TRUE), xlim = c(0, 10))


lr_fn = function(x) { -(log(x) - Zgamma[idx])^2 / (2 * tau^2) - log(x) }
curve(lr_fn, xlim = c(0, 10))

curve(raim::d_invgamma(x, kappa[idx], lambda[idx], log = TRUE), xlim = c(0, 10))
curve(raim::p_invgamma(x, kappa[idx], lambda[idx], log = FALSE), xlim = c(0, 10))

f_target = function(x, log = FALSE) {
	out = dlnorm(x, Zgamma[idx], tau, log = TRUE) +
		raim::d_invgamma(x, kappa[idx], lambda[idx], log = TRUE)
	if (log) { return(out) } else { return(exp(out)) }
}
n_target = integrate(f_target, lower = 0, upper = Inf)$value

p_target = function(x, log = FALSE) {
	integrate_out = integrate(f_target, lower = 0, upper = x)
	out = log(integrate_out$value) - log(n_target)
	if (log) { return(out) } else { return(exp(out)) }
}

xx = 0.5643967
curve(f_target)
p_target(xx)

# The conditional distribution is focused on the interval (0.1103573, 0.5643967).
# The proposal distribution is heavily focused on larger numbers:
# P(X > 0.565) = 0.9352036. The lognormal we get after the ratio cancels out is
# effectively below 0.565 as well: P(X < 0.565) = 0.9999767
raim::p_invgamma(xx, kappa[idx], lambda[idx], lower.tail = FALSE, log = FALSE)
plnorm(xx, Zgamma[idx], tau)
