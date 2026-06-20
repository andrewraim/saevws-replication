library(tidyverse)

source("../shared/functions.R", chdir = TRUE)

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

gamma = c(1.8485, -0.9029)
tau = sqrt(0.0950)

Zgamma = Z %*% gamma


# ----- Focus on one obs -----

## This obs seems to have very low df
## This is influencing the invgamma distribution so that it is positive on larger
## supports. We tend to draw large numbers from it which are not accepted in
## the MH step.
idx = 2655

tol_suff = 0.85
tol_merge = 1e-4
max_rejects = 1e6

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

f0_target = function(x, log = FALSE) {
	out = dlnorm(x, Zgamma[idx], tau, log = TRUE) +
		raim::d_invgamma(x, kappa[idx], lambda[idx], log = TRUE)
	if (log) { return(out) } else { return(exp(out)) }
}
n_target = integrate(f0_target, lower = 0, upper = Inf)$value
f_target = function(x, log = FALSE) {
	out = f0_target(x, log = TRUE) - log(n_target)
	if (log) { return(out) } else { return(exp(out)) }
}

p_target = function(x, log = FALSE) {
	integrate_out = integrate(f0_target, lower = 0, upper = x)
	out = log(integrate_out$value) - log(n_target)
	if (log) { return(out) } else { return(exp(out)) }
}

q_target = function(p) {
	ff = function(x) { p_target(x) - p }
	root_out = uniroot(ff, interval = c(0.0001, 10))
	return(root_out$root)
}

# xx = 0.5643967
# xx = q_target(0.99)
xx = qlnorm(0.999, meanlog = Zgamma[idx], sdlog = tau)
curve(f_target)
p_target(xx)


# The conditional distribution is focused on the interval (0.1103573, 0.5643967).
# The proposal distribution is heavily focused on larger numbers:
# P(X > 0.565) = 0.9352036. The lognormal we get after the ratio cancels out is
# effectively below 0.565 as well: P(X < 0.565) = 0.9999767
raim::p_invgamma(xx, kappa[idx], lambda[idx], lower.tail = FALSE, log = FALSE)
plnorm(xx, Zgamma[idx], tau)

g = ggplot() +
	geom_function(fun = dlnorm,
		args = list(meanlog = Zgamma[idx], sdlog = tau, log = FALSE),
		n = 500, lty = 1) +
	geom_function(fun = f_target, n = 500, lty = 3) +
	scale_x_continuous(limits = c(0, 1)) +
	geom_vline(xintercept = xx, col = "red", lty = 2) +
	xlab(expression(sigma[i]^2)) +
	ylab("Density") +
	theme_minimal()
ggsave("density.pdf", g, width = 3.5, height = 2.5)

g = ggplot() +
	geom_function(fun = raim::p_invgamma, args = list(a = kappa[idx], b = lambda[idx])) +
	scale_x_continuous(limits = c(0, 0.5)) +
	scale_y_continuous(breaks = seq(0, 1, 0.01)) +
	geom_vline(xintercept = xx, col = "red", lty = 2) +
	xlab(expression(sigma[i]^2)) +
	ylab("CDF") +
	theme_minimal()
ggsave("proposal.pdf", g, width = 3.5, height = 2.5)

