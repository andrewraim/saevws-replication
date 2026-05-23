library(saevws)

set.seed(1234)

# ----- Generate data from the model with known parameters -----
m = 500
X = cbind(1, rnorm(m))
sigma = rgamma(m, 1.25, 20) |> sqrt()

beta_true = c(1, -1)
Xbeta_true = X %*% beta_true
tau_true = 0.25
mu_true = rlnorm(m, Xbeta_true, tau_true)
y = rnorm(m, mu_true, tau_true)

# ----- Fit the model using IMH -----
init = init_mismatch(m, d = ncol(X), mu = mu_true)
inner = control_inner(method = "imh")
control = control_mismatch(R = 30000, burn = 20000, thin = 1, report = 1000,
	save_latent = 1:m, inner = inner)
fixed = fixed_mismatch(mu = FALSE)
gibbs_imh = gibbs_mismatch(y, sigma, X, init, control, fixed)
print(gibbs_imh)

plot(gibbs_imh$beta_hist[,1], type = "l")
plot(gibbs_imh$beta_hist[,2], type = "l")
plot(gibbs_imh$tau2_hist, type = "l")

ess_imh = ess(gibbs_imh$mu_hist)
hist(ess_imh)
i = which.min(ess_imh)

plot(gibbs_imh$mu_hist[,i], type = "l")
abline(h = mu_true[i], lty = 2, col = "red")

# ----- Fit the model using ARMS -----
init = init_mismatch(m, d = ncol(X), mu = mu_true)
inner = control_inner(method = "arms")
control = control_mismatch(R = 3000, burn = 1000, thin = 1, report = 100,
	save_latent = 1:m, inner = inner)
fixed = fixed_mismatch(mu = FALSE)
gibbs_arms = gibbs_mismatch(y, sigma, X, init, control, fixed)
print(gibbs_arms)

plot(gibbs_arms$beta_hist[,1], type = "l")
plot(gibbs_arms$beta_hist[,2], type = "l")
plot(gibbs_arms$tau2_hist, type = "l")

ess_arms = ess(gibbs_arms$mu_hist)
hist(ess_arms)
i = which.min(ess_arms)

plot(gibbs_arms$mu_hist[,i], type = "l")
abline(h = mu_true[i], lty = 2, col = "red")
