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

# ----- Fit the model -----
init = init_mismatch(m, d = ncol(X), mu = mu_true)
inner = control_inner(method = "arms")
control = control_mismatch(R = 3000, burn = 1000, thin = 1, report = 100,
	save_latent = 1:m, inner = inner)
fixed = fixed_mismatch(mu = FALSE)
gibbs_out = gibbs_mismatch(y, sigma, X, init, control, fixed)
print(gibbs_out)

plot(gibbs_out$beta_hist[,1], type = "l")
plot(gibbs_out$beta_hist[,2], type = "l")
plot(gibbs_out$tau2_hist, type = "l")

i = 4
plot(gibbs_out$mu_hist[,i], type = "l")
abline(h = mu_true[i], lty = 2, col = "red")
