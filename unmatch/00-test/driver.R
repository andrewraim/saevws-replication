library(saevws)
library(mcmcse)

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
init = init_unmatch(m, d = ncol(X), mu = mu_true)
inner = control_inner(method = "imh")
control = control_unmatch(R = 30000, burn = 20000, thin = 1, report = 5000,
	save_latent = 1:m, inner = inner)
fixed = fixed_unmatch(mu = FALSE)
gibbs_imh = gibbs_unmatch(y, sigma, X, init, control, fixed)
print(gibbs_imh)

plot(gibbs_imh$beta_hist[,1], type = "l")
plot(gibbs_imh$beta_hist[,2], type = "l")
plot(gibbs_imh$tau2_hist, type = "l")

ess_imh = ess(gibbs_imh$mu_hist)
hist(ess_imh)

i = which.min(ess_imh)
plot(gibbs_imh$mu_hist[,i], type = "l")
abline(h = mu_true[i], lty = 2, col = "red")

# ----- Fit the model using AM -----
init = init_unmatch(m, d = ncol(X), mu = mu_true)
inner = control_inner(method = "am", am_varprop_init = 25, am_varprop_eps = 0.001)
control = control_unmatch(R = 10000, burn = 8000, thin = 1, report = 5000,
	save_latent = 1:m, inner = inner)
fixed = fixed_unmatch(mu = FALSE)
gibbs_am = gibbs_unmatch(y, sigma, X, init, control, fixed)
print(gibbs_am)

plot(gibbs_am$beta_hist[,1], type = "l")
plot(gibbs_am$beta_hist[,2], type = "l")
plot(gibbs_am$tau2_hist, type = "l")

ess_am = ess(gibbs_am$mu_hist)
hist(ess_am)

i = which.min(ess_am)
plot(gibbs_am$mu_hist[,i], type = "l")
abline(h = mu_true[i], lty = 2, col = "red")

# ----- Fit the model using ARMS -----
init = init_unmatch(m, d = ncol(X), mu = mu_true)
inner = control_inner(method = "arms")
control = control_unmatch(R = 3000, burn = 1000, thin = 1, report = 100,
	save_latent = 1:m, inner = inner)
fixed = fixed_unmatch(mu = FALSE)
gibbs_arms = gibbs_unmatch(y, sigma, X, init, control, fixed)
print(gibbs_arms)

plot(gibbs_arms$beta_hist[,1], type = "l")
plot(gibbs_arms$beta_hist[,2], type = "l")
plot(gibbs_arms$tau2_hist, type = "l")

ess_arms = ess(gibbs_arms$mu_hist)
hist(ess_arms)

i = which.min(ess_arms)
plot(gibbs_arms$mu_hist[,i], type = "l")
abline(h = mu_true[i], lty = 2, col = "red")

# ----- Fit the model using basic VWS -----
init = init_unmatch(m, d = ncol(X), mu = mu_true)
inner = control_inner(method = "vws-basic", tol_suff = 0.85)
control = control_unmatch(R = 3000, burn = 1000, thin = 1, report = 100,
	save_latent = 1:m, inner = inner)
fixed = fixed_unmatch(mu = FALSE)
gibbs_vwsb = gibbs_unmatch(y, sigma, X, init, control, fixed)
print(gibbs_vwsb)

plot(gibbs_vwsb$beta_hist[,1], type = "l")
plot(gibbs_vwsb$beta_hist[,2], type = "l")
plot(gibbs_vwsb$tau2_hist, type = "l")

ess_vwsb = ess(gibbs_vwsb$mu_hist)
hist(ess_vwsb)

i = which.min(ess_vwsb)
plot(gibbs_vwsb$mu_hist[,i], type = "l")
abline(h = mu_true[i], lty = 2, col = "red")

# ----- Fit the model using VWS with self-tuning -----
init = init_unmatch(m, d = ncol(X), mu = mu_true)
inner = control_inner(method = "vws-tune", tol_suff = 0.85, tol_merge = 0.01,
	tune = 10000)
control = control_unmatch(R = 3000, burn = 1000, thin = 1, report = 100,
	save_latent = 1:m, inner = inner)
fixed = fixed_unmatch(mu = FALSE)
gibbs_vwst = gibbs_unmatch(y, sigma, X, init, control, fixed)
print(gibbs_vwst)

plot(gibbs_vwst$beta_hist[,1], type = "l")
plot(gibbs_vwst$beta_hist[,2], type = "l")
plot(gibbs_vwst$tau2_hist, type = "l")

ess_vwst = ess(gibbs_vwst$mu_hist)
hist(ess_vwst)

i = which.min(ess_vwst)
plot(gibbs_vwst$mu_hist[,i], type = "l")
abline(h = mu_true[i], lty = 2, col = "red")

# Plot number of tuned VWS proposals by iteration
data.frame(tuned = gibbs_vwst$mu_tuned_hist) %>%
	mutate(iter = row_number()) %>%
	filter(iter > 20) %>%
	ggplot() +
	geom_line(aes(iter, tuned)) +
	geom_rect(xmin = 0, xmax = 20,  ymin = 0,  ymax = Inf, fill = "red") +
	scale_y_continuous(breaks = 1:50, minor_breaks = NULL) +
	xlab("Iteration") +
	ylab("Number of Tuned VWS Proposals") +
	theme_minimal()

# Plot number of VWS rejections per area by iteration
data.frame(tuned = gibbs_vwst$mu_rejects_hist / m) %>%
	mutate(iter = row_number()) %>%
	filter(iter > 20) %>%
	ggplot() +
	geom_line(aes(iter, tuned)) +
	geom_rect(xmin = 0, xmax = 20, ymin = min(gibbs_vwst$mu_rejects_hist / m),
		ymax = Inf, fill = "red") +
	xlab("Iteration") +
	ylab("Number of VWS Rejections Per Area") +
	theme_minimal()

# ----- Same as above, but stop tuning after an initial period -----
init = init_unmatch(m, d = ncol(X), mu = mu_true)
inner = control_inner(method = "vws-tune", tol_suff = 0.85, tol_merge = 0.01,
	tune = 100)
control = control_unmatch(R = 3000, burn = 1000, thin = 1, report = 100,
	save_latent = 1:m, inner = inner)
gibbs_vwst2 = gibbs_unmatch(y, sigma, X, init, control)
print(gibbs_vwst2)

# Plot number of VWS rejections per area by iteration
data.frame(tuned = gibbs_vwst2$mu_rejects_hist / m) %>%
	mutate(iter = row_number()) %>%
	filter(iter > 20) %>%
	ggplot() +
	geom_line(aes(iter, tuned)) +
	geom_rect(xmin = 0, xmax = 20, ymin = min(gibbs_vwst$mu_rejects_hist / m),
		ymax = Inf, fill = "red") +
	xlab("Iteration") +
	ylab("Number of VWS Rejections Per Area") +
	theme_minimal()
