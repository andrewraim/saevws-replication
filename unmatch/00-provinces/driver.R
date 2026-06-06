library(saevws)
library(mcmcse)

set.seed(1234)

# ----- Generate data from the model with known parameters -----

df = read_csv("../data/undercoverage.csv")

m = nrow(df)
X = model.matrix(~ pop1991, df)
y = df$undercoverage
sigma = df$cv * y
sigma2 = sigma^2

# ----- Fit the model using IMH -----
init = init_unmatch(m, d = ncol(X))
inner = control_inner(method = "imh")
control = control_unmatch(R = 30000, burn = 20000, thin = 1, report = 5000,
	save_latent = 1:m, inner = inner)
gibbs_imh = gibbs_unmatch(y, sigma, X, init, control)
print(gibbs_imh)

plot(gibbs_imh$beta_hist[,1], type = "l")
plot(gibbs_imh$beta_hist[,2], type = "l")
plot(gibbs_imh$tau2_hist, type = "l")

ess_imh = ess(gibbs_imh$mu_hist)
hist(ess_imh)

i = which.min(ess_imh)
plot(gibbs_imh$mu_hist[,i], type = "l")
abline(h = mu_true[i], lty = 2, col = "red")

mu_hat = apply(gibbs_imh$mu_hist, 2, mean)
cv_hat = sigma / mu_hat
df %>% add_column(mu_hat, cv_hat)


# ----- Fit the model using VWS with self-tuning -----
init = init_unmatch(m, d = ncol(X))
inner = control_inner(method = "vws-tune", tol_suff = 0.85, tol_merge = 0.01,
	tune = 10000)
control = control_unmatch(R = 3000, burn = 1000, thin = 1, report = 100,
	save_latent = 1:m, inner = inner)
gibbs_vwst = gibbs_unmatch(y, sigma, X, init, control)
print(gibbs_vwst)

plot(gibbs_vwst$beta_hist[,1], type = "l")
plot(gibbs_vwst$beta_hist[,2], type = "l")
plot(gibbs_vwst$tau2_hist, type = "l")

ess_vwst = ess(gibbs_vwst$mu_hist)
hist(ess_vwst)

i = which.min(ess_vwst)
plot(gibbs_vwst$mu_hist[,i], type = "l")
abline(h = mu_true[i], lty = 2, col = "red")

mu_hat = apply(gibbs_vwst$mu_hist, 2, mean)
cv_hat = sigma / mu_hat
df %>% add_column(mu_hat, cv_hat)

