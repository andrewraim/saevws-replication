library(tidyverse)
library(xtable)
library(saevws)
library(mcmcse)
library(coda)

# Set a seed
set.seed(1234)

# Quantiles for ESS
probs = c(0.000, 0.010, 0.025)

# Significance for CI width
alpha = 0.10

# Args for VWS
tol_suff_levels = c(0.75, 0.85)
tol_merge_levels = c(0.001, 0.01)

tol_levels = expand.grid(
	idx_merge = seq_along(tol_merge_levels),
	idx_suff = seq_along(tol_suff_levels)) %>%
	mutate(tol_suff = tol_suff_levels[idx_suff]) %>%
	mutate(tol_merge = tol_merge_levels[idx_merge])

run_vws0 = FALSE

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

# ----- Initial values -----

# Pick starting values based on the observed data.
lm1_out = lm(y ~ X - 1)
lm2_out = lm(log(s2) ~ Z - 1)
beta_init = coef(lm1_out)
gamma_init = coef(lm2_out)
phi2_init = sigma(lm1_out)^2
tau2_init = sigma(lm2_out)^2
init = init_joint(m, d1, d2, beta = beta_init, gamma = gamma_init, sigma2 = s2,
	phi2 = phi2_init, tau2 = tau2_init)

# Run VWS0
inner_ctrl = control_inner(tol_suff = 0.85, tol_merge = 0.001,
	max_rejects = 1e6, method = "vws-basic", N = 50)
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 1,
	inner = inner_ctrl, save_latent = seq_len(m))
vws0_out = gibbs_joint(y, s2, X, Z, df, init, control)
print(vws0_out)
