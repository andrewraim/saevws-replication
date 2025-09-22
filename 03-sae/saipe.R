library(tidyverse)
library(xtable)
library(saevws)
library(mcmcse)

# Set a seed
set.seed(1234)

# Quantiles for ESS
probs = c(0.000, 0.010, 0.025)

# Significance for CI width
alpha = 0.10

# Args for VWG
tol1 = 0.85
tol2 = 0.0001

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
init = get_init(m, d1, d2, beta = beta_init, gamma = gamma_init, sigma2 = s2,
	phi2 = phi2_init, tau2 = tau2_init)

# ----- Self-tuned VWS within Gibbs (Version 1) -----
vws_ctrl = get_vws_control(tol1 = tol1, tol2 = tol2, max_rejects = 1e6,
	method = "vws-tune", N = 50)
control = get_gibbs_control(R = 3000, burn = 1000, thin = 1, report = 100,
	vws = vws_ctrl, save_latent = seq_len(m))
fixed = get_fixed()
vwg_out = gibbs(y, s2, X, Z, df, init, control, fixed)
print(vwg_out)

ess_vwg_sigma2 = ess(vwg_out$sigma2_hist)
ess_vwg_theta = ess(vwg_out$theta_hist)
quantile(ess_vwg_sigma2, probs)
quantile(ess_vwg_theta, probs)

par_vwg_mcmc = cbind(vwg_out$beta_hist, vwg_out$gamma_hist,
	vwg_out$phi2_hist, vwg_out$tau2_hist)
multiESS(vwg_out$beta_hist)
multiESS(vwg_out$gamma_hist)
ess(vwg_out$phi2_hist)
ess(vwg_out$tau2_hist)
multiESS(par_vwg_mcmc)

# ----- Metropolis within Gibbs -----
vws_ctrl = get_vws_control(method = "imh")
control = get_gibbs_control(R = 30000, burn = 28000, thin = 1, report = 1000,
	vws = vws_ctrl, save_latent = seq_len(m))
fixed = get_fixed()
mwg_out = gibbs(y, s2, X, Z, df, init, control, fixed)
print(mwg_out)

## Convert any NaN values of ESS (i.e., no chain movement) to zero
ess_mwg_sigma2 = ess(mwg_out$sigma2_hist)
ess_mwg_sigma2[is.na(ess_mwg_sigma2)] = 0
ess_mwg_theta = ess(mwg_out$theta_hist)
quantile(ess_mwg_sigma2, probs)
quantile(ess_mwg_theta, probs)

par_mwg_mcmc = cbind(mwg_out$beta_hist, mwg_out$gamma_hist,
	mwg_out$phi2_hist, mwg_out$tau2_hist)
multiESS(mwg_out$beta_hist)
multiESS(mwg_out$gamma_hist)
ess(mwg_out$phi2_hist)
ess(mwg_out$tau2_hist)
multiESS(par_mwg_mcmc)

# ----- Fit Fay-Herriot with Gibbs sampler -----
ctrl_fh = get_vws_control(method = "imh")
control_fh = get_gibbs_control(R = 3000, burn = 1000, thin = 1,
	report = 1000, vws = vws_ctrl, save_latent = seq_len(m))
fixed_fh = get_fixed(gamma = TRUE, tau2 = TRUE, sigma2 = TRUE)

fh_out = gibbs(y, s2, X, Z, df, init, control, fixed)
print(fh_out)

# ----- Create some plots from the results -----

# Dot plot of joint sampling variances versus estimated
g = data.frame(s2 = s2, joint = apply(vwg_out$sigma2_hist,2,mean)) %>%
	ggplot() +
	geom_point(aes(s2, joint)) +
	geom_abline(slope = 1, lty = 2, col = "red") +
	xlab(expression(s[i]^2)) +
	ylab(expression(hat(sigma)[i]^2)) +
	theme_light()
ggsave("variance-model-vs-estimated.pdf", g, width = 4, height = 4, unit="in")

# Model uncertainty in sigma2 versus area sample size
g = data.frame(df = dat_saipe$df, log_n = log(dat_saipe$hu_sampled), s2 = s2,
		lo = apply(vwg_out$sigma2_hist, 2, quantile, probs = alpha/2),
		hi = apply(vwg_out$sigma2_hist, 2, quantile, probs = 1-alpha/2)) %>%
	mutate(width = hi - lo) %>%
	ggplot() +
	geom_point(aes(x = log_n, y = width)) +
	xlab("Log of Area Sample Size") +
	ylab("Interval Width") +
	theme_light()
ggsave("variance-ci-width.pdf", g, width = 4, height = 4)

# Plot model uncertainty in theta from VWG versus FH. Ratio of CI widths versus
# area sample size.
g = data.frame(
		log_n = log(dat_saipe$hu_sampled),
		vwg_lo = apply(vwg_out$theta_hist, 2, quantile, probs = alpha/2),
		vwg_hi = apply(vwg_out$theta_hist, 2, quantile, probs = 1-alpha/2),
		fh_lo = apply(fh_out$theta_hist, 2, quantile, probs = alpha/2),
		fh_hi = apply(fh_out$theta_hist, 2, quantile, probs = 1-alpha/2)) %>%
	mutate(vwg_width = vwg_hi - vwg_lo) %>%
	mutate(fh_width = fh_hi - fh_lo) %>%
	mutate(ratio = fh_width / vwg_width) %>%
	ggplot() +
	geom_hex(aes(x = log_n, y = ratio), col = "black", lwd = 0.1, bins = 20) +
	geom_abline(intercept = 1, slope = 0, lty = 2, col = "red") +
	scale_fill_viridis_c() +
	xlab("Log of Area Sample Size") +
	ylab("Ratio") +
	theme_light() +
	scale_y_continuous(n.breaks = 15) +
	theme(legend.position = c(0.90, 0.20))
ggsave("theta-ci-ratio.pdf", g, width = 5, height = 5)

# Overlay histograms of the ESS for sigma2
g = data.frame(mwg = ess_mwg_sigma2, vwg = ess_vwg_sigma2) %>%
	ggplot() +
	geom_histogram(aes(x = mwg), col = "black", fill = "white", bins = 30, alpha = 0.4) +
	geom_histogram(aes(x = vwg),  col = "black", fill = "red2", bins = 30, alpha = 0.4) +
	xlab("ESS") +
	ylab("Count") +
	theme_light()
ggsave("sigma2-ess-hist.pdf", g, width = 5, height = 3)

# Overlay trace plot for the VWS and MWG sigma2
# Pick the 3 counties with worst MWG chains for sigma2
# and 3 counties with worst VWS chains
lowest_mwg = order(ess_mwg_sigma2)[1:3]
lowest_vws = order(ess_vwg_sigma2)[1:3]

plot_ess = c(lowest_mwg, lowest_vws)

for (ii in 1:length(plot_ess)) {
	idx = plot_ess[ii]
	g = data.frame(
			mwg = mwg_out$sigma2_hist[,idx],
			vws = vwg_out$sigma2_hist[,idx]) %>%
		mutate(x = row_number()) %>%
		ggplot() +
		geom_line(aes(x=x, y=vws), color = "red2", alpha = 0.4) +
		geom_line(aes(x=x, y=mwg), color = "blue", linewidth = 0.5, alpha = 1) +
		xlab("") +
		ylab(bquote(sigma[.(idx)]^2)) +
		theme_minimal()
	ff = sprintf("trace-%d.pdf", ii)
	ggsave(ff, g, width = 3, height = 2)
}

# Plot of estimates of sigma_i^2 for VWG vs MWG with interval widths

sigma2_mwg = apply(mwg_out$sigma2_hist, 2, mean)
sigma2_vwg = apply(vwg_out$sigma2_hist, 2, mean)
sigma2_sd_mwg = apply(mwg_out$sigma2_hist, 2, sd)
sigma2_sd_vwg = apply(vwg_out$sigma2_hist, 2, sd)
sigma2_ci_mwg = apply(mwg_out$sigma2_hist, 2, quantile, probs = c(alpha/2, 1 - alpha/2)) %>% t()
sigma2_ci_vwg = apply(vwg_out$sigma2_hist, 2, quantile, probs = c(alpha/2, 1 - alpha/2)) %>% t()
sigma2_width_mwg = apply(sigma2_ci_mwg, 1, diff)
sigma2_width_vwg = apply(sigma2_ci_vwg, 1, diff)

## Compare sigma2 between MWG and VWG using scatter/hex plots

g = data.frame(mwg = sigma2_mwg, vwg = sigma2_vwg) %>%
	add_column(log_n = log(dat_saipe$hu_sampled)) %>%
	mutate(ratio = mwg / vwg) %>%
	add_column(ess_mwg = ess_mwg_sigma2) %>%
	ggplot() +
	geom_point(aes(ess_mwg, ratio)) +
	geom_hline(yintercept = 1, lty = 2, col = "red") +
	xlab('ESS under MWG') +
	ylab('Ratio of Estimates') +
	scale_fill_viridis_c() +
	theme_light()
ggsave("sigma2-est-mwg-vs-vwg.pdf", g, width = 3.5, height = 3.5, unit="in")

g = data.frame(mwg = sigma2_width_mwg, vwg = sigma2_width_vwg) %>%
	add_column(log_n = log(dat_saipe$hu_sampled)) %>%
	add_column(ess_mwg = ess_mwg_sigma2) %>%
	mutate(ratio = mwg / vwg) %>%
	ggplot() +
	geom_point(aes(ess_mwg, ratio)) +
	geom_hline(yintercept = 1, lty = 2, col = "red") +
	xlab('ESS under MWG') +
	ylab('Ratio of Interval Widths') +
	scale_fill_viridis_c() +
	# scale_y_continuous(limits = c(0, 2)) +
	theme_light()
ggsave("sigma2-width-mwg-vs-vwg.pdf", g, width = 3.5, height = 3.5, unit="in")

g = data.frame(updates = vwg_out$sigma2_knot_updates_hist) %>%
	mutate(iter = row_number()) %>%
	ggplot() +
	geom_line(aes(iter, log10(updates + 1))) +
	xlab(NULL) +
	ylab("Count of Knot Updates (Log10)") +
	scale_x_continuous(n.breaks = 9) +
	scale_y_continuous(n.breaks = 10, expand = expansion()) +
	theme_light()
ggsave("knot-updates.pdf", g, width = 5, height = 3)

# Summaries of the regression parameters
xtable(summary(mwg_out), digits=3)
xtable(summary(vwg_out), digits=3)

if (FALSE) {
	# ----- VWS within Gibbs, no self-tuning -----
	# This takes a while to run, so do it last
	vws_ctrl = get_vws_control(tol1 = tol1, tol2 = tol2, max_rejects = 1e6,
		method = "vws-basic", N = 50)
	control = get_gibbs_control(R = 3000, burn = 1000, thin = 1, report = 1,
		vws = vws_ctrl, save_latent = seq_len(m))
	fixed = get_fixed()
	vwg_basic_out = gibbs(y, s2, X, Z, df, init, control, fixed)
	print(vwg_basic_out)

	ess_vwg_basic_sigma2 = ess(vwg_basic_out$sigma2_hist)
	ess_vwg_basic_theta = ess(vwg_basic_out$theta_hist)
	quantile(ess_vwg_basic_sigma2, probs)
	quantile(ess_vwg_basic_theta, probs)
}

save.image("results.Rdata")
