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

# Args for VWG
tol_suff = 0.85
#tol_merge = 0.0001
tol_merge = 0.01

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

fixed = fixed_joint()

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

# ----- Independent Metropolis within Gibbs -----
inner_ctrl = control_inner(method = "imh")
control = control_joint(R = 30000, burn = 28000, thin = 1, report = 1000,
	inner = inner_ctrl, save_latent = seq_len(m))
mwg_out = gibbs_joint(y, s2, X, Z, df, init, control, fixed)
print(mwg_out)

z_geweke = geweke(mwg_out$sigma2_hist)
plot(density(z_geweke))
curve(dnorm, add = TRUE, lty = 2)
sum(abs(z_geweke) > 4)
sum(z_geweke < -3.5)

# Sort areas by rejection count and plot them against some respective data
hist(mwg_out$sigma2_rejects_areas)
idx = order(mwg_out$sigma2_rejects_areas)
plot(sort(mwg_out$sigma2_rejects_areas), y[idx])
plot(sort(mwg_out$sigma2_rejects_areas), df[idx])
plot(sort(mwg_out$sigma2_rejects_areas), s2[idx])

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

# ----- Adaptive Metropolis within Gibbs -----
inner_ctrl = control_inner(method = "am", am_varprop_init = 1)
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 1000,
	inner = inner_ctrl, save_latent = seq_len(m))
am_out = gibbs_joint(y, s2, X, Z, df, init, control, fixed)
print(am_out)

plot(am_out$beta_hist[,1], type = "l")
plot(am_out$gamma_hist[,1], type = "l")
plot(am_out$phi2_hist, type = "l")
plot(am_out$tau2_hist, type = "l")

ess_am_sigma2 = ess(am_out$sigma2_hist)
ess_am_theta = ess(am_out$theta_hist)

i = which.min(ess_am_sigma2)
plot(am_out$sigma2_hist[,i], type = "l")

hist(ess_am_sigma2)

z_geweke = geweke(am_out$sigma2_hist)
plot(density(z_geweke))
curve(dnorm, add = TRUE, lty = 2)
sum(abs(z_geweke) > 4)
sum(z_geweke < -3.5)

# ----- ARMS within Gibbs -----
inner_ctrl = control_inner(method = "arms")
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
	inner = inner_ctrl, save_latent = seq_len(m))
arms_out = gibbs_joint(y, s2, X, Z, df, init, control, fixed)
print(arms_out)

plot(arms_out$beta_hist[,1], type = "l")
plot(arms_out$gamma_hist[,1], type = "l")
plot(arms_out$phi2_hist, type = "l")
plot(arms_out$tau2_hist, type = "l")

ess_arms_sigma2 = ess(arms_out$sigma2_hist)
ess_arms_theta = ess(arms_out$theta_hist)
quantile(ess_arms_sigma2, probs)
quantile(ess_arms_theta, probs)

i = which.min(ess_arms_sigma2)
plot(arms_out$sigma2_hist[,i], type = "l")
hist(ess_arms_sigma2)

z_geweke = geweke(arms_out$sigma2_hist)
plot(density(z_geweke))
curve(dnorm, add = TRUE, lty = 2)
sum(abs(z_geweke) > 4)
sum(z_geweke < -3.5)

# ----- Self-tuned VWS within Gibbs Version 1 -----
inner_ctrl = control_inner(tol_suff = tol_suff, tol_merge = tol_merge,
	max_rejects = 1e6, method = "vws-tune", N = 50)
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
	inner = inner_ctrl, save_latent = seq_len(m))
vwg_out = gibbs_joint(y, s2, X, Z, df, init, control, fixed)
print(vwg_out)

plot(vwg_out$beta_hist[,1], type = "l")
plot(vwg_out$gamma_hist[,1], type = "l")
plot(vwg_out$phi2_hist, type = "l")
plot(vwg_out$tau2_hist, type = "l")

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

i = which.min(ess_vwg_sigma2)
plot(vwg_out$sigma2_hist[,i], type = "l")

z_geweke = geweke(vwg_out$sigma2_hist)
plot(density(z_geweke))
curve(dnorm, add = TRUE, lty = 2)
sum(abs(z_geweke) > 4)
sum(z_geweke < -3.5)

# ----- Self-tuned VWS within Gibbs Version 2  -----
# Stop tuning after an initial period
inner_ctrl = control_inner(tol_suff = tol_suff, tol_merge = tol_merge,
	max_rejects = 1e6, method = "vws-tune", N = 50, tune = 100)
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
	inner = inner_ctrl, save_latent = seq_len(m))
vwg2_out = gibbs_joint(y, s2, X, Z, df, init, control, fixed)
print(vwg2_out)

ess_vwg2_sigma2 = ess(vwg2_out$sigma2_hist)

i = which.min(ess_vwg2_sigma2)
plot(vwg2_out$sigma2_hist[,i], type = "l")

# ----- Self-tuned VWS within Gibbs Version 3  -----
# Stop tuning after an initial period, then use proposal with MH algorithm
# instead of rejection sampling. It may be interesting that this does not work
# as well: it doesn't run much faster than version 2 and some of the sigma2
# chains aren't mixing that well.
tol_suff2 = 0.25
tol_merge2 = 0.001
inner_ctrl = control_inner(tol_suff = tol_suff2, tol_merge = tol_merge2,
	max_rejects = 1e6, method = "mh-vws", N = 50, tune = 400)
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
	inner = inner_ctrl, save_latent = seq_len(m))
vwg3_out = gibbs_joint(y, s2, X, Z, df, init, control, fixed)
print(vwg3_out)

ess_vwg3_sigma2 = ess(vwg3_out$sigma2_hist)

i = which.min(ess_vwg3_sigma2)
plot(vwg3_out$sigma2_hist[,i], type = "l")

# ----- Fit Fay-Herriot with Gibbs sampler -----
ctrl_fh = control_inner(method = "imh")
control_fh = control_joint(R = 3000, burn = 1000, thin = 1,
	report = 1000, inner = inner_ctrl, save_latent = seq_len(m))
fixed_fh = fixed_joint(gamma = TRUE, tau2 = TRUE, sigma2 = TRUE)

fh_out = gibbs_joint(y, s2, X, Z, df, init, control, fixed)
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

# Plot number of tuned VWS proposals by iteration
data.frame(tuned = vwg_out$sigma2_tuned_hist) %>%
	mutate(iter = row_number()) %>%
	filter(iter > 100) %>%
	ggplot() +
	geom_line(aes(iter, tuned)) +
	geom_rect(xmin = 0, xmax = 100,  ymin = 0,  ymax = Inf, fill = "red") +
	scale_y_continuous(breaks = seq(0, 95, by = 5), minor_breaks = NULL) +
	xlab("Iteration") +
	ylab("Number of Tuned VWS Proposals") +
	theme_minimal()

# Plot number of VWS rejections per area by iteration
data.frame(tuned = vwg_out$sigma2_rejects_hist / m) %>%
	mutate(iter = row_number()) %>%
	filter(iter > 100) %>%
	ggplot() +
	geom_line(aes(iter, tuned)) +
	geom_rect(xmin = 0, xmax = 20, ymin = min(vwg_out$sigma2_rejects_hist / m),
		ymax = Inf, fill = "red") +
	xlab("Iteration") +
	ylab("Number of VWS Rejections Per Area") +
	theme_minimal()

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

g = data.frame(updates = vwg_out$sigma2_tunes_hist) %>%
	mutate(iter = row_number()) %>%
	ggplot() +
	geom_line(aes(iter, log10(updates + 1))) +
	xlab(NULL) +
	ylab("Count of Updates (Log10)") +
	scale_x_continuous(n.breaks = 9) +
	scale_y_continuous(n.breaks = 10, expand = expansion()) +
	theme_light()
ggsave("region-updates-log10.pdf", g, width = 5, height = 3)

g = data.frame(updates = vwg_out$sigma2_tunes_hist) %>%
	mutate(iter = row_number()) %>%
	filter(iter > 24) %>%
	ggplot() +
	geom_line(aes(iter, updates)) +
	xlab(NULL) +
	ylab("Number of Region Updates") +
	scale_x_continuous(n.breaks = 9) +
	scale_y_continuous(n.breaks = 10, expand = expansion()) +
	theme_light()
ggsave("region-updates.pdf", g, width = 5, height = 3)

g = data.frame(count = vwg_out$sigma2_comps_hist) %>%
	mutate(iter = row_number()) %>%
	filter(iter > 24) %>%
	ggplot() +
	geom_line(aes(iter, count)) +
	xlab(NULL) +
	ylab("Number of Regions") +
	scale_x_continuous(n.breaks = 9) +
	scale_y_continuous(n.breaks = 10, expand = expansion()) +
	theme_light()
ggsave("region-counts.pdf", g, width = 5, height = 3)

g = data.frame(count = vwg_out$sigma2_rejects_hist) %>%
	mutate(iter = row_number()) %>%
	filter(iter > 24) %>%
	ggplot() +
	geom_line(aes(iter, count)) +
	xlab(NULL) +
	ylab("Number of Rejections") +
	scale_x_continuous(n.breaks = 9) +
	scale_y_continuous(n.breaks = 10, expand = expansion()) +
	theme_light()
ggsave("rejection-counts.pdf", g, width = 5, height = 3)


# Summaries of the regression parameters
xtable(summary(mwg_out), digits=3)
xtable(summary(vwg_out), digits=3)

if (FALSE) {
	# ----- VWS within Gibbs, no self-tuning -----
	# This takes a while to run, so do it last
	inner_ctrl = control_inner(tol_suff = tol_suff, tol_merge = tol_merge,
		max_rejects = 1e6, method = "vws-basic", N = 50)
	control = control_joint(R = 3000, burn = 1000, thin = 1, report = 1,
		inner = inner_ctrl, save_latent = seq_len(m))
	fixed = fixed_joint()
	vwg_basic_out = gibbs_joint(y, s2, X, Z, df, init, control, fixed)
	print(vwg_basic_out)

	ess_vwg_basic_sigma2 = ess(vwg_basic_out$sigma2_hist)
	ess_vwg_basic_theta = ess(vwg_basic_out$theta_hist)
	quantile(ess_vwg_basic_sigma2, probs)
	quantile(ess_vwg_basic_theta, probs)
}

save.image("results.Rdata")
