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

# ----- Independent Metropolis within Gibbs -----
inner_ctrl = control_inner(method = "imh")
control = control_joint(R = 30000, burn = 28000, thin = 1, report = 1000,
	inner = inner_ctrl, save_latent = seq_len(m))
imh_out = gibbs_joint(y, s2, X, Z, df, init, control)
print(imh_out)

z_geweke = geweke(imh_out$sigma2_hist)
plot(density(z_geweke))
curve(dnorm, add = TRUE, lty = 2)
sum(abs(z_geweke) > 4)
sum(z_geweke < -3.5)

# Sort areas by rejection count and plot them against some respective data
hist(imh_out$sigma2_rejects_areas)
idx = order(imh_out$sigma2_rejects_areas)
plot(sort(imh_out$sigma2_rejects_areas), y[idx])
plot(sort(imh_out$sigma2_rejects_areas), df[idx])
plot(sort(imh_out$sigma2_rejects_areas), s2[idx])

## Convert any NaN values of ESS (i.e., no chain movement) to zero
ess_imh_sigma2 = ess(imh_out$sigma2_hist)
ess_imh_sigma2[is.na(ess_imh_sigma2)] = 0
ess_imh_theta = ess(imh_out$theta_hist)
quantile(ess_imh_sigma2, probs)
quantile(ess_imh_theta, probs)

par_imh_mcmc = cbind(imh_out$beta_hist, imh_out$gamma_hist,
	imh_out$phi2_hist, imh_out$tau2_hist)
multiESS(imh_out$beta_hist)
multiESS(imh_out$gamma_hist)
ess(imh_out$phi2_hist)
ess(imh_out$tau2_hist)
multiESS(par_imh_mcmc)

tbl_ess = tibble(
	method = "IMH",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_imh_sigma2, probs[1]),
	ess2 = quantile(ess_imh_sigma2, probs[2]),
	ess3 = quantile(ess_imh_sigma2, probs[3]),
	elapsed = sum(unlist(imh_out$elapsed)),
	rejections = sum(imh_out$sigma2_rejects_hist)
)

# ----- Adaptive Metropolis within Gibbs -----
inner_ctrl = control_inner(method = "amh", am_varprop_init = 1)
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 1000,
	inner = inner_ctrl, save_latent = seq_len(m))
amh_out = gibbs_joint(y, s2, X, Z, df, init, control)
print(amh_out)

plot(amh_out$beta_hist[,1], type = "l")
plot(amh_out$gamma_hist[,1], type = "l")
plot(amh_out$phi2_hist, type = "l")
plot(amh_out$tau2_hist, type = "l")

ess_amh_sigma2 = ess(amh_out$sigma2_hist)
ess_amh_theta = ess(amh_out$theta_hist)

i = which.min(ess_amh_sigma2)
plot(amh_out$sigma2_hist[,i], type = "l")

hist(ess_amh_sigma2)

z_geweke = geweke(amh_out$sigma2_hist)
plot(density(z_geweke))
curve(dnorm, add = TRUE, lty = 2)
sum(abs(z_geweke) > 4)
sum(z_geweke < -3.5)

tbl_ess = tbl_ess %>% add_row(
	method = "AMH",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_amh_sigma2, probs[1]),
	ess2 = quantile(ess_amh_sigma2, probs[2]),
	ess3 = quantile(ess_amh_sigma2, probs[3]),
	elapsed = sum(unlist(amh_out$elapsed)),
	rejections = sum(amh_out$sigma2_rejects_hist)
)

# ----- ARMS within Gibbs -----
inner_ctrl = control_inner(method = "arms")
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
	inner = inner_ctrl, save_latent = seq_len(m))
arms_out = gibbs_joint(y, s2, X, Z, df, init, control)
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

tbl_ess = tbl_ess %>% add_row(
	method = "ARMS",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_arms_sigma2, probs[1]),
	ess2 = quantile(ess_arms_sigma2, probs[2]),
	ess3 = quantile(ess_arms_sigma2, probs[3]),
	elapsed = sum(unlist(arms_out$elapsed)),
	rejections = sum(arms_out$sigma2_rejects_hist)
)

# ----- VWS1 within Gibbs -----
vws1_out = list()

for (l in seq_along(nrow(tol_levels)))
{
	tol_suff = tol_levels$tol_suff[l]
	tol_merge = tol_levels$tol_merge[l]

	inner_ctrl = control_inner(tol_suff = tol_suff, tol_merge = tol_merge,
		max_rejects = 1e6, method = "vws-tune", N = 50)
	control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
		inner = inner_ctrl, save_latent = seq_len(m))
	gibbs_out = gibbs_joint(y, s2, X, Z, df, init, control)
	print(gibbs_out)

	plot(gibbs_out$beta_hist[,1], type = "l")
	plot(gibbs_out$gamma_hist[,1], type = "l")
	plot(gibbs_out$phi2_hist, type = "l")
	plot(gibbs_out$tau2_hist, type = "l")

	ess_sigma2 = ess(gibbs_out$sigma2_hist)
	ess_theta = ess(gibbs_out$theta_hist)
	quantile(ess_sigma2, probs)
	quantile(ess_theta, probs)

	# par_vws_mcmc = cbind(gibbs_out$beta_hist, gibbs_out$gamma_hist,
	#	gibbs_out$phi2_hist, gibbs_out$tau2_hist)
	# multiESS(gibbs_out$beta_hist)
	# multiESS(gibbs_out$gamma_hist)
	# ess(gibbs_out$phi2_hist)
	# ess(gibbs_out$tau2_hist)
	# multiESS(par_vws_mcmc)

	i = which.min(ess_vws_sigma2)
	plot(gibbs_out$sigma2_hist[,i], type = "l")

	z_geweke = geweke(gibbs_out$sigma2_hist)
	plot(density(z_geweke))
	curve(dnorm, add = TRUE, lty = 2)
	sum(abs(z_geweke) > 4)
	sum(z_geweke < -3.5)

	tbl_ess = tbl_ess %>% add_row(
		method = "VWS1",
		tol_suff = tol_suff,
		tol_merge = tol_merge,
		ess1 = quantile(ess_sigma2, probs[1]),
		ess2 = quantile(ess_sigma2, probs[2]),
		ess3 = quantile(ess_sigma2, probs[3]),
		elapsed = sum(unlist(gibbs_out$elapsed)),
		rejections = sum(gibbs_out$sigma2_rejects_hist)
	)

	vws1_out[[l]] = gibbs_out
}

# ----- VWS2 within Gibbs  -----
# Stop tuning after an initial period (100 Gibbs iterations)

vws2_out = list()

for (l in seq_along(nrow(tol_levels)))
{
	tol_suff = tol_levels$tol_suff[l]
	tol_merge = tol_levels$tol_merge[l]

	inner_ctrl = control_inner(tol_suff = tol_suff, tol_merge = tol_merge,
		max_rejects = 1e6, method = "vws-tune", N = 50, tune = 100)
	control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
		inner = inner_ctrl, save_latent = seq_len(m))
	gibbs_out = gibbs_joint(y, s2, X, Z, df, init, control)
	print(gibbs_out)

	ess_sigma2 = ess(gibbs_out$sigma2_hist)

	i = which.min(ess_sigma2)
	plot(gibbs_out$sigma2_hist[,i], type = "l")

	tbl_ess = tbl_ess %>% add_row(
		method = "VWS2",
		tol_suff = tol_suff,
		tol_merge = tol_merge,
		ess1 = quantile(ess_sigma2, probs[1]),
		ess2 = quantile(ess_sigma2, probs[2]),
		ess3 = quantile(ess_sigma2, probs[3]),
		elapsed = sum(unlist(gibbs_out$elapsed)),
		rejections = sum(gibbs_out$sigma2_rejects_hist)
	)

	vws2_out[[l]] = gibbs_out
}

# ----- VWS3 within Gibbs -----
# Stop tuning after an initial period, then use proposal with MH algorithm
# instead of rejection sampling. It may be interesting that this does not work
# as well: it doesn't run much faster than version 2 and some of the sigma2
# chains aren't mixing that well.

tol_suff = 0.25
tol_merge = 0.001

inner_ctrl = control_inner(tol_suff = tol_suff, tol_merge = tol_merge,
	max_rejects = 1e6, method = "mh-vws", N = 50, tune = 400)
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
	inner = inner_ctrl, save_latent = seq_len(m))
gibbs_out = gibbs_joint(y, s2, X, Z, df, init, control)
print(gibbs_out)

ess_sigma2 = ess(gibbs_out$sigma2_hist)

i = which.min(ess_sigma2)
plot(gibbs_out$sigma2_hist[,i], type = "l")

# ----- Fit Fay-Herriot with Gibbs sampler -----
ctrl_fh = control_inner(method = "imh")
control_fh = control_joint(R = 3000, burn = 1000, thin = 1,
	report = 1000, inner = inner_ctrl, save_latent = seq_len(m))
fixed_fh = fixed_joint(gamma = TRUE, tau2 = TRUE, sigma2 = TRUE)

fh_out = gibbs_joint(y, s2, X, Z, df, init, control, fixed_fh)
print(fh_out)

# ----- VWS0 within Gibbs -----
# Construct a new proposal whenever a sigma2 step is encountered.
# This takes a while to run, so do it last, and only if switch is enabled.

if (run_vws0)
{
	inner_ctrl = control_inner(tol_suff = 0.85, tol_merge = 0.001,
		max_rejects = 1e6, method = "vws-basic", N = 50)
	control = control_joint(R = 3000, burn = 1000, thin = 1, report = 1,
		inner = inner_ctrl, save_latent = seq_len(m))
	vws0_out = gibbs_joint(y, s2, X, Z, df, init, control)
	print(vws0_out)

	ess_sigma2 = ess(vws0_out$sigma2_hist)
	ess_theta = ess(vws0_out$theta_hist)
	quantile(ess_sigma2, probs)
	quantile(ess_theta, probs)

	tbl_ess = tbl_ess %>% add_row(
		method = "VWS2",
		tol_suff = tol_suff,
		tol_merge = tol_merge,
		ess1 = quantile(ess_sigma2, probs[1]),
		ess2 = quantile(ess_sigma2, probs[2]),
		ess3 = quantile(ess_sigma2, probs[3]),
		elapsed = sum(unlist(vws0_out$elapsed)),
		rejections = sum(vws0_out$sigma2_rejects_hist)
	)

	vws2_out[[l]] = vws0_out
}


# ----- Create some plots from the results -----

# Dot plot of joint sampling variances versus estimated
g = data.frame(s2 = s2, joint = apply(vws_out$sigma2_hist,2,mean)) %>%
	ggplot() +
	geom_point(aes(s2, joint)) +
	geom_abline(slope = 1, lty = 2, col = "red") +
	xlab(expression(s[i]^2)) +
	ylab(expression(hat(sigma)[i]^2)) +
	theme_light()
ggsave("variance-model-vs-estimated.pdf", g, width = 4, height = 4, unit="in")

# Model uncertainty in sigma2 versus area sample size
g = data.frame(df = dat_saipe$df, log_n = log(dat_saipe$hu_sampled), s2 = s2,
		lo = apply(vws_out$sigma2_hist, 2, quantile, probs = alpha/2),
		hi = apply(vws_out$sigma2_hist, 2, quantile, probs = 1-alpha/2)) %>%
	mutate(width = hi - lo) %>%
	ggplot() +
	geom_point(aes(x = log_n, y = width)) +
	xlab("Log of Area Sample Size") +
	ylab("Interval Width") +
	theme_light()
ggsave("variance-ci-width.pdf", g, width = 4, height = 4)

# Plot model uncertainty in theta from VWS versus FH. Ratio of CI widths versus
# area sample size.
g = data.frame(
		log_n = log(dat_saipe$hu_sampled),
		vws_lo = apply(vws_out$theta_hist, 2, quantile, probs = alpha/2),
		vws_hi = apply(vws_out$theta_hist, 2, quantile, probs = 1-alpha/2),
		fh_lo = apply(fh_out$theta_hist, 2, quantile, probs = alpha/2),
		fh_hi = apply(fh_out$theta_hist, 2, quantile, probs = 1-alpha/2)) %>%
	mutate(vws_width = vws_hi - vws_lo) %>%
	mutate(fh_width = fh_hi - fh_lo) %>%
	mutate(ratio = fh_width / vws_width) %>%
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
g = data.frame(imh = ess_imh_sigma2, vws = ess_vws_sigma2) %>%
	ggplot() +
	geom_histogram(aes(x = imh), col = "black", fill = "white", bins = 30, alpha = 0.4) +
	geom_histogram(aes(x = vws),  col = "black", fill = "red2", bins = 30, alpha = 0.4) +
	xlab("ESS") +
	ylab("Count") +
	theme_light()
ggsave("sigma2-ess-hist.pdf", g, width = 5, height = 3)

# Overlay trace plot for the VWS and IMH sigma2
# Pick the 3 counties with worst IMH chains for sigma2
# and 3 counties with worst VWS chains
lowest_imh = order(ess_imh_sigma2)[1:3]
lowest_vws = order(ess_vws_sigma2)[1:3]

plot_ess = c(lowest_imh, lowest_vws)

for (ii in 1:length(plot_ess)) {
	idx = plot_ess[ii]
	g = data.frame(
			imh = imh_out$sigma2_hist[,idx],
			vws = vws_out$sigma2_hist[,idx]) %>%
		mutate(x = row_number()) %>%
		ggplot() +
		geom_line(aes(x=x, y=vws), color = "red2", alpha = 0.4) +
		geom_line(aes(x=x, y=imh), color = "blue", linewidth = 0.5, alpha = 1) +
		xlab("") +
		ylab(bquote(sigma[.(idx)]^2)) +
		theme_minimal()
	ff = sprintf("trace-%d.pdf", ii)
	ggsave(ff, g, width = 3, height = 2)
}

# Plot of estimates of sigma_i^2 for VWS vs IMH with interval widths

sigma2_imh = apply(imh_out$sigma2_hist, 2, mean)
sigma2_vws = apply(vws_out$sigma2_hist, 2, mean)
sigma2_sd_imh = apply(imh_out$sigma2_hist, 2, sd)
sigma2_sd_vws = apply(vws_out$sigma2_hist, 2, sd)
sigma2_ci_imh = apply(imh_out$sigma2_hist, 2, quantile, probs = c(alpha/2, 1 - alpha/2)) %>% t()
sigma2_ci_vws = apply(vws_out$sigma2_hist, 2, quantile, probs = c(alpha/2, 1 - alpha/2)) %>% t()
sigma2_width_imh = apply(sigma2_ci_imh, 1, diff)
sigma2_width_vws = apply(sigma2_ci_vws, 1, diff)

# Plot number of tuned VWS proposals by iteration
data.frame(tuned = vws_out$sigma2_tuned_hist) %>%
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
data.frame(tuned = vws_out$sigma2_rejects_hist / m) %>%
	mutate(iter = row_number()) %>%
	filter(iter > 100) %>%
	ggplot() +
	geom_line(aes(iter, tuned)) +
	geom_rect(xmin = 0, xmax = 20, ymin = min(vws_out$sigma2_rejects_hist / m),
		ymax = Inf, fill = "red") +
	xlab("Iteration") +
	ylab("Number of VWS Rejections Per Area") +
	theme_minimal()

## Compare sigma2 between IMH and VWS using scatter/hex plots

g = data.frame(imh = sigma2_imh, vws = sigma2_vws) %>%
	add_column(log_n = log(dat_saipe$hu_sampled)) %>%
	mutate(ratio = imh / vws) %>%
	add_column(ess_imh = ess_imh_sigma2) %>%
	ggplot() +
	geom_point(aes(ess_imh, ratio)) +
	geom_hline(yintercept = 1, lty = 2, col = "red") +
	xlab('ESS under IMH') +
	ylab('Ratio of Estimates') +
	scale_fill_viridis_c() +
	theme_light()
ggsave("sigma2-est-imh-vs-vws.pdf", g, width = 3.5, height = 3.5, unit="in")

g = data.frame(imh = sigma2_width_imh, vws = sigma2_width_vws) %>%
	add_column(log_n = log(dat_saipe$hu_sampled)) %>%
	add_column(ess_imh = ess_imh_sigma2) %>%
	mutate(ratio = imh / vws) %>%
	ggplot() +
	geom_point(aes(ess_imh, ratio)) +
	geom_hline(yintercept = 1, lty = 2, col = "red") +
	xlab('ESS under IMH') +
	ylab('Ratio of Interval Widths') +
	scale_fill_viridis_c() +
	# scale_y_continuous(limits = c(0, 2)) +
	theme_light()
ggsave("sigma2-width-imh-vs-vws.pdf", g, width = 3.5, height = 3.5, unit="in")

g = data.frame(updates = vws_out$sigma2_tunes_hist) %>%
	mutate(iter = row_number()) %>%
	ggplot() +
	geom_line(aes(iter, log10(updates + 1))) +
	xlab(NULL) +
	ylab("Count of Updates (Log10)") +
	scale_x_continuous(n.breaks = 9) +
	scale_y_continuous(n.breaks = 10, expand = expansion()) +
	theme_light()
ggsave("region-updates-log10.pdf", g, width = 5, height = 3)

g = data.frame(updates = vws_out$sigma2_tunes_hist) %>%
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

g = data.frame(count = vws_out$sigma2_comps_hist) %>%
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

g = data.frame(count = vws_out$sigma2_rejects_hist) %>%
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
xtable(summary(imh_out), digits=3)
xtable(summary(vws_out), digits=3)

save.image("results.Rdata")
