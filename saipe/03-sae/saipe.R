library(tidyverse)
library(saevws)
library(mcmcse)
library(xtable)
library(coda)
library(knitr)

source("../shared/functions.R", chdir = TRUE)

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

run_vws0 = TRUE

# ----- Data setup -----
ff = file.path("..", "data", "saipe.csv")
saipe = read_csv(ff) %>%
	mutate(df = 0.36 * sqrt(hu_sampled)) %>%
	mutate(y = log(pov_count)) %>%
	mutate(s2 = pov_se^2 / pov_count^2) %>%
	filter(pov_count > 1 & df > 1) %>%
	filter(!is.na(snap)) # Removing one county with missing snap

m = nrow(saipe)
X = model.matrix(~ log1p(snap) + log1p(pep), data = saipe)
Z = model.matrix(~ log(hu_sampled) , data = saipe)
df = saipe$df
y = saipe$y
s2 = saipe$s2
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

# ----- IMH within Gibbs -----
inner_ctrl = control_inner(method = "imh")
control = control_joint(R = 30000, burn = 28000, thin = 1, report = 1000,
	inner = inner_ctrl, save_latent = seq_len(m))
imh_out = gibbs_joint(y, s2, X, Z, df, init, control)
print(imh_out)

# Sort areas by rejection count and plot them against some respective data
# hist(imh_out$sigma2_rejects_areas)
# idx = order(imh_out$sigma2_rejects_areas)
# plot(sort(imh_out$sigma2_rejects_areas), y[idx])
# plot(sort(imh_out$sigma2_rejects_areas), df[idx])
# plot(sort(imh_out$sigma2_rejects_areas), s2[idx])

## A few of the sigma2 chains are likely not to move. The ESS function produces
## an NaN for them. Convert them to zero.
ess_sigma2 = ess(imh_out$sigma2_hist)
ess_sigma2[is.na(ess_sigma2)] = 0
# ess_theta = ess(imh_out$theta_hist)
# quantile(ess_sigma2, probs)
# quantile(ess_theta, probs)

# par_mcmc = cbind(imh_out$beta_hist, imh_out$gamma_hist,
#	imh_out$phi2_hist, imh_out$tau2_hist)
# multiESS(imh_out$beta_hist)
# multiESS(imh_out$gamma_hist)
# ess(imh_out$phi2_hist)
# ess(imh_out$tau2_hist)
# multiESS(par_mcmc)

tbl_ess = tibble(
	method = "IMH",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_sigma2, probs[1]),
	ess2 = quantile(ess_sigma2, probs[2]),
	ess3 = quantile(ess_sigma2, probs[3]),
	elapsed = sum(unlist(imh_out$elapsed)),
	rejections = sum(imh_out$sigma2_rejects_hist)
)

# ----- AMH within Gibbs -----
inner_ctrl = control_inner(method = "amh", amh_varprop_init = 1)
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 1000,
	inner = inner_ctrl, save_latent = seq_len(m))
amh_out = gibbs_joint(y, s2, X, Z, df, init, control)
print(amh_out)

# plot(amh_out$beta_hist[,1], type = "l")
# plot(amh_out$gamma_hist[,1], type = "l")
# plot(amh_out$phi2_hist, type = "l")
# plot(amh_out$tau2_hist, type = "l")

ess_sigma2 = ess(amh_out$sigma2_hist)
# ess_theta = ess(amh_out$theta_hist)

# i = which.min(ess_sigma2)
# plot(amh_out$sigma2_hist[,i], type = "l")

# hist(ess_sigma2)

tbl_ess = tbl_ess %>% add_row(
	method = "AMH",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_sigma2, probs[1]),
	ess2 = quantile(ess_sigma2, probs[2]),
	ess3 = quantile(ess_sigma2, probs[3]),
	elapsed = sum(unlist(amh_out$elapsed)),
	rejections = sum(amh_out$sigma2_rejects_hist)
)

# ----- ARMS within Gibbs -----
inner_ctrl = control_inner(method = "arms")
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
	inner = inner_ctrl, save_latent = seq_len(m))
arms_out = gibbs_joint(y, s2, X, Z, df, init, control)
print(arms_out)

# plot(arms_out$beta_hist[,1], type = "l")
# plot(arms_out$gamma_hist[,1], type = "l")
# plot(arms_out$phi2_hist, type = "l")
# plot(arms_out$tau2_hist, type = "l")

ess_sigma2 = ess(arms_out$sigma2_hist)
# ess_theta = ess(arms_out$theta_hist)
# quantile(ess_sigma2, probs)
# quantile(ess_theta, probs)

# i = which.min(ess_sigma2)
# plot(arms_out$sigma2_hist[,i], type = "l")
# hist(ess_sigma2)

tbl_ess = tbl_ess %>% add_row(
	method = "ARMS",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_sigma2, probs[1]),
	ess2 = quantile(ess_sigma2, probs[2]),
	ess3 = quantile(ess_sigma2, probs[3]),
	elapsed = sum(unlist(arms_out$elapsed)),
	rejections = sum(arms_out$sigma2_rejects_hist)
)

# ----- VWS0 within Gibbs -----
# Construct a new proposal whenever a sigma2 step is encountered.
inner_ctrl = control_inner(tol_suff = 0.85, tol_merge = 0.001,
	max_rejects = 1e6, method = "vws-basic", N = 50)
control = control_joint(R = 3000, burn = 1000, thin = 1, report = 50,
	inner = inner_ctrl, save_latent = seq_len(m))
vws0_out = gibbs_joint(y, s2, X, Z, df, init, control)
print(vws0_out)

ess_sigma2 = ess(vws0_out$sigma2_hist)
# ess_theta = ess(vws0_out$theta_hist)
# quantile(ess_sigma2, probs)
# quantile(ess_theta, probs)

g = plot_rejects(vws0_out$sigma2_rejects_hist, burn = 0, tol = 1.0)
ggsave("rejects-vws0.pdf", g, width = 3, height = 2)

tbl_ess = tbl_ess %>% add_row(
	method = "VWS0",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_sigma2, probs[1]),
	ess2 = quantile(ess_sigma2, probs[2]),
	ess3 = quantile(ess_sigma2, probs[3]),
	elapsed = sum(unlist(vws0_out$elapsed)),
	rejections = sum(vws0_out$sigma2_rejects_hist)
)

# ----- VWS1 within Gibbs -----
vws1_out = list()

for (l in seq_len(nrow(tol_levels)))
{
	tol_suff = tol_levels$tol_suff[l]
	tol_merge = tol_levels$tol_merge[l]

	inner_ctrl = control_inner(tol_suff = tol_suff, tol_merge = tol_merge,
		max_rejects = 1e6, method = "vws-tune", N = 50)
	control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
		inner = inner_ctrl, save_latent = seq_len(m))
	gibbs_out = gibbs_joint(y, s2, X, Z, df, init, control)
	print(gibbs_out)

	# plot(gibbs_out$beta_hist[,1], type = "l")
	# plot(gibbs_out$gamma_hist[,1], type = "l")
	# plot(gibbs_out$phi2_hist, type = "l")
	# plot(gibbs_out$tau2_hist, type = "l")

	ess_sigma2 = ess(gibbs_out$sigma2_hist)
	# ess_theta = ess(gibbs_out$theta_hist)
	# quantile(ess_sigma2, probs)
	# quantile(ess_theta, probs)

	# par_vws_mcmc = cbind(gibbs_out$beta_hist, gibbs_out$gamma_hist,
	#	gibbs_out$phi2_hist, gibbs_out$tau2_hist)
	# multiESS(gibbs_out$beta_hist)
	# multiESS(gibbs_out$gamma_hist)
	# ess(gibbs_out$phi2_hist)
	# ess(gibbs_out$tau2_hist)
	# multiESS(par_vws_mcmc)

	# i = which.min(ess_sigma2)
	# plot(gibbs_out$sigma2_hist[,i], type = "l")

	g = plot_tunes(gibbs_out$sigma2_tunes_hist, burn = 500, tol = 0.05)
	ff = sprintf("tunes-vws1-%d.pdf", l)
	ggsave(ff, g, width = 3, height = 2)

	g = plot_comps(gibbs_out$sigma2_comps_hist, burn = 500, tol = 0.05)
	ff = sprintf("comps-vws1-%d.pdf", l)
	ggsave(ff, g, width = 3, height = 2)

	g = plot_rejects(gibbs_out$sigma2_rejects_hist, burn = 500, tol = 0.05)
	ff = sprintf("rejects-vws1-%d.pdf", l)
	ggsave(ff, g, width = 3, height = 2)

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

for (l in seq_len(nrow(tol_levels)))
{
	tol_suff = tol_levels$tol_suff[l]
	tol_merge = tol_levels$tol_merge[l]
	tune = 100

	inner_ctrl = control_inner(tol_suff = tol_suff, tol_merge = tol_merge,
		max_rejects = 1e6, method = "vws-tune", N = 50, tune = tune)
	control = control_joint(R = 3000, burn = 1000, thin = 1, report = 100,
		inner = inner_ctrl, save_latent = seq_len(m))
	gibbs_out = gibbs_joint(y, s2, X, Z, df, init, control)
	print(gibbs_out)

	ess_sigma2 = ess(gibbs_out$sigma2_hist)

	# i = which.min(ess_sigma2)
	# plot(gibbs_out$sigma2_hist[,i], type = "l")

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
# instead of rejection sampling. Seems interesting that this does not work
# as well: it doesn't run much faster than version 2 and some of the sigma2
# chains aren't mixing that well.

if (FALSE) {
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
}

# ----- Fit Fay-Herriot with Gibbs sampler -----
# This is just for reference, not for comparison to other sampling methods.
inner_ctrl = control_inner(method = "imh")
control = control_joint(R = 3000, burn = 1000, thin = 1,
	report = 1000, inner = inner_ctrl, save_latent = seq_len(m))
fixed_fh = fixed_joint(gamma = TRUE, tau2 = TRUE, sigma2 = TRUE)
fh_out = gibbs_joint(y, s2, X, Z, df, init, control, fixed_fh)
print(fh_out)

save.image("results.Rdata")

# ----- Additional Plots -----

# Dot plot of joint sampling variances versus estimated
g = data.frame(s2 = s2, joint = apply(vws2_out[[4]]$sigma2_hist, 2, mean)) %>%
	ggplot() +
	geom_point(aes(s2, joint)) +
	geom_abline(slope = 1, lty = 2, col = "red") +
	xlab(expression(s[i]^2)) +
	ylab(expression(hat(sigma)[i]^2)) +
	theme_light()
ggsave("variance-model-vs-estimated.pdf", g, width = 4, height = 4, unit="in")

# Model uncertainty in sigma2 versus area sample size
g = data.frame(df = saipe$df, log_n = log(saipe$hu_sampled), s2 = s2,
		lo = apply(vws2_out[[4]]$sigma2_hist, 2, quantile, probs = alpha/2),
		hi = apply(vws2_out[[4]]$sigma2_hist, 2, quantile, probs = 1-alpha/2)) %>%
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
		log_n = log(saipe$hu_sampled),
		vws_lo = apply(vws2_out[[4]]$theta_hist, 2, quantile, probs = alpha/2),
		vws_hi = apply(vws2_out[[4]]$theta_hist, 2, quantile, probs = 1-alpha/2),
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
# g = data.frame(
#		imh = ess(imh_out$sigma2_hist),
#		vws = ess(vws2_out[[4]]$sigma2_hist)) %>%
#	ggplot() +
#	geom_histogram(aes(x = imh), col = "black", fill = "white", bins = 30, alpha = 0.4) +
#	geom_histogram(aes(x = vws),  col = "black", fill = "red2", bins = 30, alpha = 0.4) +
#	xlab("ESS") +
#	ylab("Count") +
#	theme_light()
# ggsave("sigma2-ess-hist.pdf", g, width = 5, height = 3)

df_plot = data.frame(
		imh = ess(imh_out$sigma2_hist),
		arms = ess(arms_out$sigma2_hist),
		amh = ess(amh_out$sigma2_hist),
		vws = ess(vws2_out[[4]]$sigma2_hist)) %>%
	mutate(iter = row_number())
# df_quantiles = data.frame(probs = c(0.20, 0.40, 0.60, 0.80)) %>%
#	mutate(imh = quantile(df_plot$imh, probs, na.rm = TRUE)) %>%
#	mutate(arms = quantile(df_plot$arms, probs, na.rm = TRUE)) %>%
#	mutate(amh = quantile(df_plot$amh, probs, na.rm = TRUE)) %>%
#	mutate(vws = quantile(df_plot$vws, probs, na.rm = TRUE))
df_points = data.frame(x = c(500, 1000, 1500)) %>%
	mutate(imh = ecdf(df_plot$imh)(x)) %>%
	mutate(arms = ecdf(df_plot$arms)(x)) %>%
	mutate(amh = ecdf(df_plot$amh)(x)) %>%
	mutate(vws = ecdf(df_plot$vws)(x))
g = pivot_longer(df_plot, cols = c("imh", "arms", "amh", "vws")) %>%
	ggplot() +
	stat_ecdf(aes(x = value, group = name)) +
	geom_point(data = df_points, aes(x, imh), pch = 16, cex = 3) +
	geom_point(data = df_points, aes(x, amh), pch = 17, cex = 3) +
	geom_point(data = df_points, aes(x, arms), pch = 18, cex = 3) +
	geom_point(data = df_points, aes(x, vws), pch = 15, cex = 3) +
	scale_x_continuous(expand = expansion(0,0)) +
	xlab("ESS") +
	ylab("Empirical CDF") +
	theme_light()
ggsave("sigma2-ess-ecdf.pdf", g, width = 4, height = 3)

# Overlay trace plot for the VWS and IMH sigma2
# Pick the 3 counties with worst IMH chains for sigma2
# and 3 counties with worst VWS chains
ess_imh = ess(imh_out$sigma2_hist)
ess_imh[is.na(ess_imh)] = 0
ess_vws = ess(vws2_out[[4]]$sigma2_hist)
lowest_imh = order(ess_imh)[1:3]
lowest_vws = order(ess_vws)[1:3]
plot_ess = c(lowest_imh, lowest_vws)

print(saipe[lowest_imh,])
print(saipe[lowest_vws,])

for (ii in 1:length(plot_ess)) {
	idx = plot_ess[ii]

	g = data.frame(
			imh = imh_out$sigma2_hist[,idx],
			vws = vws2_out[[4]]$sigma2_hist[,idx]) %>%
		mutate(x = row_number()) %>%
		ggplot() +
		geom_line(aes(x=x, y=vws), color = "red2", alpha = 0.4) +
		geom_line(aes(x=x, y=imh), color = "blue", linewidth = 0.5, alpha = 1) +
		xlab("") +
		ylab(bquote(sigma[.(idx)]^2)) +
		theme_light()
	ff = sprintf("imh-trace-%d.pdf", ii)
	ggsave(ff, g, width = 3, height = 2)
}

# Also look at the three worst mixing chains under AMH
ess_amh = ess(amh_out$sigma2_hist)
lowest_amh = order(ess_amh)[1:3]
for (ii in 1:length(lowest_amh)) {
	idx = lowest_amh[ii]
	g = data.frame(amh = amh_out$sigma2_hist[,idx]) %>%
		mutate(x = row_number()) %>%
		ggplot() +
		geom_line(aes(x=x, y=amh), linewidth = 0.5) +
		xlab("") +
		ylab(bquote(sigma[.(idx)]^2)) +
		theme_light()
	ff = sprintf("amh-trace-%d.pdf", ii)
	ggsave(ff, g, width = 3, height = 2)
}
print(saipe[lowest_amh,])

# Also look at the three worst mixing chains under ARMS
ess_arms = ess(arms_out$sigma2_hist)
lowest_arms = order(ess_arms)[1:3]
for (ii in 1:length(lowest_arms)) {
	idx = lowest_arms[ii]
	g = data.frame(arms = arms_out$sigma2_hist[,idx]) %>%
		mutate(x = row_number()) %>%
		ggplot() +
		geom_line(aes(x=x, y=arms), linewidth = 0.5) +
		xlab("") +
		ylab(bquote(sigma[.(idx)]^2)) +
		theme_light()
	ff = sprintf("arms-trace-%d.pdf", ii)
	ggsave(ff, g, width = 3, height = 2)
}
print(saipe[lowest_arms,])

# Plot of estimates of sigma_i^2 for VWS vs IMH with interval widths

sigma2_imh = apply(imh_out$sigma2_hist, 2, mean)
sigma2_vws = apply(vws2_out[[4]]$sigma2_hist, 2, mean)
sigma2_sd_imh = apply(imh_out$sigma2_hist, 2, sd)
sigma2_sd_vws = apply(vws2_out[[4]]$sigma2_hist, 2, sd)
sigma2_ci_imh = apply(imh_out$sigma2_hist, 2, quantile, probs = c(alpha/2, 1 - alpha/2)) %>% t()
sigma2_ci_vws = apply(vws2_out[[4]]$sigma2_hist, 2, quantile, probs = c(alpha/2, 1 - alpha/2)) %>% t()
sigma2_width_imh = apply(sigma2_ci_imh, 1, diff)
sigma2_width_vws = apply(sigma2_ci_vws, 1, diff)
ess_imh_sigma2 = ess(imh_out$sigma2_hist)

# Plot number of tuned VWS proposals by iteration
# data.frame(tuned = vws_out$sigma2_tuned_hist) %>%
# 	mutate(iter = row_number()) %>%
# 	filter(iter > 100) %>%
# 	ggplot() +
# 	geom_line(aes(iter, tuned)) +
# 	geom_rect(xmin = 0, xmax = 100,  ymin = 0,  ymax = Inf, fill = "red") +
# 	scale_y_continuous(breaks = seq(0, 95, by = 5), minor_breaks = NULL) +
# 	xlab("Iteration") +
# 	ylab("Number of Tuned VWS Proposals") +
# 	theme_minimal()

# Plot number of VWS rejections per area by iteration
# data.frame(tuned = vws_out$sigma2_rejects_hist / m) %>%
# 	mutate(iter = row_number()) %>%
# 	filter(iter > 100) %>%
# 	ggplot() +
# 	geom_line(aes(iter, tuned)) +
# 	geom_rect(xmin = 0, xmax = 20, ymin = min(vws_out$sigma2_rejects_hist / m),
# 		ymax = Inf, fill = "red") +
# 	xlab("Iteration") +
# 	ylab("Number of VWS Rejections Per Area") +
# 	theme_minimal()

## Compare sigma2 between IMH and VWS using scatter/hex plots

g = data.frame(imh = sigma2_imh, vws = sigma2_vws) %>%
	add_column(log_n = log(saipe$hu_sampled)) %>%
	mutate(ratio = imh / vws) %>%
	add_column(ess_imh = ess_imh_sigma2) %>%
	ggplot() +
	geom_point(aes(ess_imh, ratio)) +
	geom_hline(yintercept = 1, lty = 2, col = "red") +
	xlab('ESS under IMH') +
	ylab('Ratio of Estimates') +
	scale_fill_viridis_c() +
	theme_light()
ggsave("sigma2-est-imh-vs-vws.pdf", g, width = 4, height = 3, unit="in")

g = data.frame(imh = sigma2_width_imh, vws = sigma2_width_vws) %>%
	add_column(log_n = log(saipe$hu_sampled)) %>%
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
ggsave("sigma2-width-imh-vs-vws.pdf", g, width = 4, height = 3, unit="in")

# g = data.frame(updates = vws_out$sigma2_tunes_hist) %>%
#	mutate(iter = row_number()) %>%
#	ggplot() +
#	geom_line(aes(iter, log10(updates + 1))) +
#	xlab(NULL) +
#	ylab("Count of Updates (Log10)") +
#	scale_x_continuous(n.breaks = 9) +
#	scale_y_continuous(n.breaks = 10, expand = expansion()) +
#	theme_light()
# ggsave("region-updates-log10.pdf", g, width = 5, height = 3)

# Plot rejections for VWS with limited tuning period. Omit iterations during
# the tuning period.
for (l in 1:nrow(tol_levels)) {
	g = vws2_out[[l]]$sigma2_rejects_hist %>%
		tail(-tune) %>%
		plot_rejects(burn = 0, tol = 100)
	ff = sprintf("rejects-vws2-%d-tail.pdf", l)
	ggsave(ff, g, width = 3, height = 2)
}

# ----- Tables to summarize MCMC results -----
tbl_ess %>%
	mutate(tol_merge = format(tol_merge, scientific = FALSE)) %>%
	mutate(ess1_sec = format(ess1 / elapsed, digits = 2, big.mark = ",")) %>%
	mutate(ess2_sec = format(ess2 / elapsed, digits = 2, big.mark = ",")) %>%
	mutate(ess3_sec = format(ess3 / elapsed, digits = 2, big.mark = ",")) %>%
	mutate(ess1 = format(ess1, digits = 2, big.mark = ",")) %>%
	mutate(ess2 = format(ess2, digits = 2, big.mark = ",")) %>%
	mutate(ess3 = format(ess3, digits = 2, big.mark = ",")) %>%
	mutate(rejections = format(rejections, digits = 2, big.mark = ",", scientific = FALSE)) %>%
	mutate(elapsed = sprintf("%0.2f", elapsed)) %>%
	select(method, tol_suff, tol_merge, rejections, elapsed, ess1, ess2, ess3,
		ess1_sec, ess2_sec, ess3_sec) %>%
	kable(format = "latex", linesep = "")

# Summaries of the regression parameters
xtable(summary(imh_out), digits = 4)
xtable(summary(vws2_out[[4]]), digits = 4)

s_imh = summary(imh_out)
s_amh = summary(arms_out)
s_arms = summary(arms_out)
s_vws0 = summary(vws0_out)
s_vws1_1 = summary(vws1_out[[1]])
s_vws1_2 = summary(vws1_out[[2]])
s_vws1_3 = summary(vws1_out[[3]])
s_vws1_4 = summary(vws1_out[[4]])
s_vws2_1 = summary(vws2_out[[1]])
s_vws2_2 = summary(vws2_out[[2]])
s_vws2_3 = summary(vws2_out[[3]])
s_vws2_4 = summary(vws2_out[[4]])

res_theta_ess = rbind(
	s_imh$ess,
	s_amh$ess,
	s_arms$ess,
	s_vws0$ess,
	s_vws1_1$ess,
	s_vws1_2$ess,
	s_vws1_3$ess,
	s_vws1_4$ess,
	s_vws2_1$ess,
	s_vws2_2$ess,
	s_vws2_3$ess,
	s_vws2_4$ess
) %>% round() %>% format(big.mark = ",")
colnames(res_theta_ess) = rownames(s_imh)
method = c("IMH", "AMH", "ARMS", "VWS0", rep("VWS1", 4), rep("VWS2", 4))
tol_suff = c(rep(NA, 4), tol_levels$tol_suff, tol_levels$tol_suff)
tol_merge = c(rep(NA, 4), tol_levels$tol_merge, tol_levels$tol_merge)
tbl_theta_ess = res_theta_ess %>%
	as.data.frame() %>%
	add_column(method) %>%
	add_column(tol_suff) %>%
	add_column(tol_merge) %>%
	select(method, tol_suff, tol_merge, everything())
kable(tbl_theta_ess, format = "latex", linesep = "")

# ----- Experimental: Gelman-Rubin Diagnostic -----
# Run three additional chains with IMH and then diagnose the four together.
imh2_out = gibbs_joint(y, s2, X, Z, df, init, control)
imh3_out = gibbs_joint(y, s2, X, Z, df, init, control)
imh4_out = gibbs_joint(y, s2, X, Z, df, init, control)

# TBD: can we get the statistic for all m sigma2 entries, or is that too much?
ess_sigma2 = ess(imh_out$sigma2_hist)
lowest_imh = order(ess_sigma2)[1:3]
# lowest_imh = 1:m

# Try Gelman-Rubin with coda package on IMH.
mcmc_list = mcmc.list(
	as.mcmc(imh_out$sigma2_hist[,lowest_imh]),
	as.mcmc(imh2_out$sigma2_hist[,lowest_imh]),
	as.mcmc(imh3_out$sigma2_hist[,lowest_imh]),
	as.mcmc(imh4_out$sigma2_hist[,lowest_imh])
)
gr_imh = gelman.diag(mcmc_list, confidence = 0.95, autoburnin = FALSE)

# Try Gelman-Rubin with coda package on VWS. Base it on the four runs of VWS1.
# It shouldn't matter that they are based on different tunings of the rejection
# sampler.

mcmc_list = mcmc.list(
	as.mcmc(vws2_out[[1]]$sigma2_hist[,lowest_imh]),
	as.mcmc(vws2_out[[2]]$sigma2_hist[,lowest_imh]),
	as.mcmc(vws2_out[[3]]$sigma2_hist[,lowest_imh]),
	as.mcmc(vws2_out[[4]]$sigma2_hist[,lowest_imh])
)
gr_vws = gelman.diag(mcmc_list, confidence = 0.95, autoburnin = FALSE)


# ----- Experimental: Geweke Diagnostic -----
# TBD: Geweke diagnostic seems about to detect the worst mixing chains. But how
# to summarize it in a table for results?

z_geweke = geweke(imh_out$sigma2_hist)
pval = pnorm(2 * abs(z_geweke), lower.tail = FALSE)
idx = order(pval)[1:6]
pval[idx]
plot(imh_out$sigma2_hist[,idx[6]], type = "l")

plot(density(z_geweke))
curve(dnorm, add = TRUE, lty = 2)
sum(abs(z_geweke) > 4)
sum(z_geweke < -3.5)

