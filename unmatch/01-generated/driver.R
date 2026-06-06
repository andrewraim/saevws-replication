library(tidyverse)
library(saevws)
library(mcmcse)
library(xtable)
library(coda)
library(knitr)

source("../../saipe/shared/functions.R", chdir = TRUE)

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

# ----- Read data -----
unmatch = read_csv("../data/unmatch.csv")

m = nrow(unmatch)
X = model.matrix(~ x, data = unmatch)
sigma = unmatch$sigma
y = unmatch$y

# ----- IMH within Gibbs -----
init = init_unmatch(m, d = ncol(X))
inner = control_inner(method = "imh")
control = control_unmatch(R = 30000, burn = 28000, thin = 1, report = 5000,
	save_latent = 1:m, inner = inner)
imh_out = gibbs_unmatch(y, sigma, X, init, control)
print(imh_out)

ess_mu = ess(imh_out$mu_hist)
ess_mu[is.na(ess_mu)] = 0

# i = which.min(ess_mu)
# plot(imh_out$mu_hist[,i], type = "l")
# abline(h = mu_true[i], lty = 2, col = "red")

tbl_ess = tibble(
	method = "IMH",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_mu, probs[1]),
	ess2 = quantile(ess_mu, probs[2]),
	ess3 = quantile(ess_mu, probs[3]),
	elapsed = sum(unlist(imh_out$elapsed)),
	rejections = sum(imh_out$mu_rejects_hist)
)

# ----- AMH within Gibbs -----
init = init_unmatch(m, d = ncol(X))
inner = control_inner(method = "am", am_varprop_init = 25, am_varprop_eps = 1e-4)
control = control_unmatch(R = 10000, burn = 8000, thin = 1, report = 5000,
	save_latent = 1:m, inner = inner)
amh_out = gibbs_unmatch(y, sigma, X, init, control)
print(amh_out)

# plot(amh_out$beta_hist[,1], type = "l")
# plot(amh_out$beta_hist[,2], type = "l")
# plot(amh_out$tau2_hist, type = "l")

ess_mu = ess(amh_out$mu_hist)
# hist(ess_mu)

# i = which.min(ess_mu)
# plot(amh_out$mu_hist[,i], type = "l")
# abline(h = mu_true[i], lty = 2, col = "red")

tbl_ess = tbl_ess %>% add_row(
	method = "AMH",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_mu, probs[1]),
	ess2 = quantile(ess_mu, probs[2]),
	ess3 = quantile(ess_mu, probs[3]),
	elapsed = sum(unlist(amh_out$elapsed)),
	rejections = sum(amh_out$mu_rejects_hist)
)

# ----- ARMS within Gibbs -----
init = init_unmatch(m, d = ncol(X))
inner = control_inner(method = "arms")
control = control_unmatch(R = 3000, burn = 1000, thin = 1, report = 100,
	save_latent = 1:m, inner = inner)
arms_out = gibbs_unmatch(y, sigma, X, init, control)
print(arms_out)

# plot(arms_out$beta_hist[,1], type = "l")
# plot(arms_out$beta_hist[,2], type = "l")
# plot(arms_out$tau2_hist, type = "l")

ess_mu = ess(arms_out$mu_hist)
# hist(ess_mu)

# i = which.min(ess_mu)
# plot(arms_out$mu_hist[,i], type = "l")
# abline(h = mu_true[i], lty = 2, col = "red")

tbl_ess = tbl_ess %>% add_row(
	method = "ARMS",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_mu, probs[1]),
	ess2 = quantile(ess_mu, probs[2]),
	ess3 = quantile(ess_mu, probs[3]),
	elapsed = sum(unlist(arms_out$elapsed)),
	rejections = sum(arms_out$mu_rejects_hist)
)

# ----- VWS0 within Gibbs -----
init = init_unmatch(m, d = ncol(X))
inner = control_inner(method = "vws-basic", tol_suff = 0.85)
control = control_unmatch(R = 3000, burn = 1000, thin = 1, report = 100,
	save_latent = 1:m, inner = inner)
vws0_out = gibbs_unmatch(y, sigma, X, init, control)
print(vws0_out)

# plot(vws0_out$beta_hist[,1], type = "l")
# plot(vws0_out$beta_hist[,2], type = "l")
# plot(vws0_out$tau2_hist, type = "l")

ess_mu = ess(vws0_out$mu_hist)
# hist(ess_mu)

g = plot_rejects(vws0_out$mu_rejects_hist, burn = 0, tol = 1.0)
ggsave("rejects-vws0.pdf", g, width = 3, height = 2)

# i = which.min(ess_mu)
# plot(vws0_out$mu_hist[,i], type = "l")
# abline(h = mu_true[i], lty = 2, col = "red")

tbl_ess = tbl_ess %>% add_row(
	method = "VWS0",
	tol_suff = NA,
	tol_merge = NA,
	ess1 = quantile(ess_mu, probs[1]),
	ess2 = quantile(ess_mu, probs[2]),
	ess3 = quantile(ess_mu, probs[3]),
	elapsed = sum(unlist(vws0_out$elapsed)),
	rejections = sum(vws0_out$mu_rejects_hist)
)

# ----- VWS1 within Gibbs -----
vws1_out = list()

for (l in seq_len(nrow(tol_levels)))
{
	tol_suff = tol_levels$tol_suff[l]
	tol_merge = tol_levels$tol_merge[l]

	init = init_unmatch(m, d = ncol(X))
	inner = control_inner(method = "vws-tune", tol_suff = tol_suff,
		tol_merge = tol_merge)
	control = control_unmatch(R = 3000, burn = 1000, thin = 1, report = 100,
		save_latent = 1:m, inner = inner)
	gibbs_out = gibbs_unmatch(y, sigma, X, init, control)
	print(gibbs_out)

	# plot(gibbs_out$beta_hist[,1], type = "l")
	# plot(gibbs_out$beta_hist[,2], type = "l")
	# plot(gibbs_out$tau2_hist, type = "l")

	ess_mu = ess(gibbs_out$mu_hist)
	# hist(ess_mu)

	# i = which.min(ess_mu)
	# plot(gibbs_out$mu_hist[,i], type = "l")
	# abline(h = mu_true[i], lty = 2, col = "red")

	g = plot_tunes(gibbs_out$mu_tunes_hist, burn = 500, tol = 0.05)
	ff = sprintf("tunes-vws1-%d.pdf", l)
	ggsave(ff, g, width = 3, height = 2)

	g = plot_comps(gibbs_out$mu_comps_hist, burn = 500, tol = 0.01)
	ff = sprintf("comps-vws1-%d.pdf", l)
	ggsave(ff, g, width = 3, height = 2)

	g = plot_rejects(gibbs_out$mu_rejects_hist, burn = 500, tol = 0.05)
	ff = sprintf("rejects-vws1-%d.pdf", l)
	ggsave(ff, g, width = 3, height = 2)

	tbl_ess = tbl_ess %>% add_row(
		method = "VWS1",
		tol_suff = tol_suff,
		tol_merge = tol_merge,
		ess1 = quantile(ess_mu, probs[1]),
		ess2 = quantile(ess_mu, probs[2]),
		ess3 = quantile(ess_mu, probs[3]),
		elapsed = sum(unlist(gibbs_out$elapsed)),
		rejections = sum(gibbs_out$mu_rejects_hist)
	)

	vws1_out[[l]] = gibbs_out
}

# ----- VWS2 within Gibbs -----
vws2_out = list()

for (l in seq_len(nrow(tol_levels)))
{
	tol_suff = tol_levels$tol_suff[l]
	tol_merge = tol_levels$tol_merge[l]

	init = init_unmatch(m, d = ncol(X))
	inner = control_inner(method = "vws-tune", tol_suff = tol_suff,
		tol_merge = tol_merge, tune = 100)
	control = control_unmatch(R = 3000, burn = 1000, thin = 1, report = 100,
		save_latent = 1:m, inner = inner)
	gibbs_out = gibbs_unmatch(y, sigma, X, init, control)
	print(gibbs_out)

	ess_mu = ess(gibbs_out$mu_hist)
	# hist(ess_mu)

	# i = which.min(ess_mu)
	# plot(gibbs_out$mu_hist[,i], type = "l")
	# abline(h = mu_true[i], lty = 2, col = "red")

	tbl_ess = tbl_ess %>% add_row(
		method = "VWS2",
		tol_suff = tol_suff,
		tol_merge = tol_merge,
		ess1 = quantile(ess_mu, probs[1]),
		ess2 = quantile(ess_mu, probs[2]),
		ess3 = quantile(ess_mu, probs[3]),
		elapsed = sum(unlist(gibbs_out$elapsed)),
		rejections = sum(gibbs_out$mu_rejects_hist)
	)

	vws2_out[[l]] = gibbs_out
}

save.image("results.Rdata")

# ----- Additional Plots -----
df_plot = data.frame(
		imh = ess(imh_out$mu_hist),
		arms = ess(arms_out$mu_hist),
		amh = ess(amh_out$mu_hist),
		vws = ess(vws1_out[[4]]$mu_hist)) %>%
	mutate(iter = row_number())
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
	scale_x_continuous(expand = expansion(0, 0)) +
	scale_y_continuous(expand = expansion(0, 0.025)) +
	xlab("ESS") +
	ylab("Empirical CDF") +
	theme_light()
ggsave("mu-ess-ecdf.pdf", g, width = 4, height = 3)

# Overlay trace plot for the VWS and IMH sigma2
# Pick the 3 counties with worst IMH chains for sigma2
# and 3 counties with worst VWS chains
ess_imh = ess(imh_out$mu_hist)
ess_imh[is.na(ess_imh)] = 0
ess_vws = ess(vws1_out[[4]]$mu_hist)
lowest_imh = order(ess_imh)[1:3]
lowest_vws = order(ess_vws)[1:3]
plot_ess = c(lowest_imh, lowest_vws)

for (ii in 1:length(plot_ess)) {
	idx = plot_ess[ii]

	g = data.frame(
			imh = imh_out$mu_hist[,idx],
			vws = vws1_out[[4]]$mu_hist[,idx]) %>%
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
ess_amh = ess(amh_out$mu_hist)
lowest_amh = order(ess_amh)[1:3]
for (ii in 1:length(lowest_amh)) {
	idx = lowest_amh[ii]
	g = data.frame(amh = amh_out$mu_hist[,idx]) %>%
		mutate(x = row_number()) %>%
		ggplot() +
		geom_line(aes(x=x, y=amh), linewidth = 0.5) +
		xlab("") +
		ylab(bquote(sigma[.(idx)]^2)) +
		theme_light()
	ff = sprintf("amh-trace-%d.pdf", ii)
	ggsave(ff, g, width = 3, height = 2)
}

# Also look at the three worst mixing chains under ARMS
ess_arms = ess(arms_out$mu_hist)
lowest_arms = order(ess_arms)[1:3]
for (ii in 1:length(lowest_arms)) {
	idx = lowest_arms[ii]
	g = data.frame(arms = arms_out$mu_hist[,idx]) %>%
		mutate(x = row_number()) %>%
		ggplot() +
		geom_line(aes(x=x, y=arms), linewidth = 0.5) +
		xlab("") +
		ylab(bquote(sigma[.(idx)]^2)) +
		theme_light()
	ff = sprintf("arms-trace-%d.pdf", ii)
	ggsave(ff, g, width = 3, height = 2)
}

# Plot rejections for VWS with limited tuning period. Omit iterations during
# the tuning period.
tune = 100
for (l in 1:nrow(tol_levels)) {
	g = vws2_out[[l]]$mu_rejects_hist %>%
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
xtable(summary(vws1_out[[4]]), digits = 4)

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

