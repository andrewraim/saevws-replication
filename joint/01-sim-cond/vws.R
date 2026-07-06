library(vws)
library(statmod)
library(tidyverse)

source("../shared/functions.R", chdir = TRUE)

set.seed(1234)

R = 10000
n = 20
max_rejects = 1e6

mu = 0
lambda = 1
tol_suff_levels = c(0.75, 0.50)
tol_merge_levels = c(0.01, 0.001)
kappa_levels = c(10, 50)
tau_levels = c(0.5, 1.0)

tbl = expand.grid(kappa = kappa_levels, tau = tau_levels,
	tol_suff = tol_suff_levels, tol_merge = tol_merge_levels)
lb_list = list()
regions_list = list()
elapsed_list = list()
rejections_list = list()

# ----- Run the study -----
S = nrow(tbl)
for (s in seq_len(S))
{
	kappa = tbl$kappa[s]
	tau = tbl$tau[s]
	tol_suff = tbl$tol_suff[s]
	tol_merge = tbl$tol_merge[s]

	logger("[%d of %d] kappa %d  tau %0.1f  tol_suff %0.2f  tol_merge %0.3f\n",
		s, S, kappa, tau, tol_suff, tol_merge)

	res_lb = matrix(NA, R, n)
	res_regions = matrix(NA, R, n)
	res_elapsed = numeric(R)
	res_rejects = numeric(R)

	# Transform back to original parameterization

	for (r in 1:R) {
		st = Sys.time()
		out = r_target(n, mu, tau, kappa, lambda, tol_suff, tol_merge, max_rejects)
		et = Sys.time()
		res_lb[r,] = out$log_bounds
		res_regions[r,] = out$regions
		res_elapsed[r] = as.numeric(et - st, units = "secs")
		res_rejects[r] = sum(out$rejects)
	}

	lb_med = apply(res_lb, 2, quantile, probs = 0.5)
	regions_med = apply(res_regions, 2, quantile, probs = 0.5)

	lb_list[[s]] = lb_med
	regions_list[[s]] = regions_med
	elapsed_list[[s]] = sum(res_elapsed)
	rejections_list[[s]] = sum(res_rejects)
}

# ----- Make 2^2 x 2^2 crosstabs -----

# Elapsed times
tbl |>
	add_column(elapsed = unlist(elapsed_list)) |>
	mutate(par = sprintf("kappa=%g, tau=%g", kappa, tau)) |>
	mutate(tol = sprintf("tol_suff=%g, tol_merge=%g", tol_suff, tol_merge)) |>
	mutate(par = as.factor(par)) |>
	mutate(tol = as.factor(tol)) |>
	mutate(elapsed = round(elapsed, 3)) |>
	xtabs(elapsed ~ par + tol, data = _)

# Rejection counts
tbl |>
	add_column(rejections = unlist(rejections_list)) |>
	mutate(par = sprintf("kappa=%g, tau=%g", kappa, tau)) |>
	mutate(tol = sprintf("tol_suff=%g, tol_merge=%g", tol_suff, tol_merge)) |>
	mutate(par = as.factor(par)) |>
	mutate(tol = as.factor(tol)) |>
	xtabs(rejections ~ par + tol, data = _)


# ----- Make plots of bounds and knots -----
# Plots are grouped by (tol_suff, tol_merge) values. Series within each plot vary with
# (kappa, tau).
for (idx1 in seq_along(tol_suff_levels)) {
for (idx2 in seq_along(tol_merge_levels)) {
	tol_suff = tol_suff_levels[idx1]
	tol_merge = tol_merge_levels[idx2]

	g1 = ggplot() +
		geom_hline(yintercept = log(tol_suff), lty = 2, col = "blue") +
		xlab("Iteration") +
		ylab("Log of Bound") +
		scale_y_continuous(n.breaks = 6) +
		theme_light()

	g2 = ggplot() +
		xlab("Iteration") +
		ylab("Number of Regions") +
		scale_y_continuous(n.breaks = 6) +
		theme_light()

	ltype = 0
	for (idx3 in seq_along(kappa_levels)) {
	for (idx4 in seq_along(tau_levels)) {
		kappa = kappa_levels[idx3]
		tau = tau_levels[idx4]
		idx_row = which(
			tbl$tol_suff == tol_suff &
			tbl$tol_merge == tol_merge &
			tbl$kappa == kappa &
			tbl$tau == tau)
		ltype = ltype + 1

		df = data.frame(iter = seq_len(n), x = lb_list[[idx_row]])
		g1 = g1 +
			geom_line(data = df, aes(iter, x)) +
			geom_point(data = df |> filter(row_number() %% 3 == 0),
				aes(iter, x), pch = ltype)

		df = data.frame(iter = seq_len(n), x = regions_list[[idx_row]])
		g2 = g2 +
			geom_line(data = df, aes(iter, x)) +
			geom_point(data = df |> filter(row_number() %% 3 == 0),
				aes(iter, x), pch = ltype)

		# Print this to double check which series are which in the results
		printf("kappa %0.0f  tau %0.1f  tol1 %0.2f  tol2 %0.3f  ltype %d  maxregions %d\n",
			kappa, tau, tol_suff, tol_merge, ltype, max(df$x))

	}
	}

	sprintf("bound-%d-%d.pdf", idx1, idx2) |>
		ggsave(plot = g1, width = 3, height = 2)
	sprintf("regions-%d-%d.pdf", idx1, idx2) |>
		ggsave(plot = g2, width = 3, height = 2)
}
}
