library(vws)
library(statmod)
library(tidyverse)

source("functions.R")

set.seed(1234)

R = 10000
n = 20
max_rejects = 1e6

mu = 0
lambda = 1
tol1_levels = c(0.75, 0.50)
tol2_levels = c(0.01, 0.001)
kappa_levels = c(10, 50)
tau_levels = c(0.5, 1.0)

tbl = expand.grid(kappa = kappa_levels, tau = tau_levels, tol1 = tol1_levels,
		tol2 = tol2_levels)
lb_list = list()
knots_list = list()
elapsed_list = list()
rejections_list = list()

# ----- Run the study -----
S = nrow(tbl)
for (s in seq_len(S)) {
	logger("Setting %d of %d\n", s, S)

	kappa = tbl$kappa[s]
	tau = tbl$tau[s]
	tol1 = tbl$tol1[s]
	tol2 = tbl$tol2[s]

	res_lb = matrix(NA, R, n)
	res_knots = matrix(NA, R, n)
	res_elapsed = numeric(R)
	res_rejections = numeric(R)

	# Transform back to original parameterization

	for (r in 1:R) {
		st = Sys.time()
		out = r_target(n, mu, tau, kappa, lambda, tol1, tol2, max_rejects)
		et = Sys.time()
		res_lb[r,] = out$log_bounds
		res_knots[r,] = out$knots
		res_elapsed[r] = as.numeric(et - st, units = "secs")
		res_rejections[r] = sum(out$rejections)
	}

	lb_med = apply(res_lb, 2, quantile, probs = 0.5)
	knots_med = apply(res_knots, 2, quantile, probs = 0.5)

	lb_list[[s]] = lb_med
	knots_list[[s]] = knots_med
	elapsed_list[[s]] = sum(res_elapsed)
	rejections_list[[s]] = sum(res_rejections)
}

# ----- Make 2^2 x 2^2 crosstabs -----

# Elapsed times
tbl |>
	add_column(elapsed = unlist(elapsed_list)) |>
	mutate(par = sprintf("kappa=%g, tau=%g", kappa, tau)) |>
	mutate(tol = sprintf("tol1=%g, tol2=%g", tol1, tol2)) |>
	mutate(par = as.factor(par)) |>
	mutate(tol = as.factor(tol)) |>
	mutate(elapsed = round(elapsed, 3)) |>
	xtabs(elapsed ~ par + tol, data = _)

# Rejection counts
tbl |>
	add_column(rejections = unlist(rejections_list)) |>
	mutate(par = sprintf("kappa=%g, tau=%g", kappa, tau)) |>
	mutate(tol = sprintf("tol1=%g, tol2=%g", tol1, tol2)) |>
	mutate(par = as.factor(par)) |>
	mutate(tol = as.factor(tol)) |>
	xtabs(rejections ~ par + tol, data = _)


# ----- Make plots of bounds and knots -----
# Plots are grouped by (tol1, tol2) values. Series within each plot vary with
# (kappa, tau).
for (idx1 in seq_along(tol1_levels)) {
for (idx2 in seq_along(tol2_levels)) {
	tol1 = tol1_levels[idx1]
	tol2 = tol2_levels[idx2]

	g1 = ggplot() +
		geom_hline(yintercept = log(tol1), lty = 2, col = "blue") +
		xlab("Iteration") +
		ylab("Log of Bound") +
		scale_y_continuous(n.breaks = 6) +
		theme_light()

	g2 = ggplot() +
		xlab("Iteration") +
		ylab("Number of Knots") +
		scale_y_continuous(n.breaks = 6) +
		theme_light()

	ltype = 0
	for (idx3 in seq_along(kappa_levels)) {
	for (idx4 in seq_along(tau_levels)) {
		kappa = kappa_levels[idx3]
		tau = tau_levels[idx4]
		idx_row = which(
			tbl$tol1 == tol1 &
			tbl$tol2 == tol2 &
			tbl$kappa == kappa &
			tbl$tau == tau)
		ltype = ltype + 1

		df = data.frame(iter = seq_len(n), x = lb_list[[idx_row]])
		g1 = g1 +
			geom_line(data = df, aes(iter, x)) +
			geom_point(data = df |> filter(row_number() %% 3 == 0),
				aes(iter, x), pch = ltype)

		df = data.frame(iter = seq_len(n), x = knots_list[[idx_row]])
		g2 = g2 +
			geom_line(data = df, aes(iter, x)) +
			geom_point(data = df |> filter(row_number() %% 3 == 0),
				aes(iter, x), pch = ltype)
	}
	}

	sprintf("bound-%d-%d.pdf", idx1, idx2) |>
		ggsave(plot = g1, width = 3, height = 2)
	sprintf("knots-%d-%d.pdf", idx1, idx2) |>
		ggsave(plot = g2, width = 3, height = 2)
}
}
