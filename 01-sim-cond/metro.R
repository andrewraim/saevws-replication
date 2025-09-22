library(vws)
library(saevws)
library(statmod)
library(tidyverse)
library(knitr)
library(mcmcse)

source("functions.R")
Rcpp::sourceCpp("samplers.cpp")

set.seed(1234)

n = 200000

mu = 0
lambda = 1
kappa_levels = c(10, 50)
tau_levels = c(0.5, 1.0)

tbl = expand.grid(kappa = kappa_levels, tau = tau_levels) |>
	mutate(ess = NA) |>
	mutate(elapsed = NA) |>
	mutate(rejects = NA) |>
	mutate(autocor = NA)

# ----- Run the study -----
S = nrow(tbl)
for (s in seq_len(S)) {
	logger("Setting %d of %d\n", s, S)

	kappa = tbl$kappa[s]
	tau = tbl$tau[s]
	idx_kappa = which(kappa == kappa_levels)
	idx_tau = which(tau == tau_levels)

	# Transform back to original parameterization
	opt_out = optimize(d_target_unnorm, interval = c(0, 10), maximum = TRUE,
		mu = mu, lambda = lambda, tau = tau, kappa = kappa, log = TRUE)
	init = opt_out$maximum

	st = Sys.time()
	out = r_metro(n, init, mu, tau, kappa, lambda)
	et = Sys.time()

	tbl$ess[s] = ess(out$draws)
	tbl$elapsed[s] = as.numeric(et - st, units = "secs")
	tbl$rejects[s] = out$rejects
	tbl$autocor[s] = acf(out$draws, lag.max = 1, plot = F)$acf[2,1,1]

	g = data.frame(x = out$draws) |>
		mutate(iter = row_number()) |>
		filter(iter > 50000) |>
		ggplot() +
		geom_line(aes(x = iter, y = x)) +
		xlab("") +
		ylab("") +
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 25, hjust = 1))

	ff = sprintf("trace-%d-%d.png", idx_kappa, idx_tau)
	ggsave(ff, g, width = 2.5, height = 1.75)

	g = data.frame(x = out$draws) |>
		mutate(iter = row_number()) |>
		filter(iter > n - 1000) |>
		ggplot() +
		geom_line(aes(x = iter, y = x)) +
		xlab("") +
		ylab("") +
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 25, hjust = 1))
	ff = sprintf("trace-%d-%d-closeup.pdf", idx_kappa, idx_tau)
	ggsave(ff, g, width = 3, height = 2)
}

# ----- Make a table -----

tbl |>
	select(-elapsed) |>
	kable(format = "latex")
