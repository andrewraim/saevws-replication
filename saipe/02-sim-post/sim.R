library(vws)
library(saevws)
library(mcmcse)
library(tidyverse)
library(jsonlite)

set.seed(1234)

config = fromJSON("config.json")
m = config$m
tol1 = config$tol1
tol2 = config$tol2

# Set fixed arguments
N_sim = 500
mwg_draws = 30000
mwg_burn = 28000
mwg_thin = 1
mwg_report = 5000
vwg_draws = 3000
vwg_burn = 1000
vwg_thin = 1
vwg_report = 1000

# If the following is true, save individual fits from the three methods. This
# will use a lot of storage if the number of reps is large.
save_fits = FALSE

# Probabilities for ESS and ACF Quantiles
probs = c(0.000, 0.010, 0.025)

ff_mwg = sprintf("mwg.csv")
ff_vwg = sprintf("vwg.csv")

# ----- Check for arguments set by the launch script -----
stopifnot(exists("m"))
stopifnot(exists("tol1"))
stopifnot(exists("tol2"))

# ----- Set up design matrices and data-generating parameters -----
tau_true = sqrt(0.25)
phi_true =  sqrt(0.20)

beta_true = c(1.5, 0.85)
x = rnorm(m, mean = 8, sd = 2)
X = cbind(1, x)
Xbeta_true = X %*% beta_true

n = rchisq(m, df = 16)
df = n - 1
Z = cbind(1, log(n))
gamma_true = c(2.6, -1)
Zgamma_true = Z %*% gamma_true

d1 = ncol(X)
d2 = ncol(Z)

save.image("setup.Rdata")

# ----- Prepare the results table -----
# We are trying to make the code restartable, so that it can continue where it
# left off if it gets halted.

if (file.exists(ff_mwg)) {
	df_mwg = read_csv(ff_mwg, show_col_types = F)
	stopifnot(nrow(df_mwg) == N_sim)
	a = c("rep", "essQ1", "essQ2", "essQ3", "essNA", "par_mess", "theta_essQ1",
		"theta_essQ2", "theta_essQ3", "elapsed", "rejections")
	stopifnot(all(a %in% colnames(df_mwg)))
} else {
	df_mwg = data.frame(rep = seq_len(N_sim)) |>
		add_column(essQ1 = NA) |>
		add_column(essQ2 = NA) |>
		add_column(essQ3 = NA) |>
		add_column(acfQ1 = NA) |>
		add_column(acfQ2 = NA) |>
		add_column(acfQ3 = NA) |>
		add_column(essNA = NA) |>
		add_column(par_mess = NA) |>
		add_column(theta_essQ1 = NA) |>
		add_column(theta_essQ2 = NA) |>
		add_column(theta_essQ3 = NA) |>
		add_column(elapsed = NA) |>
		add_column(rejections = NA)
}

if (file.exists(ff_vwg)) {
	df_vwg = read_csv(ff_vwg, show_col_types = F)
	stopifnot(nrow(df_vwg) == N_sim)
	a = c("rep", "essQ1", "essQ2", "essQ3", "essNA", "par_mess", "theta_essQ1",
		"theta_essQ2", "theta_essQ3", "elapsed", "rejections",
		"knot_updates_burn", "knot_updates_keep")
	stopifnot(all(a %in% colnames(df_vwg)))
} else {
	df_vwg = data.frame(rep = seq_len(N_sim)) |>
		add_column(essQ1 = NA) |>
		add_column(essQ2 = NA) |>
		add_column(essQ3 = NA) |>
		add_column(acfQ1 = NA) |>
		add_column(acfQ2 = NA) |>
		add_column(acfQ3 = NA) |>
		add_column(essNA = NA) |>
		add_column(par_mess = NA) |>
		add_column(theta_essQ1 = NA) |>
		add_column(theta_essQ2 = NA) |>
		add_column(theta_essQ3 = NA) |>
		add_column(elapsed = NA) |>
		add_column(rejections = NA) |>
		add_column(knot_updates_burn = NA) |>
		add_column(knot_updates_keep = NA)
}

# Check the tables from each method to see how many reps they completed. In
# case these values are different, take the smallest as the checkpoint and
# continue from there.
last_mwg = df_mwg |> filter(!is.na(elapsed)) |> pull(rep) |> max(0)
last_vwg = df_vwg |> filter(!is.na(!elapsed)) |> pull(rep) |> max(0)
last_rep = min(last_mwg, last_vwg)

if (!dir.exists("fits") && save_fits) { dir.create("fits") }

# ----- Run the given level of the simulation -----

for (s in setdiff(seq_len(N_sim), seq_len(last_rep))) {
	logger("Simulation rep %d\n", s)

	# Generate data
	sigma2_true = rlnorm(m, Zgamma_true, tau_true)
	theta_true = rnorm(m, Xbeta_true, phi_true)
	s2 = sigma2_true / df * rchisq(m, df)
	y = rnorm(m, theta_true, sqrt(sigma2_true))

	# Pick starting values based on the observed data.
	lm1_out = lm(y ~ X - 1)
	lm2_out = lm(log(s2) ~ Z - 1)
	beta_init = coef(lm1_out)
	phi2_init = sigma(lm1_out)^2
	gamma_init = coef(lm2_out)
	tau2_init = sigma(lm2_out)^2
	init = get_init(m, d1, d2, beta = beta_init, gamma = gamma_init,
		sigma2 = s2, phi2 = phi2_init, tau2 = tau2_init)

	# ----- MWG -----
	logger("Running MWG\n")
	vws_ctrl = get_vws_control(method = "imh")
	control = get_gibbs_control(R = mwg_draws, burn = mwg_burn,
		thin = mwg_thin, report = mwg_report, vws = vws_ctrl,
		save_latent = seq_len(m))
	fixed = get_fixed()
	gibbs0_out = gibbs(y, s2, X, Z, df, init, control, fixed)

	if (save_fits) {
		saveRDS(gibbs0_out, sprintf("fits/mwg-%04d.rds", s))
	}

	autocorr0 = numeric(m)
	for (i in 1:m) {
		acf_out = acf(gibbs0_out$sigma2_hist[,i], lag.max = 1, plot = F)
		autocorr0[i] = acf_out$acf[2,1,1]
	}

	ess0_out = ess(gibbs0_out$sigma2_hist)
	theta_ess0_out = ess(gibbs0_out$theta_hist)
	par0_mcmc = cbind(gibbs0_out$beta_hist, gibbs0_out$gamma_hist,
		gibbs0_out$phi2_hist, gibbs0_out$tau2_hist)
	df_mwg$essQ1[s] = quantile(ess0_out, probs[1], na.rm = TRUE)
	df_mwg$essQ2[s] = quantile(ess0_out, probs[2], na.rm = TRUE)
	df_mwg$essQ3[s] = quantile(ess0_out, probs[3], na.rm = TRUE)
	df_mwg$acfQ1[s] = quantile(autocorr0, 1 - probs[1], na.rm = TRUE)
	df_mwg$acfQ2[s] = quantile(autocorr0, 1 - probs[2], na.rm = TRUE)
	df_mwg$acfQ3[s] = quantile(autocorr0, 1 - probs[3], na.rm = TRUE)
	df_mwg$essNA[s] = sum(is.na(ess0_out))
	df_mwg$par_mess[s] = multiESS(par0_mcmc)
	df_mwg$theta_essQ1[s] = quantile(theta_ess0_out, probs[1], na.rm = TRUE)
	df_mwg$theta_essQ2[s] = quantile(theta_ess0_out, probs[2], na.rm = TRUE)
	df_mwg$theta_essQ3[s] = quantile(theta_ess0_out, probs[3], na.rm = TRUE)

	df_mwg$rejections[s] = sum(gibbs0_out$sigma2_rejections_hist)
	df_mwg$elapsed[s] = sum(unlist(gibbs0_out$elapsed))

	write_csv(df_mwg, file = ff_mwg)
	gc()

	# ----- VWG -----
	logger("Running VWG\n")
	vws_ctrl = get_vws_control(tol1 = tol1, tol2 = tol2, max_rejects = 1e6,
		method = "vws-tune")
	control = get_gibbs_control(R = vwg_draws, burn = vwg_burn,
		thin = vwg_thin, report = vwg_report, vws = vws_ctrl,
		save_latent = seq_len(m))
	fixed = get_fixed()
	gibbs1_out = tryCatch({
		gibbs(y, s2, X, Z, df, init, control, fixed)
	}, error = function(e) {
		printf("VWG failed\n")
		NULL
	})

	if (save_fits) {
		saveRDS(gibbs1_out, sprintf("fits/vwg-%04d.rds", s))
	}

	if (!is.null(gibbs1_out)) {
		autocorr1 = numeric(m)
		for (i in 1:m) {
			acf_out = acf(gibbs1_out$sigma2_hist[,i], lag.max = 1, plot = F)
			autocorr1[i] = acf_out$acf[2,1,1]
		}

		ess1_out = ess(gibbs1_out$sigma2_hist)
		theta_ess1_out = ess(gibbs1_out$theta_hist)
		par1_mcmc = cbind(gibbs1_out$beta_hist, gibbs1_out$gamma_hist,
			gibbs1_out$phi2_hist, gibbs1_out$tau2_hist)
		df_vwg$essQ1[s] = quantile(ess1_out, probs[1], na.rm = TRUE)
		df_vwg$essQ2[s] = quantile(ess1_out, probs[2], na.rm = TRUE)
		df_vwg$essQ3[s] = quantile(ess1_out, probs[3], na.rm = TRUE)
		df_vwg$acfQ1[s] = quantile(autocorr1, 1 - probs[1], na.rm = TRUE)
		df_vwg$acfQ2[s] = quantile(autocorr1, 1 - probs[2], na.rm = TRUE)
		df_vwg$acfQ3[s] = quantile(autocorr1, 1 - probs[3], na.rm = TRUE)
		df_vwg$essNA[s] = sum(is.na(ess1_out))
		df_vwg$par_mess[s] = multiESS(par1_mcmc)
		df_vwg$theta_essQ1[s] = quantile(theta_ess1_out, probs[1], na.rm = TRUE)
		df_vwg$theta_essQ2[s] = quantile(theta_ess1_out, probs[2], na.rm = TRUE)
		df_vwg$theta_essQ3[s] = quantile(theta_ess1_out, probs[3], na.rm = TRUE)

		df_vwg$rejections[s] = sum(gibbs1_out$sigma2_rejections_hist)
		df_vwg$elapsed[s] = sum(unlist(gibbs1_out$elapsed))

		df_vwg$knot_updates_burn[s] =
			data.frame(updates = gibbs1_out$sigma2_knot_updates_hist) |>
			filter(row_number() <= vwg_burn) |>
			pull(updates) |>
			sum()

		df_vwg$knot_updates_keep[s] =
			data.frame(updates = gibbs1_out$sigma2_knot_updates_hist) |>
			filter(row_number() > vwg_burn) |>
			pull(updates) |>
			sum()

	}
	write_csv(df_vwg, file = ff_vwg)
	gc()
}

