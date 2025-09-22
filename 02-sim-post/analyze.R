library(tidyverse)
library(jsonlite)
library(knitr)

options(knitr.table.format = "latex")

# ---- Begin config -----
m_levels = c(500, 2000)
tol1_levels = c(0.5, 0.75, 0.85)
tol2_levels = c(0.0001, 0.001, 0.01)
# ---- End config -----

tbl = expand_grid(idx_m = seq_along(m_levels),
		idx_tol1 = seq_along(tol1_levels),
		idx_tol2 = seq_along(tol2_levels)) |>
	mutate(m = m_levels[idx_m]) |>
	mutate(tol1 = tol1_levels[idx_tol1]) |>
	mutate(tol2 = tol2_levels[idx_tol2])

tbl_mwg = tbl |>
	select(-idx_m, -idx_tol1, -idx_tol2) |>
	mutate(essQ1 = NA_real_, essQ2 = NA_real_, essQ3 = NA_real_,
		acfQ1 = NA_real_, acfQ2 = NA_real_, acfQ3 = NA_real_,
		essNA = NA_real_, par_mess = NA_real_, theta_essQ1 = NA_real_,
		theta_essQ2 = NA_real_, theta_essQ3 = NA_real_,
		elapsed = NA_real_, rejections = NA_real_)

tbl_vwg = tbl |>
	select(-idx_m, -idx_tol1, -idx_tol2) |>
	mutate(essQ1 = NA_real_, essQ2 = NA_real_, essQ3 = NA_real_,
		acfQ1 = NA_real_, acfQ2 = NA_real_, acfQ3 = NA_real_,
		essNA = NA_real_, par_mess = NA_real_, theta_essQ1 = NA_real_,
		theta_essQ2 = NA_real_, theta_essQ3 = NA_real_,
		elapsed = NA_real_, rejections = NA_real_,
		knot_updates_burn = NA_real_, knot_updates_keep = NA_real_)

for (i in seq_len(nrow(tbl))) {
	# Create a folder for this level of the simulation
	# If it already exists, skip it
	dd = sprintf("results/m%d_tol%d_%d", tbl$idx_m[i], tbl$idx_tol1[i],
		tbl$idx_tol2[i])
	ff_mwg = file.path(dd, "mwg.csv")
	ff_vwg = file.path(dd, "vwg.csv")
	if (!dir.exists(dd)) { next }
	if (!file.exists(ff_mwg)) { next }
	if (!file.exists(ff_vwg)) { next }

	df_mwg = ff_mwg |> read_csv(show_col_types = FALSE)
	x = df_mwg |>
		select(essQ1, essQ2, essQ3, acfQ1, acfQ2, acfQ3, essNA, par_mess,
			theta_essQ1, theta_essQ2, theta_essQ3, elapsed, rejections) |>
		colMeans()
	tbl_mwg$essQ1[i] = x["essQ1"]
	tbl_mwg$essQ2[i] = x["essQ2"]
	tbl_mwg$essQ3[i] = x["essQ3"]
	tbl_mwg$acfQ1[i] = x["acfQ1"]
	tbl_mwg$acfQ2[i] = x["acfQ2"]
	tbl_mwg$acfQ3[i] = x["acfQ3"]
	tbl_mwg$essNA[i] = x["essNA"]
	tbl_mwg$par_mess[i] = x["par_mess"]
	tbl_mwg$theta_essQ1[i] = x["theta_essQ1"]
	tbl_mwg$theta_essQ2[i] = x["theta_essQ2"]
	tbl_mwg$theta_essQ3[i] = x["theta_essQ3"]
	tbl_mwg$elapsed[i] = x["elapsed"]
	tbl_mwg$rejections[i] = x["rejections"]

	df_vwg = ff_vwg |> read_csv(show_col_types = FALSE)
	x = df_vwg |>
		select(essQ1, essQ2, essQ3, acfQ1, acfQ2, acfQ3, essNA, par_mess,
			theta_essQ1, theta_essQ2, theta_essQ3, elapsed, rejections,
			knot_updates_burn, knot_updates_keep) |>
		colMeans()
	tbl_vwg$essQ1[i] = x["essQ1"]
	tbl_vwg$essQ2[i] = x["essQ2"]
	tbl_vwg$essQ3[i] = x["essQ3"]
	tbl_vwg$acfQ1[i] = x["acfQ1"]
	tbl_vwg$acfQ2[i] = x["acfQ2"]
	tbl_vwg$acfQ3[i] = x["acfQ3"]
	tbl_vwg$essNA[i] = x["essNA"]
	tbl_vwg$par_mess[i] = x["par_mess"]
	tbl_vwg$theta_essQ1[i] = x["theta_essQ1"]
	tbl_vwg$theta_essQ2[i] = x["theta_essQ2"]
	tbl_vwg$theta_essQ3[i] = x["theta_essQ3"]
	tbl_vwg$elapsed[i] = x["elapsed"]
	tbl_vwg$rejections[i] = x["rejections"]
	tbl_vwg$knot_updates_burn[i] = x["knot_updates_burn"]
	tbl_vwg$knot_updates_keep[i] = x["knot_updates_keep"]
}

# ---- Table of sigma2 results -----

dfmt = function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")

tbl_mwg |>
	group_by(m) |>
	summarize(essQ1 = mean(essQ1), essQ2 = mean(essQ2), essQ3 = mean(essQ3),
		elapsed = mean(elapsed), rejections = mean(rejections)) |>
	mutate(m = dfmt(m)) |>
	mutate(essQ1 = dfmt(essQ1)) |>
	mutate(essQ2 = dfmt(essQ2)) |>
	mutate(essQ3 = dfmt(essQ3)) |>
	mutate(elapsed = dfmt(elapsed, 2)) |>
	mutate(tol1 = NA) |>
	mutate(tol2 = NA) |>
	mutate(knot_updates_burn = NA) |>
	mutate(knot_updates_keep = NA) |>
	select(m, tol1, tol2, essQ1, essQ2, essQ3, elapsed, rejections,
		knot_updates_burn, knot_updates_keep) |>
	kable(booktabs = T, linesep = "")

tbl_vwg |>
	mutate(m = dfmt(m)) |>
	mutate(essQ1 = dfmt(essQ1)) |>
	mutate(essQ2 = dfmt(essQ2)) |>
	mutate(essQ3 = dfmt(essQ3)) |>
	mutate(elapsed = dfmt(elapsed, 2)) |>
	mutate(rejections = dfmt(rejections)) |>
	mutate(knot_updates_burn = dfmt(knot_updates_burn)) |>
	mutate(knot_updates_keep = dfmt(knot_updates_keep)) |>
	select(m, tol1, tol2, essQ1, essQ2, essQ3, elapsed, rejections,
		knot_updates_burn, knot_updates_keep) |>
	kable(booktabs = T, linesep = "")

# ---- Table of other results -----

tbl_mwg |>
	group_by(m) |>
	summarize(essQ1 = mean(theta_essQ1), essQ2 = mean(theta_essQ2),
		essQ3 = mean(theta_essQ3), par_mess = mean(par_mess)) |>
	mutate(m = dfmt(m)) |>
	mutate(essQ1 = dfmt(essQ1)) |>
	mutate(essQ2 = dfmt(essQ2)) |>
	mutate(essQ3 = dfmt(essQ3)) |>
	mutate(par_mess = dfmt(par_mess)) |>
	mutate(tol1 = NA) |>
	mutate(tol2 = NA) |>
	select(m, tol1, tol2, essQ1, essQ2, essQ3, par_mess) |>
	kable(booktabs = T, linesep = "")

tbl_vwg |>
	mutate(m = dfmt(m)) |>
	mutate(essQ1 = dfmt(theta_essQ1)) |>
	mutate(essQ2 = dfmt(theta_essQ2)) |>
	mutate(essQ3 = dfmt(theta_essQ3)) |>
	mutate(par_mess = dfmt(par_mess)) |>
	select(m, tol1, tol2, essQ1, essQ2, essQ3, par_mess) |>
	kable(booktabs = T, linesep = "")
