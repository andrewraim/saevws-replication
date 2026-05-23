library(tidyverse)
library(jsonlite)

# ---- Begin config -----
m_levels = c(500, 2000)
#tol1_levels = c(0.60, 0.75, 0.90)
tol1_levels = c(0.5, 0.75, 0.85)
tol2_levels = c(0.0001, 0.001, 0.01)
# ---- End config -----

tbl = expand_grid(idx_m = seq_along(m_levels),
		idx_tol1 = seq_along(tol1_levels),
		idx_tol2 = seq_along(tol2_levels)) %>%
	mutate(m = m_levels[idx_m]) %>%
	mutate(tol1 = tol1_levels[idx_tol1]) %>%
	mutate(tol2 = tol2_levels[idx_tol2])

if (!dir.exists("results")) {
	dir.create("results")
}

for (i in seq_len(nrow(tbl))) {
	# Create a folder for this level of the simulation
	# If it already exists, skip it
	dd = sprintf("results/m%d_tol%d_%d", tbl$idx_m[i], tbl$idx_tol1[i],
		tbl$idx_tol2[i])
	if (dir.exists(dd)) { next }
	dir.create(dd)

	# Generate launch script and save to file launch.R in run folder
	# The script sets variables particular to this level of the simulation,
	# then calls sim.R to run the simulation.
	list(m = tbl$m[i],
		tol1 = tbl$tol1[i],
		tol2 = tbl$tol2[i]) |>
		toJSON(pretty = TRUE) |>
		cat(file = sprintf("%s/config.json", dd))

	paste("source(\"../../sim.R\")\n") |>
		cat(file = sprintf("%s/launch.R", dd))
}

