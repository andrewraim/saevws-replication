Rcpp::sourceCpp("samplers.cpp")

d_invgamma = function(x, shape, scale, log = FALSE)
{
	n = length(x)
	out = rep(-Inf, n)

	lnc = base::lgamma(shape) - shape * base::log(scale)

	for (i in 1:n) {
		if (x[i] > 0) {
			out[i] = -(shape + 1)*base::log(x[i]) - scale / x[i] - lnc
		}
	}
	if (log) { return(out) } else { return(exp(out)) }
}

d_target_unnorm = function(x, mu, tau, kappa, lambda, log = FALSE)
{
	n = length(x)
	out = rep(-Inf, n)
	for (i in 1:n) {
		if (x[i] > 0) {
			out[i] = d_invgamma(x[i], kappa, lambda, log = TRUE) +
				dlnorm(x[i], mu, tau, log = TRUE)
		}
	}
	if (log) { return(out) } else { return(exp(out)) }
}

n_target = function(mu, tau, kappa, lambda, log = FALSE)
{
	f = function(x) { d_target_unnorm(x, mu, tau, kappa, lambda) }
	int_out = integrate(f, lower = 0, upper = Inf)
	out = base::log(int_out$value)
	if (log) { return(out) } else { return(exp(out)) }
}

d_target = function(x, mu, tau, kappa, lambda, log = FALSE)
{
	lnc = n_target(mu, tau, kappa, lambda, log = TRUE)
	out = d_target_unnorm(x, mu, tau, kappa, lambda, log = TRUE) - lnc
	if (log) { return(out) } else { return(exp(out)) }
}

w_target = function(x, kappa, lambda, log = TRUE)
{
	n = length(x)
	out = rep(-Inf, n)
	for (i in 1:n) {
		if (x[i] > 0) {
			out[i] = d_invgamma(x[i], kappa, lambda, log = TRUE)
		}
	}
	if (log) { return(out) } else { return(exp(out)) }
}

# For plotting a function that is generally increasing, find the first index
# where we cross tol percent of the max. If the function is generally
# increasing, find the first index where the min is larger than tol percent of
# the series.
first_cross = function(x, burn, tol, increasing = TRUE)
{
	idx_burn = seq(burn + 1, length(x))
	if (increasing) {
		idx = which(x > (1 - tol) * min(x[idx_burn]))
	} else {
		idx = which(x < (1 + tol) * max(x[idx_burn]))
	}

	return(idx[1])
}

plot_tunes = function(x, burn, tol = 0.10)
{
	R = length(x)
	iter_first = first_cross(x, burn, tol, increasing = FALSE)
	ymin = min(x[iter_first:R])
	ymax = max(x[iter_first:R])

	out = data.frame(updates = x) %>%
		mutate(iter = row_number()) %>%
		filter(iter > iter_first) %>%
		ggplot() +
		geom_line(aes(iter, updates)) +
		xlab("Iteration") +
		ylab("Number of Region Updates") +
		scale_x_continuous(
			n.breaks = 9,
			expand = expansion(mult = c(0, 0.05), add = c(1, 0))) +
		scale_y_continuous(
			n.breaks = 10,
			expand = expansion(mult = 0, add = 1)) +
		theme_light()

	if (iter_first > 1) {
		out = out +
			annotate("rect", xmin = 1, xmax = iter_first, ymin = ymin,
				ymax = ymax, fill = "gray", alpha = 0.5)
	}

	return(out)
}

plot_comps = function(x, burn, tol = 0.01)
{
	R = length(x)
	iter_first = first_cross(x, burn, tol, increasing = TRUE)
	ymin = min(x[iter_first:R])
	ymax = max(x[iter_first:R])

	out = data.frame(count = x) %>%
		mutate(iter = row_number()) %>%
		filter(iter >= iter_first) %>%
		ggplot() +
		geom_line(aes(iter, count)) +
		xlab("Iteration") +
		ylab("Number of Regions") +
		scale_x_continuous(
			n.breaks = 9,
			expand = expansion(mult = c(0, 0.05), add = c(1, 0))) +
		scale_y_continuous(
			n.breaks = 10,
			expand = expansion(mult = 0, add = 1)) +
		theme_light()

	if (iter_first > 1) {
		out = out +
			annotate("rect", xmin = 1, xmax = iter_first, ymin = ymin,
				ymax = ymax, fill = "gray", alpha = 0.5)
	}

	return(out)
}

plot_rejects = function(x, burn, tol = 0.20)
{
	R = length(x)
	iter_first = first_cross(x, burn, tol, increasing = FALSE)
	ymin = min(x[iter_first:R])
	ymax = max(x[iter_first:R])

	out = data.frame(count = x) %>%
		mutate(iter = row_number()) %>%
		filter(iter > iter_first) %>%
		ggplot() +
		geom_line(aes(iter, count)) +
		xlab("Iteration") +
		ylab("Number of Rejections") +
		scale_x_continuous(
			n.breaks = 9,
			expand = expansion(mult = c(0, 0.05), add = c(1, 0))) +
		scale_y_continuous(
			n.breaks = 10,
			expand = expansion(mult = 0, add = 1)) +
		theme_light()

	if (iter_first > 1) {
		out = out +
			annotate("rect", xmin = 1, xmax = iter_first, ymin = ymin,
				ymax = ymax, fill = "gray", alpha = 0.5)
	}

	return(out)
}

