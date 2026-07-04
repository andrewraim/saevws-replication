#' Where Leveled
#'
#' A utility to plot a series which increases dramatically at the beginning and
#' then tends to level off but still contain variability. The objective is to
#' find an index where the dramatic increase is beginning to level off. The
#' caller can use this to zoom into the leveled-off data. Similarly, it can be
#' used for series which decrease.
#'
#' @param x Numeric vector containing the series of interest.
#' @param burn Integer scalar; burn-in period before considering the series to
#' have leveled off.
#' @param tol Numeric scalar that determines how much of the initial dramatic
#' @param increasing Logical scalar; `TRUE` if the series is expected to
#' increase and then level off, `FALSE` to instead handle series which are
#' expected to decrease and then level off.
#'
#' @returns An integer
#'
#' @details
#' Suppose the series `x` is expected to increase dramatically and then level
#' off (the caller should confirm this). Let `x0` be the minimum value of `x`
#' that occurs after index `burn`. Return the first index of `x` whose value is
#' larger than `(1 - tol) * x0)`.
#'
#' Now suppose `x` expected to decrease dramatically and then level off. Let
#' `x0` be the maximum value of `x` that occurs after index `burn`. Return the
#' first index of `x` whose value is smaller than `(1 + tol) * x0)`.
#'
#' @examples
#' set.seed(1234)
#'
#' # Generate a series with that decreases to a stable value, with variability
#' n = 1000
#' iter = 1:n
#' x = exp(-(iter - 10) / 10) + rnorm(n, 0, 0.05)
#'
#' # idx1 gives an index where the series has mostly descended to stability.
#' # idx2 gives an earlier index.
#' idx1 = where_leveled(x, burn = 200, tol = 0.01, increasing = FALSE)
#' idx2 = where_leveled(x, burn = 200, tol = 10, increasing = FALSE)
#'
#' plot(x, type = "l")
#' abline(v = idx1, lty = 2, col = "red")
#' abline(v = idx2, lty = 2, col = "red")
#'
#'
#' @export
where_leveled = function(x, burn, tol, increasing = TRUE)
{
	idx_burn = seq(burn + 1, length(x))
	if (increasing) {
		x0 = min(x[idx_burn])
		idx = which(x > (1 - tol) * x0)
	} else {
		x0 = max(x[idx_burn])
		idx = which(x < (1 + tol) * x0)
	}

	out = max(idx[1], 0)
	return(out)
}

#' Plot Number of Tunes
#'
#' Make a formatted ggplot with the number of tunes (adjustments to adaptive
#' proposals) per iteration.
#'
#' @param x Numeric vector; a series with the number of tunes per iteration.
#' @param burn Integer scalar; used to call [where_leveled].
#' @param tol Numeric scalar; used to call [where_leveled].
#'
#' @return A ggplot.
#'
#' @export
plot_tunes = function(x, burn, tol = 0.10)
{
	R = length(x)
	iter_first = where_leveled(x, burn, tol, increasing = FALSE)
	ymin = min(x[iter_first:R])
	ymax = max(x[iter_first:R])

	out = data.frame(updates = x) %>%
		mutate(iter = row_number()) %>%
		filter(iter > iter_first) %>%
		ggplot() +
		geom_line(aes(iter, updates)) +
		xlab("Iteration") +
		ylab("Tunes") +
		scale_x_continuous(
			n.breaks = 5,
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

#' Plot Number of Components
#'
#' Make a formatted ggplot with the number of components (mixture components
#' in finite mixture proposals) per iteration.
#'
#' @param x Numeric vector; a series with the number of components per iteration.
#' @param burn Integer scalar; used to call [where_leveled].
#' @param tol Numeric scalar; used to call [where_leveled].
#'
#' @return A ggplot.
#'
#' @export
plot_comps = function(x, burn, tol = 0.01)
{
	R = length(x)
	iter_first = where_leveled(x, burn, tol, increasing = TRUE)
	ymin = min(x[iter_first:R])
	ymax = max(x[iter_first:R])

	out = data.frame(count = x) %>%
		mutate(iter = row_number()) %>%
		filter(iter >= iter_first) %>%
		ggplot() +
		geom_line(aes(iter, count)) +
		xlab("Iteration") +
		ylab("Regions") +
		scale_x_continuous(
			n.breaks = 5,
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

#' Plot Number of Rejections
#'
#' Make a formatted ggplot with the number of rejections per iteration.
#'
#' @param x Numeric vector; a series with the number of rejections per iteration.
#' @param burn Integer scalar; used to call [where_leveled].
#' @param tol Numeric scalar; used to call [where_leveled].
#'
#' @return A ggplot.
#'
#' @export
plot_rejects = function(x, burn, tol = 0.20)
{
	R = length(x)
	iter_first = where_leveled(x, burn, tol, increasing = FALSE)
	ymin = min(x[iter_first:R])
	ymax = max(x[iter_first:R])

	out = data.frame(count = x) %>%
		mutate(iter = row_number()) %>%
		filter(iter > iter_first) %>%
		ggplot() +
		geom_line(aes(iter, count)) +
		xlab("Iteration") +
		ylab("Rejections") +
		scale_x_continuous(
			n.breaks = 5,
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
