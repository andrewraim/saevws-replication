#' traceplots
#'
#' Prepare traceplots for a matrix with draws.
#'
#' @param x An \eqn{R \times k} matrix where each column contains \eqn{R}
#' draws.
#' @param xlab The label for the x-coordinate. Default is an empty string.
#'
#' @return A list of \eqn{k} `ggplot` objects; the \eqn{j}th element is a
#' traceplot corresponding to the \eqn{j}th column of `x`.
#'
#' @details
#' If `x` has non-`NULL` column names, they will be used as the y-axes names in
#' the plots.
#'
#' @examples
#' set.seed(1234)
#' x = cbind(x1 = rnorm(1000), x2 = rnorm(1000))
#' out = traceplots(x)
#' print(out[[1]])
#' print(out[[2]])
#'
#' colnames(x) = c("Variable 1", "Variable 2")
#' out = traceplots(x)
#' print(out[[1]])
#' print(out[[2]])
#'
#' @export
traceplots = function(x, xlab = "")
{
	k = ncol(x)
	R = nrow(x)
	nn = colnames(x)
	if (is.null(nn)) { nn = rep("", k) }

	out = list()

	for (j in 1:k) {
		out[[j]] = data.frame(step = 1:R, draw = x[,j]) |>
			ggplot() +
			geom_line(aes(x = step, y = draw)) +
			xlab(xlab) +
			ylab(nn[j]) +
			theme_minimal()
	}

	return(out)
}

