#' Equally Spaced Univariate Knots
#'
#' Prepare knots which define \eqn{N} equally-spaced intervals between finite
#' endpoints `lo` and `hi`.
#'
#' @param N Number of desired intervals.
#' @param lo Left endpoint; must be finite.
#' @param hi Right endpoint; must be finite.
#'
#' @return A vector of \eqn{N-1} knots.
#'
#' @examples
#' equal_spaced_knots(2, lo = 0, hi = 1)
#' equal_spaced_knots(3, lo = 0, hi = 1)
#'
#' @export
equal_spaced_knots = function(N, lo, hi)
{
	stopifnot(N > 0)
	x = seq(lo, hi, length.out = N + 1)
	x |> head(n = -1) |> tail(n = -1)
}
