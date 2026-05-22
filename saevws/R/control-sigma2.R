#' Gibbs Sampler VWS Control
#'
#' @param tol_suff The tolerance \eqn{\epsilon_1}.
#' @param tol_merge The tolerance \eqn{\epsilon_2}.
#' @param max_rejects Maximum number of rejections allowed per \eqn{\sigma_i^2}
#' conditional per Gibbs step. The sampler halts when this is exceeded.
#' @param method Can be independent Metropolis step (`"imh"`), self-tuned VWS
#' (`"vws-tune"`), or basic VWS without self-tuning (`"vws-basic"`) for the
#' \eqn{\sigma_i^2} draws.
#'
#' @return A list with results.
#'
#' @examples
#' ctrl = control_sigma2()
#'
#' @export
control_sigma2 = function(tol_suff = 1e-2, tol_merge = exp(-100),
	max_rejects = 1e6, method = c("vws-tune", "vws-basic", "imh", "arms"),
	N = 50)
{
	ret = list(tol_suff = tol_suff, tol_merge = tol_merge,
		max_rejects = max_rejects, method = match.arg(method), N = N)
	class(ret) = "control_vws"
	return(ret)
}

