#' Gibbs Sampler Inner Control
#'
#' @param tol_suff The tolerance \eqn{\epsilon_1}.
#' @param tol_merge The tolerance \eqn{\epsilon_2}.
#' @param max_rejects Maximum number of rejections allowed per \eqn{\sigma_i^2}
#' conditional per Gibbs step. The sampler halts when this is exceeded.
#' @param method Can be independent Metropolis step (`"imh"`), self-tuned VWS
#' (`"vws-tune"`), basic VWS without self-tuning (`"vws-basic"`) for the
#' \eqn{\sigma_i^2} draws, Adaptive Rejection Metropolis Sampling (`"arms"`),
#' or Adaptive Metropolis (`"am"`).
#' @param N Maximum number of refinements for VWS without self-tuning.
#' @param am_varprop_init TBD
#' @param am_varprop_c TBD
#' @param am_varprop_eps TBD
#' @param tune TBD
#'
#' @return A list with results.
#'
#' @examples
#' ctrl = control_inner()
#'
#' @export
control_inner = function(tol_suff = 1e-2, tol_merge = exp(-100),
	max_rejects = 1e6,
	method = c("vws-tune", "vws-basic", "imh", "arms", "am", "mh-vws"),
	N = 50, tune = 1e6, am_varprop_init = 25,
	am_varprop_c = 2.4, am_varprop_eps = 0.05)
{
	ret = list(tol_suff = tol_suff, tol_merge = tol_merge,
		max_rejects = max_rejects, method = match.arg(method), N = N,
		tune = tune,
		am_varprop_init = am_varprop_init, am_varprop_c = am_varprop_c,
		am_varprop_eps = am_varprop_eps)
	class(ret) = "control_inner"
	return(ret)
}
