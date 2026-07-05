#' Gibbs Sampler Inner Control
#'
#' Parameters for within-Gibbs samplers: independent Metropolis-Hastings (IMH),
#' Adaptive Metropolis-Hastings (AMH), vertical weighted strips (VWS) with
#' and without self-tuning, and Adaptive Rejection Metropolis Sampling (ARMS).
#'
#' @param tol_suff The tolerance \eqn{\epsilon_1}.
#' @param tol_merge The tolerance \eqn{\epsilon_2}.
#' @param max_rejects Maximum number of rejections allowed per within-Gibbs
#' step. The Gibbs sampler halts when this is exceeded.
#' @param method Can be IMH (`"imh"`), AMH (`"amh"`), self-tuned VWS
#' (`"vws-tune"`), basic VWS without self-tuning (`"vws-basic"`), and
#' AMRS (`"arms"`).
#' @param N Maximum number of refinements for a VWS proposal without self-tuning.
#' @param tune Maximum number of iterations in the Gibbs sampler to continue
#' tuning for self-tuned VWS.
#' @param amh_varprop_init Initial value \eqn{v_j^{(0)}} of proposal variance
#' for AMH.
#' @param amh_varprop_c Multiplier \eqn{s} for proposal variance in AMH. The
#' default 2.4 is recommended by Roberts & Rosenthal (2009).
#' @param amh_varprop_eps A value \eqn{\epsilon} to add to proposal variance
#' for AMH to ensure positivity.
#'
#' @return A list with the settings.
#'
#' @details
#' See the replication materials guide document for details about AMH.
#'
#' @examples
#' ctrl = control_inner()
#'
#' @references
#' Roberts, Gareth O., and Jeffrey S. Rosenthal. 2009. Examples of Adaptive
#' MCMC. Journal of Computational and Graphical Statistics 18 (2): 349-67.
#'
#' @export
control_inner = function(tol_suff = 1e-2, tol_merge = exp(-100),
	max_rejects = 1e6,
	method = c("vws-tune", "vws-basic", "imh", "arms", "amh"),
	N = 50, tune = 1e6, amh_varprop_init = 25,  amh_varprop_c = 2.4,
	amh_varprop_eps = 0.05)
{
	ret = list(tol_suff = tol_suff, tol_merge = tol_merge,
		max_rejects = max_rejects, method = match.arg(method), N = N,
		tune = tune,
		amh_varprop_init = amh_varprop_init, amh_varprop_c = amh_varprop_c,
		amh_varprop_eps = amh_varprop_eps)
	class(ret) = "control_inner"
	return(ret)
}
