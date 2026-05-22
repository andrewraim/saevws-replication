#' Control for Joint Model Gibbs Sampler
#'
#' @param R Desired length of MCMC chain.
#' @param burn Number of draws to burn.
#' @param thin Thinning factor for draws to save.
#' @param report Determines how often progress of the sampler is
#' reported.
#' @param save_latent integer vector; specify indices of the \eqn{m}
#' observations whose latent draws will be saved (\eqn{\sigma^2} and
#' \eqn{\vartheta}). Values should be 1-based, corresponding to a subset of
#' \eqn{\{1, \ldots, m\}}. Default is an empty vector. Saving many observations
#' over many draws can use a lot of memory.
#' @param sigma2 An control object obtained from [control_sigma2].
#'
#' @return A list with results.
#'
#' @examples
#' ctrl = control_joint()
#'
#' @export
control_joint = function(R = 1000, burn = 0, thin = 1, report = R+1,
	save_latent = integer(0), sigma2 = control_sigma2())
{
	ret = list(R = R, burn = burn, thin = thin, report = report,
		save_latent = save_latent, sigma2 = sigma2)
	class(ret) = "control_joint"
	return(ret)
}

#' Gibbs Sampler Fixed Components
#'
#' @param beta logical; if `TRUE`, Gibbs sampler will leave \eqn{\beta} fixed
#' in MCMC.
#' @param gamma logical; if `TRUE`, Gibbs sampler will leave \eqn{\gamma} fixed
#' in MCMC.
#' @param phi2 logical; if `TRUE`, Gibbs sampler will leave \eqn{\phi^2} fixed
#' in MCMC.
#' @param tau2 logical; if `TRUE`, Gibbs sampler will leave \eqn{\tau^2} fixed
#' in MCMC.
#' @param sigma2 logical; if `TRUE`, Gibbs sampler will leave \eqn{\sigma^2}
#' fixed in MCMC.
#' @param theta logical; if `TRUE`, Gibbs sampler will leave \eqn{\vartheta}
#' fixed in MCMC.
#'
#' @return A list with results.
#'
#' @examples
#' fixed = fixed_joint()
#'
#' @export
fixed_joint = function(beta = FALSE, gamma = FALSE, phi2 = FALSE, tau2 = FALSE,
	sigma2 = FALSE, theta = FALSE)
{
	ret = list(beta = beta, gamma = gamma, phi2 = phi2, tau2 = tau2,
		sigma2 = sigma2, theta = theta)
	class(ret) = "fixed_joint"
	return(ret)
}

#' Gibbs Sampler Initial Values
#'
#' @param m Number of subjects.
#' @param d1 Dimension of \eqn{X} matrix.
#' @param d2 Dimension of \eqn{Z} matrix.
#' @param beta Initial value for \eqn{\beta}.
#' @param gamma Initial value for \eqn{\gamma}.
#' @param phi2 Initial value for \eqn{\phi^2}.
#' @param tau2 Initial value for \eqn{\tau^2}.
#' @param sigma2 Initial value for \eqn{\sigma^2}.
#' @param theta Initial value for \eqn{\vartheta^2}.
#'
#' @return A list with results.
#'
#' @examples
#' init = init_joint(500, d1 = 5, d2 = 2)
#'
#' @export
init_joint = function(m, d1, d2, beta = NULL, gamma = NULL, phi2 = NULL,
	tau2 = NULL, sigma2 = NULL, theta = NULL)
{
	if (is.null(beta)) { beta = numeric(d1)	}
	if (is.null(gamma)) { gamma = numeric(d2) }
	if (is.null(phi2)) { phi2 = 1 }
	if (is.null(tau2)) { tau2 = 1 }
	if (is.null(sigma2)) { sigma2 = rep(1, m) }
	if (is.null(theta)) { theta = numeric(m) }

	stopifnot(length(beta) == d1)
	stopifnot(length(gamma) == d2)
	stopifnot(length(phi2) == 1)
	stopifnot(length(tau2) == 1)
	stopifnot(length(sigma2) == m)
	stopifnot(length(theta) == m)

	ret = list(beta = beta, gamma = gamma, phi2 = phi2, tau2 = tau2,
		sigma2 = sigma2, theta = theta)
	class(ret) = "init_joint"
	return(ret)
}

#' Gibbs Sampler
#'
#' Run the Gibbs sampler.
#'
#' @param y Observed point estimates.
#' @param s2 Observed variance estimates.
#' @param X Design matrix for regression on point estimates.
#' @param Z Design matrix for regression on variance estimates.
#' @param df Degrees of freedom to use for variance estimates.
#' @param init Initial values from [init_joint].
#' @param control Control object from [control_joint].
#' @param fixed Fixed value indicators from [fixed_joint].
#'
#' @return A list with results from the sampler.
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(1234)
#'
#' m = 500
#' tau_true = sqrt(0.25)
#' phi_true = sqrt(0.2)
#' beta_true = c(1.5, 0.85)
#' gamma_true = c(2.6, -1)
#' df = rchisq(m, 16)
#' X = cbind(1, rnorm(m, 8, 2))
#' Z = cbind(1, rnorm(m, 7, 1.25))
#'
#' sigma2_true = rlnorm(m, Z %*% gamma_true, tau_true)
#' theta_true = rnorm(m, X %*% beta_true, phi_true)
#' s2 = sigma2_true / df * rchisq(m, df)
#' y = rnorm(m, theta_true, sqrt(s2))
#'
#' ctrl = control_joint(R = 100, report = 20)
#' gibbs_out = gibbs_joint(y, s2, X, Z, df, control = ctrl)
#' }
#'
#' @export
gibbs_joint = function(y, s2, X, Z, df,
	init = init_joint(m = length(y), d1 = ncol(X), d2 = ncol(Z)),
	control = control_joint(), fixed = fixed_joint())
{
	m = length(y)
	stopifnot(m == length(s2))
	stopifnot(m == length(df))
	stopifnot(m == nrow(X))
	stopifnot(m == nrow(Z))

	## Ensure several conditions for indices:
	## - No duplicate values,
	## - All values are in 1:n.
	## Convert to 0-based indices to pass to C++
	save_latent = control$save_latent
	stopifnot(table(save_latent) == 1)
	stopifnot(all(save_latent %in% 1:m))
	control$save_latent = save_latent - 1

	out = gibbs_joint_cpp(y, s2, X, Z, df, init, control, fixed)
	class(out) = "joint_fit"
	return(out)
}

#' Gibbs Sampler Summary
#'
#' @param object A result from [gibbs_joint].
#' @param pr Vector of quantiles to present in summary.
#' @param ... Additional arguments.
#'
#' @return A data frame with results.
#'
#' @export
summary.joint_fit = function(object, pr = c(0.05, 0.95), ...)
{
	d1 = ncol(object$beta_hist)
	d2 = ncol(object$gamma_hist)

	df_beta = as.data.frame(cbind(
		apply(object$beta_hist, 2, mean),
		apply(object$beta_hist, 2, sd),
		t(apply(object$beta_hist, 2, quantile, probs = pr)),
		apply(object$beta_hist, 2, ess)
	))
	rownames(df_beta) = sprintf("beta%d", 1:d1)

	df_gamma = as.data.frame(cbind(
		apply(object$gamma_hist, 2, mean),
		apply(object$gamma_hist, 2, sd),
		t(apply(object$gamma_hist, 2, quantile, probs = pr)),
		apply(object$gamma_hist, 2, ess)
	))
	rownames(df_gamma) = sprintf("gamma%d", 1:d2)

	df_phi2 = as.data.frame(cbind(
		mean(object$phi2_hist),
		sd(object$phi2_hist),
		t(quantile(object$phi2_hist, probs = pr)),
		ess(object$phi2_hist)
	))
	rownames(df_phi2) = sprintf("phi2")

	df_tau2 = as.data.frame(cbind(
		mean(object$tau2_hist),
		sd(object$tau2_hist),
		t(quantile(object$tau2_hist, probs = pr)),
		ess(object$tau2_hist)
	))
	rownames(df_tau2) = sprintf("tau2")

	df = rbind(df_beta, df_gamma, df_phi2, df_tau2)
	quantile_names = sprintf("%g%%", 100 * pr)
	colnames(df) = c("mean", "sd", quantile_names, "ess")
	idx1 = 1:2
	idx2 = seq_along(quantile_names) + 2
	idx3 = length(quantile_names) + 3
	df[,idx1] = round(df[,idx1], 4)
	df[,idx2] = round(df[,idx2], 4)
	df[,idx3] = round(df[,idx3], 2)
	return(df)
}

#' Gibbs Sampler Print Summary
#'
#' @param x A result from [gibbs_joint].
#' @param pr Vector of quantiles to present in summary.
#' @param ... Additional arguments.
#'
#' @export
print.joint_fit = function(x, pr = c(0.05, 0.95), ...)
{
	cat("Summary of fit for Joint SAE model\n")
	print(summary(x, pr))

	cat("----\n")
	printf("Total iterations R: %d   Burn: %d   Thin: %d   Saved draws: %d\n",
		x$R, x$burn, x$thin, x$R_keep)

	printf("Rejections in sigma2 step: %d  Total proposals: %d  Rejection rate: %g%%\n",
		sum(x$sigma2_rejections_hist),
		sum(x$sigma2_rejections_hist) + x$R * x$m,
		100 * sum(x$sigma2_rejections_hist) / (sum(x$sigma2_rejections_hist) + x$R*x$m))

	printf("Avg regions in sigma2 step: %g\n", sum(x$sigma2_knots_hist) / (x$R * x$m))

	cat("----\n")
	printf("Elapsed time (Seconds):\n")
	x$elapsed$total = sum(unlist(x$elapsed))
	tab = round(as.data.frame(x$elapsed), 4)
	rownames(tab) = ""
	print(tab)
}

