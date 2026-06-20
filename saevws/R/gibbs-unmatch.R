#' Control for Unmatched Model Gibbs Sampler
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
#' @param inner An control object obtained from [control_inner].
#'
#' @return A list with results.
#'
#' @examples
#' ctrl = control_unmatch()
#'
#' @export
control_unmatch = function(R = 1000, burn = 0, thin = 1, report = R+1,
	save_latent = integer(0), inner = control_inner())
{
	ret = list(R = R, burn = burn, thin = thin, report = report,
		save_latent = save_latent, inner = inner)
	class(ret) = "control_unmatch"
	return(ret)
}

#' Gibbs Sampler Fixed Components for Unmatched Model
#'
#' @param beta logical; if `TRUE`, Gibbs sampler will leave \eqn{\beta} fixed
#' in MCMC.
#' @param tau2 logical; if `TRUE`, Gibbs sampler will leave \eqn{\tau^2} fixed
#' in MCMC.
#' @param mu logical; if `TRUE`, Gibbs sampler will leave \eqn{\mu}
#' fixed in MCMC.
#'
#' @return A list with results.
#'
#' @examples
#' fixed = fixed_unmatch()
#'
#' @export
fixed_unmatch = function(beta = FALSE, tau2 = FALSE, mu = FALSE)
{
	ret = list(beta = beta, tau2 = tau2, mu = mu)
	class(ret) = "fixed_unmatch"
	return(ret)
}

#' Gibbs Sampler Initial Values for Unmatched Model
#'
#' @param m Number of subjects.
#' @param d Dimension of \eqn{X} matrix.
#' @param beta Initial value for \eqn{\beta}.
#' @param tau2 Initial value for \eqn{\tau^2}.
#' @param mu Initial value for \eqn{\mu}.
#'
#' @return A list with results.
#'
#' @examples
#' init = init_unmatch(500, d1 = 5, d2 = 2)
#'
#' @export
init_unmatch = function(m, d, beta = NULL, tau2 = NULL, mu = NULL)
{
	if (is.null(beta)) { beta = numeric(d)	}
	if (is.null(tau2)) { tau2 = 1 }
	if (is.null(mu)) { mu = rep(1, m) }

	stopifnot(length(beta) == d)
	stopifnot(length(tau2) == 1)
	stopifnot(length(mu) == m)

	ret = list(beta = beta, tau2 = tau2, mu = mu)
	class(ret) = "init_unmatch"
	return(ret)
}

#' Gibbs Sampler for Unmatched Model
#'
#' Run the Gibbs sampler.
#'
#' @param y Observed point estimates.
#' @param X Design matrix for regression on point estimates.
#' @param sigma Fixed variances.
#' @param init Initial values from [init_unmatch].
#' @param control Control object from [control_unmatch].
#' @param fixed Fixed value indicators from [fixed_unmatch].
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
gibbs_unmatch = function(y, sigma, X,
	init = init_unmatch(m = length(y), d = ncol(X)),
	control = control_unmatch(), fixed = fixed_unmatch())
{
	m = length(y)
	stopifnot(m == length(sigma))
	stopifnot(m == nrow(X))

	## Ensure several conditions for indices:
	## - No duplicate values,
	## - All values are in 1:n.
	## Convert to 0-based indices to pass to C++
	save_latent = control$save_latent
	stopifnot(table(save_latent) == 1)
	stopifnot(all(save_latent %in% 1:m))
	control$save_latent = save_latent - 1

	out = gibbs_unmatch_cpp(y, sigma, X, init, control, fixed)
	class(out) = "gibbs_unmatch"
	return(out)
}

#' Gibbs Sampler Summary for Unmatched Model
#'
#' @param object A result from [gibbs_unmatch].
#' @param pr Vector of quantiles to present in summary.
#' @param ... Additional arguments.
#'
#' @return A data frame with results.
#'
#' @export
summary.gibbs_unmatch = function(object, pr = c(0.05, 0.95), ...)
{
	d = ncol(object$beta)

	df_beta = as.data.frame(cbind(
		apply(object$beta, 2, mean),
		apply(object$beta, 2, sd),
		t(apply(object$beta, 2, quantile, probs = pr)),
		apply(object$beta, 2, ess)
	))
	rownames(df_beta) = sprintf("beta%d", 1:d)

	df_tau2 = as.data.frame(cbind(
		mean(object$tau2),
		sd(object$tau2),
		t(quantile(object$tau2, probs = pr)),
		ess(object$tau2)
	))
	rownames(df_tau2) = sprintf("tau2")

	df = rbind(df_beta, df_tau2)
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

#' Gibbs Sampler Print Summary for Unmatched Model
#'
#' @param x A result from [gibbs_unmatch].
#' @param pr Vector of quantiles to present in summary.
#' @param ... Additional arguments.
#'
#' @export
print.gibbs_unmatch = function(x, pr = c(0.05, 0.95), ...)
{
	cat("Summary of fit for Unmatched SAE model\n")
	print(summary(x, pr))

	cat("----\n")
	printf("Iterations: %d  Burn: %d  Thin: %d  Saved draws: %d\n",
		x$R, x$burn, x$thin, x$R_keep)

	cat("----\n")
	printf("mu step\n")
	printf("   Proposed: %d  Rejected: %d\n",
		sum(x$mu_rejects) + x$R * x$m,
		sum(x$mu_rejects))

	printf("   Rejection rate: %g%%\n",
		100 * sum(x$mu_rejects) / (sum(x$mu_rejects) + x$R*x$m))

	if (x$inner_method == "vws-tune") {
		printf("   Avg regions: %g\n", sum(x$mu_comps) / (x$R * x$m))
	}

	cat("----\n")
	printf("Elapsed time (sec):\n")
	x$elapsed$total = sum(unlist(x$elapsed))
	tab = round(as.data.frame(x$elapsed), 4)
	rownames(tab) = ""
	print(tab)
}
