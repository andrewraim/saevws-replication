#' Univariate Distributions Truncated to an Interval
#'
#' Functions to truncate a given univariate distribution to an interval of the
#' form \eqn{(\text{lo}, \text{hi}]}.
#'
#' @param n Desired sample size.
#' @param x Quantile or argument of density.
#' @param log logical; if `TRUE`, probabilities `p` are given as
#' \eqn{\log(p)}.
#' @param log.p logical; if `TRUE`, probabilities `p` are given as
#' \eqn{\log(p)}.
#' @param p Probability.
#' @param lo Lower limit.
#' @param hi Upper limit.
#' @param df Density function for the untruncated distribution.
#' @param pf Cumulative distribution function for the untruncated distribution.
#' @param qf Quantile function for the untruncated distribution.
#' @param tol A tolerance used to decide whether to upper or lower tailed
#' probabilities internally.
#' @param ... Additional arguments.
#'
#' @details
#' This code was initially copied from the `rtruncated` function in the
#' LearnBayes package.
#'
#' @examples
#' # Take a large sample from truncated standard normal
#' x = r_truncated(n = 100000, lo = 1, hi = 2, pf = pnorm, qf = qnorm)
#'
#' # Compare density to histogram of draws
#' hist(x, probability = TRUE)
#' curve(d_truncated(x, lo = 1, hi = 2, df = dnorm, pf = pnorm), add = TRUE)
#'
#' # Compare CDF to empirical CDF of draws
#' plot(ecdf(x))
#' curve(p_truncated(x, lo = 1, hi = 2, pf = pnorm), col = "blue", lwd = 2, add = TRUE)
#'
#' # Compare quantile function to quantiles of draws
#' pr = seq(0, 1, length.out = 50)
#' q_emp = quantile(x, probs = pr)
#' q_thry = q_truncated(pr, lo = 1, hi = 2, pnorm, qf = qnorm)
#' plot(q_thry, q_emp)
#' abline(c(0,1), lty = 2, col = "red")
#'
#' @name TruncatedUnivariate
NULL

#' @name TruncatedUnivariate
#' @export
d_truncated = function(x, lo, hi, df, pf, log = FALSE, tol = 1e-6, ...)
{
	n = length(x)
	out = numeric(n)
	if (length(lo) == 1) { lo = rep(lo, n) }
	if (length(hi) == 1) { hi = rep(hi, n) }

	log_p_lo = pf(lo, log.p = TRUE, ...)
	idx1 = which(log_p_lo > log1p(-tol))
	idx2 = setdiff(1:n, idx1)

	# If we are truncating way into the upper tail of a distribution, working
	# with the complement of the CDF helps to retain precision.
	log_cp_lo = pf(lo[idx1], log.p = TRUE, lower.tail = FALSE, ...)
	log_cp_hi = pf(hi[idx1], log.p = TRUE, lower.tail = FALSE, ...)
	out[idx1] = df(x[idx1], log = TRUE, ...) +
		log(lo[idx1] < x[idx1] & x[idx1] <= hi[idx1]) -
		log_sub2_exp(log_cp_lo, log_cp_hi)

	# Otherwise, work with the CDF function.
	log_p_hi = pf(hi[idx2], log.p = TRUE, ...)
	out[idx2] = df(x[idx2], log = TRUE, ...) +
		log(lo[idx2] < x[idx2] & x[idx2] <= hi[idx2]) -
		log_sub2_exp(log_p_hi, log_p_lo[idx2])

	if (log) { return(out) } else { return(exp(out)) }
}

#' @name TruncatedUnivariate
#' @export
p_truncated = function(x, lo, hi, pf, log.p = FALSE, tol = 1e-6, ...)
{
	n = length(x)
	out = numeric(n)
	if (length(lo) == 1) { lo = rep(lo, n) }
	if (length(hi) == 1) { hi = rep(hi, n) }

	log_p_lo = pf(lo, log.p = TRUE, ...)
	idx1 = which(log_p_lo > log1p(-tol))
	idx2 = setdiff(1:n, idx1)

	# If we are truncating way into the upper tail of a distribution, working
	# with the complement of the CDF helps to retain precision.
	log_cp_lo = pf(lo[idx1], log.p = TRUE, lower.tail = FALSE, ...)
	log_cp_hi = pf(hi[idx1], log.p = TRUE, lower.tail = FALSE, ...)
	log_cp_x = pf(x[idx1], log.p = TRUE, lower.tail = FALSE, ...)
	out[idx1] = log_sub2_exp(log_cp_lo, log_cp_x) - log_sub2_exp(log_cp_lo, log_cp_hi)

	# Otherwise, work with the CDF function.
	log_p_hi = pf(hi[idx2], log.p = TRUE, ...)
	log_p_x = pf(x[idx2], log.p = TRUE, ...)
	out[idx2] = log_sub2_exp(log_p_x, log_p_lo[idx2]) - log_sub2_exp(log_p_hi, log_p_lo)

	if (log.p) { return(out) } else { return(exp(out)) }
}

#' @name TruncatedUnivariate
#' @export
q_truncated = function(p, lo, hi, pf, qf, log.p = FALSE, tol = 1e-6, ...)
{
	if (log.p) { lp = p } else { lp = log(p) }

	n = length(p)
	out = numeric(n)
	if (length(lo) == 1) { lo = rep(lo, n) }
	if (length(hi) == 1) { hi = rep(hi, n) }

	log_p_lo = pf(lo, log.p = TRUE, ...)
	idx1 = which(log_p_lo > log1p(-tol))
	idx2 = setdiff(1:n, idx1)

	# If we are truncating way into the upper tail of a distribution, working
	# with the complement of the CDF helps to retain precision.
	log_cp_lo = pf(lo[idx1], log.p = TRUE, lower.tail = FALSE, ...)
	log_cp_hi = pf(hi[idx1], log.p = TRUE, lower.tail = FALSE, ...)
	log_p_adj = log_sub2_exp(log_cp_lo, lp[idx1] + log_sub2_exp(log_cp_lo, log_cp_hi))
	out[idx1] = qf(log_p_adj, log.p = TRUE, lower.tail = FALSE, ...)

	# Otherwise, work with the CDF function.
	log_p_hi = pf(hi[idx2], log.p = TRUE, ...)
	log_p_adj = log_add2_exp(log_p_lo[idx2], lp[idx2] + log_sub2_exp(log_p_hi, log_p_lo[idx2]))
	out[idx2] = qf(log_p_adj, log.p = TRUE, ...)

	return(out)
}

#' @name TruncatedUnivariate
#' @export
r_truncated = function(n, lo, hi, pf, qf, ...)
{
	u = runif(n)
	q_truncated(log(u), lo, hi, pf, qf, log.p = TRUE, ...)
}
