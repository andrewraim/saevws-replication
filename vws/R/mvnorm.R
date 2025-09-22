#' Multivariate Normal distribution
#'
#' k-dimensional Multivariate Normal distribution.
#'
#' @param n Number of observations.
#' @param x An \eqn{n \times k} matrix with n arguments.
#' @param mu mean vector parameter.
#' @param Sigma Covariance matrix parameter.
#' @param Sigma_chol Cholesky factor of covariance matrix parameter.
#' @param Omega Precision matrix parameter.
#' @param log If `TRUE`, return densities and probabilities on the log-scale.
#'
#' @return
#' \item{`dmvnorm`}{gives the density using specified covariance matrix,}
#' \item{`dmvnorm_prec`}{gives the density using specified precision matrix,}
#' \item{`rmvnorm`}{generates random deviates using specified covariance
#' matrix,}
#' \item{`rmvnorm_prec`}{generates random deviates using specified precision
#' matrix,}
#' \item{`r_singular_mvnorm`}{generates random deviates using specified singular
#' covariance matrix (and mean zero).}
#'
#' @details
#' A use case for the function `r_singular_mvnorm` would be to start with
#' a singular precision matrix `Q`, and use something like
#' `Sigma = ginverse(Q)` as the covariance matrix.
#'
#' @name MultivariateNormal
NULL

#' @name MultivariateNormal
#' @export
r_mvnorm = function(n, mu, Sigma)
{
	k = length(mu)
	stopifnot(k == nrow(Sigma) && k == ncol(Sigma))
	Z = matrix(rnorm(n*k), k, n)
	A = t(chol(Sigma))
	A %*% Z + mu %x% matrix(1, 1, n)
}

#' @name MultivariateNormal
#' @export
r_mvnorm_chol = function(n, mu, Sigma_chol)
{
	k = length(mu)
	stopifnot(k == nrow(Sigma_chol) && k == ncol(Sigma_chol))
	Z = matrix(rnorm(n*k), k, n)
	Sigma_chol %*% Z + mu %x% matrix(1, 1, n)
}

#' @name MultivariateNormal
#' @export
r_mvnorm_prec = function(n, mu, Omega)
{
	k = length(mu)
	stopifnot(k == nrow(Omega) && k == ncol(Omega))
	Z = matrix(rnorm(n*k), k, n)
	A = chol(Omega)
	solve(A, Z) + mu %x% matrix(1, 1, n)
}

#' @name MultivariateNormal
#' @export
d_mvnorm = function(x, mu, Sigma, log = FALSE)
{
	n = nrow(x)
	k = length(mu)
	stopifnot(k == nrow(Sigma) && k == ncol(Sigma) && k == ncol(x))

	# Use QR decomposition to efficiently compute determinant and inverse
	xc = x - t(mu) %x% matrix(1,n,1)
	qr_out = qr(Sigma)
	Omega_xc = solve(Sigma, t(xc))
	logdetA = sum(log(abs(diag(qr.R(qr_out)))))
	out = -k/2*log(2*pi) - logdetA / 2 - rowSums(xc * t(Omega_xc)) / 2

	if (log) { return(out) } else { return(exp(out)) }
}

#' @name MultivariateNormal
#' @export
d_mvnorm_prec = function(x, mu, Omega, log = FALSE)
{
	n = nrow(x)
	k = length(mu)
	stopifnot(k == nrow(Omega) && k == ncol(Omega) && k == ncol(x))

	xc = x - t(mu) %x% matrix(1,n,1)
	Omega_xc = Omega %*% t(xc)
	logdetA = -as.numeric(determinant(Omega)$modulus)
	out = -k/2*log(2*pi) - logdetA / 2 - rowSums(xc * t(Omega_xc)) / 2

	if (log) { return(out) } else { return(exp(out)) }
}

