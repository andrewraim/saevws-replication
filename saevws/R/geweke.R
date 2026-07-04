#' Geweke Diagnostic
#'
#' Compute the Geweke diagnostic with the coda package and (optionally) replace
#' `NA` values with `-Inf`.
#'
#' @param x Vector representing a series of draws
#' @param frac1 Proportion of draws to consider as the head of the sequence.
#' @param frac2 Proportion of draws to consider as the tail of the sequence.
#' @param na.rm logical; if `TRUE`, replace `NA` with `-Inf`.
#'
#' @returns A vector of Geweke z-values
#'
#' @export
geweke = function(x, frac1 = 0.1, frac2 = 0.5, na.rm = TRUE)
{
	mcmc_out = mcmc(x)
	geweke_out = geweke.diag(mcmc_out, frac1, frac2)
	out = geweke_out$z
	out[is.na(out)] = -Inf
	return(out)
}
