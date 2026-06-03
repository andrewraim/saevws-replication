#' Geweke diagnostic
#'
#' Compute the Geweke diagnostic with the coda package and (optionally) replace
#' `NA` values with `-Inf`.
#'
#' @param x TBD
#' @param frac1 TBD
#' @param frac2 TBD
#' @param na.rm logical; if `TRUE`, replace `NA` with `-Inf`.
#'
#' @returns A vector of Geweke z-values
#' @export
geweke = function(x, frac1 = 0.1, frac2 = 0.5, na.rm = TRUE)
{
	mcmc_out = mcmc(x)
	geweke_out = geweke.diag(mcmc_out, frac1, frac2)
	out = geweke_out$z
	out[is.na(out)] = -Inf
	return(out)
}
