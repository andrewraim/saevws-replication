#' @useDynLib saevws, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
#' @import stats
#' @import ggplot2
#' @importFrom vws printf
#' @importFrom mcmcse ess
#' @importFrom coda mcmc geweke.diag
NULL

utils::globalVariables(c("eff", "draw"))
