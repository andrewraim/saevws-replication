#' @useDynLib saevws, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom stats sd step quantile
#' @importFrom ggplot2 ggplot aes annotate expansion geom_line
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_light
#' @importFrom ggplot2 theme_minimal xlab ylab
#' @importFrom vws printf
#' @importFrom mcmcse ess
#' @importFrom coda mcmc geweke.diag
#' @importFrom dplyr filter mutate %>% row_number
NULL

utils::globalVariables(c("eff", "draw", "iter", "count", "updates"))
