#' Printing
#'
#' Functions to print messages using an `sprintf` syntax.
#'
#' @param fmt Format string which can be processed by `sprintf`
#' @param dt_fmt Format string which can be processed by `format.POSIXct`
#' @param file A connection, or a character string naming the file to print to
#' @param join A string to place between the timestamp and the message.
#' @param ... Additional arguments
#'
#' @examples
#' printf("Hello world %f %d\n", 0.1, 5)
#' logger("Hello world\n")
#' logger("Hello world %f %d\n", 0.1, 5)
#' logger("Hello world %f %d\n", 0.1, 5, dt_fmt = "%H:%M:%S")
#' logger("Hello world %f %d\n", 0.1, 5, join = " >> ")
#' logger("Hello world %f %d\n", 0.1, 5, join = " ")
#'
#' @name Print
NULL

#' @name Print
#' @export
printf = function(fmt, ...)
{
	cat(sprintf(fmt, ...))
}

#' @name Print
#' @export
logger = function(fmt, ..., dt_fmt = "%Y-%m-%d %H:%M:%S", join = " - ")
{
	sys.time = format(Sys.time(), format = dt_fmt)
	cat(sys.time, join, sprintf(fmt, ...), sep = "")
}

#' @name Print
#' @export
fprintf = function(file, fmt, ...)
{
	cat(sprintf(fmt, ...), file = file)
}

