#ifndef VWS_LOGGER_H
#define VWS_LOGGER_H

#include <Rcpp.h>

namespace vws {

/*
* Logger
*
* Print a message prepended with a time stamp.
*
* - `fmt`: format string for message; should be compatible with `printf`.
* - `...`: variable argument list; is passed to `printf`
*
* Returns a status code from the function `vprintf`.
*/
inline int logger(const char* fmt, ...)
{
	// Prepare a string with the time stamp
	time_t timer;
	char buffer[64];
	struct tm* tm_info;

	timer = time(NULL);
	tm_info = localtime(&timer);

	strftime(buffer, 64, "%Y-%m-%d %H:%M:%S", tm_info);

	// Print time stamp and separator
	Rprintf("%s%s", buffer, " - ");

	// Call printf to handle the rest of the message
	va_list arg;
	va_start(arg, fmt);
	int out = vprintf(fmt, arg);
	va_end(arg);
	return out;
}

}

#endif
