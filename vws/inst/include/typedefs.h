#ifndef VWS_TYPEDEFS_H
#define VWS_TYPEDEFS_H

#include <Rcpp.h>

namespace vws {

/*
* A constant that represents the maximum value of unsigned int.
*/
static const unsigned int uint_max = std::numeric_limits<unsigned int>::max();

/*
* Typedefs for several functions used in VWS programming.
*/
typedef std::function<double(double, bool)> uv_weight_function;
typedef std::function<double(const Rcpp::NumericVector&, bool)> mv_weight_function;

/*
* Enumeration that represents actions to be taken when an error condition
* arises. The general meanings are:
*
* - `STOP`: an exception should be thrown;
* - `WARNING`: a warning should be emitted but execution should continue;
* - `MESSAGE`: a message should be emitted but execution should continue;
* - `NONE`: the condition should be ignored.
*/
enum class error_action {
	STOP,
	WARNING,
	MESSAGE,
	NONE
};

}

#endif
