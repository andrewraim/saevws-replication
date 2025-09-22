#ifndef VWS_LOG_SUM_EXP_H
#define VWS_LOG_SUM_EXP_H

#include <Rcpp.h>

/*
* C++ implementations of log-sum-exp functions. These carry out addition and
* subtraction for arguments given on the log-scale, without explicitly
* exponentiating and taking logarithms. They can be used to avoid precision
* issues for numbers with very small or large magnitudes.
*/

namespace vws {

/*
* Returns `log(exp(x) + exp(y))` for scalars `x` and `y`.
*/
inline double log_add2_exp(double x, double y)
{
	double s = std::min(x,y);
	double t = std::max(x,y);
	return t + std::log1p(exp(s - t));
}

/*
* Returns `log(exp(x) - exp(y))` for scalars `x` and `y`. Result is `NaN` when
* $x < y$.
*/
inline double log_sub2_exp(double x, double y)
{
	if (std::isinf(x) && std::isinf(y) && x < 0 && y < 0) {
		return R_NegInf;
	}

	return x + std::log1p(-exp(y - x));
}

/*
* Returns elementwise `log(exp(x) + exp(y))` for $n$-dimensional vectors `x`
* and `y`. This version operates on STL vectors.
*/
inline std::vector<double> log_add2_exp(const std::vector<double>& x, const std::vector<double>& y)
{
	unsigned int n = x.size();
	if (n != y.size()) {
		Rcpp::stop("Dimensions do not match");
	}

	std::vector<double> out(n);
	for (unsigned int i = 0; i < n; i++) {
		out[i] =  log_add2_exp(x[i], y[i]);
	}

	return out;
}

/*
* Returns elementwise `log(exp(x) - exp(y))` for $n$-dimensional vectors `x`
* and `y`. This version operates on STL vectors. Result is `NaN` for $x < y$.
*/
inline std::vector<double> log_sub2_exp(const std::vector<double>& x, const std::vector<double>& y)
{
	unsigned int n = x.size();
	if (n != y.size()) {
		Rcpp::stop("Dimensions do not match");
	}

	std::vector<double> out(n);
	for (unsigned int i = 0; i < n; i++) {
		out[i] =  log_sub2_exp(x[i], y[i]);
	}

	return out;
}

/*
* Returns elementwise `log(exp(x) + exp(y))` for $n$-dimensional vectors `x`
* and `y`. This version operates on Rcpp vectors.
*/
inline Rcpp::NumericVector log_add2_exp(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y)
{
	unsigned int n = x.size();
	if (n != y.size()) {
		Rcpp::stop("Dimensions do not match");
	}

	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) =  log_add2_exp(x(i), y(i));
	}

	return out;
}

/*
* Returns elementwise `log(exp(x) - exp(y))` for $n$-dimensional vectors `x`
* and `y`. This version operates on Rcpp vectors. Result is `NaN` for $x < y$.
*/

inline Rcpp::NumericVector log_sub2_exp(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y)
{
	unsigned int n = x.size();
	if (n != y.size()) {
		Rcpp::stop("Dimensions do not match");
	}

	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) =  log_sub2_exp(x(i), y(i));
	}

	return out;
}

/*
* Returns scalar `log(sum(exp(x))` for vector `x`. Computed using the method
* described in StackExchange post <https://stats.stackexchange.com/a/381937>.
*/
inline double log_sum_exp(const Rcpp::NumericVector& x)
{
    unsigned int k = x.size();

    Rcpp::NumericVector v = Rcpp::clone(x);
    std::sort(v.begin(), v.end(), std::greater<double>());

    double s = v(0);

    for (unsigned int j = 1; j < k; j++) {
        double ind = s > R_NegInf;
        double dd = v(j) - std::pow(s, ind);
        s = std::max(v(j), s) + std::log1p(std::exp(-std::fabs(dd)));
    }

    return s;
}

/*
* Apply log-sum-exp operation to rows of a matrix
*
* - `x`: an $n \times k$ matrix
*
* Returns vector `log(sum(exp(x[i,]))`for rows $i = 1, \ldots, n$.
*/
inline Rcpp::NumericVector log_sum_exp_mat(const Rcpp::NumericMatrix& x)
{
    unsigned int n = x.nrow();
    Rcpp::NumericVector out(n);

    for (unsigned int i = 0; i < n; i++) {
        out(i) = log_sum_exp(x.row(i));
    }

    return out;
}

}

#endif
