#ifndef VWS_STEP_OUTPUT_H
#define VWS_STEP_OUTPUT_H

#include <RcppArmadillo.h>

struct VWSStepOutput {
	arma::vec sigma2;
	arma::uvec rejects;
	arma::uvec updates;
	arma::vec log_bound;
};

#endif
