#ifndef JOINT_VWS_OUTPUT_H
#define JOINT_VWS_OUTPUT_H

#include <RcppArmadillo.h>

struct JointVWSOutput {
	arma::vec sigma2;
	arma::uvec rejects;
	arma::uvec updates;
	arma::vec log_bound;
};

#endif
