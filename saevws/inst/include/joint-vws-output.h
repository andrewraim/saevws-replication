#ifndef SAEVWS_JOINT_VWS_OUTPUT_H
#define SAEVWS_JOINT_VWS_OUTPUT_H

#include <RcppArmadillo.h>

struct joint_vws_output {
	arma::vec sigma2;
	arma::uvec rejects;
	arma::uvec updates;
	arma::vec log_bound;
};

#endif
