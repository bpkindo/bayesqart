#ifndef GUARD_bd_h
#define GUARD_bd_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include "info.h"
#include "tree.h"
#include "rgig.h"


bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi,double GAMMA_QUANT,
            double TAU_SQ_QUANT, size_t min_obs_node, std::vector<int>& varused, 
			size_t maxdepth);
#endif
