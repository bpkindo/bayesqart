#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <RcppArmadillo.h>
#include <Rcpp.h>
//pi and log(2*pi)

#define ZTOL 0.00000000001

//compute gig mod
double _gig_mode(double lambda, double omega);
double _rgig_ROU_shift_alt (double lambda, double lambda_old, double omega, double alpha);
double  _rgig_ROU_noshift (double lambda, double lambda_old, double omega, double alpha);
double _rgig_newapproach1 (double lambda, double lambda_old, double omega, double alpha);

//generalized inverse gaussian generation
double do_rgig(double lambda, double chi, double psi);
/*   lambda .. parameter for gig distribution                                    */
/*   chi   ... parameter for gig  distribution                                    */
/*   psi   ... parameter for gig  distribution                                    */

