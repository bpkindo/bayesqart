#ifndef GUARD_info_h
#define GUARD_info_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

//data
class dinfo {
public:
   dinfo() {p=0;n=0;}
   size_t p;  //number of vars
   size_t n;  //number of observations
   arma::mat x; // jth var of ith obs is x(i,j)
   arma::vec y; // ith y  y[i]
   arma::vec vv; // ith v  v[i]
};

//prior and mcmc
class pinfo
{
public:
   pinfo() {pertalpha=0.10; pbd=1.0;pb=.5;alpha=.95;betap=.5;tau=1.0;phi=1.0;sigma0=1.0,taccept=0;tavgd=0;tmaxd=0;}
//mcmc info
	double pertalpha; //alpha for scaling pert proposal
   double pbd; //prob of birth/death
   double pb;  //prob of birth
//prior info
   double alpha;
   double betap;
   double tau;
//sigma
   double phi;
   double sigma0;
   std::vector< std::vector<double> > cvpm;  // Change of variable proposal probability matrix
   unsigned int taccept; //acceptance count of tree proposals
   unsigned int tproposal; //number of tree proposals
   unsigned int tavgd; //average tree depth
   unsigned int tmaxd; //maximum tree depth

};





//sufficient statistics for 1 node
class sinfo
{
public:
   sinfo() {n=0;sum_y=0.0;sum_y_div_v=0.0;sum_v_inv = 0.0;sy2=0.0;}
   size_t n;
   double sum_y;
   double sum_y_div_v;
   double sum_v_inv;
   double sy2;
};

#endif
