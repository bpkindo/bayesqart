// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <ctime>
#include "RcppArmadillo.h"
#include "info.h"
#include "tree.h"
#include "functions.h"
#include "rgig.h"

using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Rcpp::List BayesQArt(arma::vec const& y,
					arma::mat const& X, 
					arma::mat const& Xtest, 
					double quantile,
					int burn, 
					int nd, 
					int m,
					int min_obs_node, 
					double aa_parm,
					double bb_parm,
					int nc, 
					double pbd, 
					double pb, 
					double alpha, 
					double betap, 
					double kappa, 
					int maxdepth) {
					
  RNGScope scope;         // Initialize Random number generator
  size_t n = y.size();
  size_t p = X.size()/n;
  const size_t minperbot=(size_t)min_obs_node;	
   double GAMMA_QUANT = (1.0 - 2.0*quantile)/(quantile *(1-quantile));
   double TAU_SQ_QUANT = 2.0/(quantile *(1-quantile));
   

  Rcpp::Rcout <<" number of observations: " << n << std::endl;
  Rcpp::Rcout <<" number of covariates: " << p << std::endl;
  // scale min and max are -0.5 and 0.5
  arma::vec yscaled;
  //center and scale it
  arma::vec ysorted = sort(y);
  arma::vec ycentered = y - (ysorted(floor(quantile*n)) - GAMMA_QUANT);
  yscaled = -0.5 + (ycentered - arma::min(ycentered))/(arma::max(ycentered) - arma::min(ycentered));
  arma::vec v(n); //storage for latent vector v
  sinfo allys_vs;       //sufficient stats for all of y, use to initialize the trees.
  v.fill(1.0);
  allys_vs.sum_y = arma::as_scalar(arma::sum(yscaled));
  allys_vs.sum_v_inv = arma::as_scalar(arma::sum(1/v));
  allys_vs.sy2 = arma::as_scalar(yscaled.t() * yscaled);
  allys_vs.n = n;
  allys_vs.sum_y_div_v = arma::as_scalar(arma::sum(yscaled/v));
  
  //x for predictions
  dinfo dip; //data information for prediction
  dip.x = Xtest;
  dip.n=0;
  size_t np = dip.x.size()/p;

  if(np>0){
	dip.n=np; dip.p=p;
  }


   xinfo xi;
   makexinfo(p,n,X,xi,nc); 
   std::vector<tree> t(m);
 //for now initialize it at zero since I'm using the "centered" y, there could be a better value that  I need to figure out. Probably to depend on the quantile being calculated.
   for(size_t j=0;j<(size_t)m;j++) t[j].setm(0.0);
  //--------------------------------------------------
  //prior and mcmc
   pinfo pi;
   pi.pbd=pbd; //prob of birth/death move
   pi.pb=pb; //prob of birth given  birth/death

   pi.alpha=alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
   pi.betap=betap; //2 for bart means it is harder to build big trees.	
   pi.sigma0 = sqrt(0.5/(kappa*sqrt((double)m)));
   
   pi.phi = 1.0;
   arma::vec y_current(n);
   y_current.fill(0.0);
   arma::vec latent_current(n);   
   latent_current.fill(1.0);
   arma::vec ftemp(n);
   arma::vec allfit(n);
   allfit.fill(0.0);

   dinfo di;
   di.n=n; di.p=p; di.x = X; di.y=y_current;di.vv=latent_current;
   double rgig_chi, rgig_psi;

 //--------------------------------------------------
 //storage for ouput
 //in sample fit
  arma::vec pmean(n); //posterior mean of in-sample fit, sum draws,then divide
       //out of sample fit
  pmean.fill(0.0);
  arma::vec ppredmean(1); //posterior mean for prediction. If there is test data, resize latar.
  ppredmean.fill(0.0);
  arma::vec fpredtemp(1); //temporary fit vector to compute prediction. If there is test data, resize latar.
  fpredtemp.fill(0.0);
  if(dip.n>0) {
      ppredmean.resize(dip.n);
	  ppredmean.fill(0.00);
      fpredtemp.resize(dip.n);
	  fpredtemp.fill(0.00);
   }
   
   
   double cvptemp = 1.0/(double)di.p;
	for(size_t i=0;i<di.p;i++) {
	  std::vector<double> cvvtemp;
	  for(size_t j=0;j<di.p;j++) {
	        
	    cvvtemp.push_back(cvptemp);
	  }
	  pi.cvpm.push_back(cvvtemp);
	}

	
   
	/*for(size_t i=0;i<di.p;i++) {
		for(size_t j=0;j<di.p;j++)
			Rcpp:Rcout  << "(" << i << "," << j << ")" << pi.cvpm[i][j] << std::endl;
	}*/


//--------------------------------------------------
//mcmc

  
time_t tp;
int time1 = time(&tp);
arma::vec restemp;
std::vector<int> splitvars;
for(size_t i=0;i<(nd+burn);i++) {
if(i%100==0) Rcpp::Rcout << "mcmc iteration: " << i << std::endl;
		for(size_t j=0; j<(size_t)m; j++){
			rgig_psi = (2 + GAMMA_QUANT*GAMMA_QUANT/TAU_SQ_QUANT)/pi.phi;
			for(size_t k=0;k<n;k++) {
			rgig_chi = (yscaled(k) - allfit(k))*(yscaled(k) - allfit(k))/(pi.phi*TAU_SQ_QUANT);
			latent_current(k) = do_rgig(0.5,rgig_chi,rgig_psi);
		}
		di.vv=latent_current;
		fit(t[j],xi,di,ftemp);//current tree fit
		allfit -= ftemp;
		y_current = yscaled-allfit- GAMMA_QUANT*latent_current; //current y is the residual
		//sample latent to itegrate out

		di.y=y_current;
		
		if(unif_rand() > pi.pbd){
			tree::tree_p tnew;
			tnew=new tree(t[j]);
			calirotp(tnew, t[j], xi, di, pi, GAMMA_QUANT, TAU_SQ_QUANT, (size_t) min_obs_node);
			delete tnew;
		} else {
			bd(t[j],xi,di,pi,GAMMA_QUANT,TAU_SQ_QUANT, min_obs_node,splitvars,maxdepth);
		}
		
		//calichgv(t[j],xi,di,pi);
		calichgv(t[j], xi, di, pi, GAMMA_QUANT, TAU_SQ_QUANT, minperbot);

		drmu(t[j],xi,di,pi, GAMMA_QUANT,TAU_SQ_QUANT);
		fit(t[j],xi,di,ftemp);
		allfit += ftemp;
	}
	restemp = (yscaled-allfit - GAMMA_QUANT*latent_current);
	restemp %= (yscaled-allfit - GAMMA_QUANT*latent_current);
	restemp %= arma::pow(2*TAU_SQ_QUANT*latent_current, -1);
    pi.phi = 1/R::rgamma(n/2+aa_parm/2, (bb_parm + arma::sum(restemp))/2);
	

	if(i>=(size_t)burn) {
        pmean += allfit;
        if(dip.n) {
         for(size_t j=0;j<(size_t)m;j++) {
			fit(t[j],xi,dip,fpredtemp);
			ppredmean += fpredtemp;
          }
        }
	}
}
   int time2 = time(&tp);
   Rcpp::Rcout << "time for MCMC loop: " << time2-time1 << " seconds" << std::endl;
   pmean /= (double)nd;
   ppredmean /= (double)nd;
   
   arma::ivec varsused = arma::conv_to< arma::icolvec >::from(splitvars);
   
return Rcpp::List::create(
            Rcpp::Named("pred_train")=(ysorted(floor(quantile*n)) - GAMMA_QUANT) + arma::min(ycentered) + (pmean + 0.5)*(arma::max(ycentered) - arma::min(ycentered)), 
			Rcpp::Named("pred_test") = (ysorted(floor(quantile*n)) - GAMMA_QUANT) + arma::min(ycentered) + (ppredmean + 0.5)*(arma::max(ycentered) - arma::min(ycentered)), 
			Rcpp::Named("vars_used")=varsused
			);
}




