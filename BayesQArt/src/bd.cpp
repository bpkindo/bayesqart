#include <iostream>

#include "info.h"
#include "tree.h"
#include "bd.h"
#include "functions.h"


bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi,double GAMMA_QUANT,
            double TAU_SQ_QUANT, size_t min_obs_node, std::vector<int>& varused, 
			size_t maxdepth)
{

   tree::npv goodbots;  //nodes we could birth at (split on)
    double PBx = getpb(x,xi,pi,goodbots,maxdepth); //prob of a birth at x

     if(R::runif(0,1) < PBx) { //do birth or death
      //--------------------------------------------------
      //draw proposal

      //draw bottom node, choose node index ni from list in goodbots

      size_t ni = floor(R::runif(0,1)*goodbots.size());
      tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at

      //draw v,  the variable
      std::vector<size_t> goodvars; //variables nx can split on
      getgoodvars(nx,xi,goodvars);
      size_t vi = floor(R::runif(0,1)*goodvars.size()); //index of chosen split variable
      size_t v = goodvars[vi];

      //draw c, the cutpoint
      int L,U;
      L=0; U = xi[v].size()-1;
      nx->rg(v,&L,&U);

      size_t c = L + floor(R::runif(0,1)*(U-L+1)); //U-L+1 is number of available split points

//      Rcpp::Rcout << "chosen variable and cutpoint are: " << v << ", " << xi[v][c] << std::endl;

            //--------------------------------------------------
      //compute things needed for metropolis ratio

      double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
      size_t dnx = nx->depth();
      double PGnx = pi.alpha/pow(1.0 + dnx,pi.betap); //prior prob of growing at nx
      double PGly, PGry; //prior probs of growing at new children (l and r) of proposal

         if(goodvars.size()>1) { //know there are variables we could split l and r on
         PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.betap); //depth of new nodes would be one more
         PGry = PGly;
      } else { //only had v to work with, if it is exhausted at either child need PG=0
         if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
            PGly = 0.0;
         } else {
            PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.betap);
         }
         if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
            PGry = 0.0;
         } else {
            PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.betap);
         }
      }

      double PDy; //prob of proposing death at y
      if(goodbots.size()>1) { //can birth at y because splittable nodes left
         PDy = 1.0 - pi.pb;
      } else { //nx was the only node you could split on
         if((PGry==0) && (PGly==0)) { //cannot birth at y
            PDy=1.0;
         } else { //y can birth at either l or r
            PDy = 1.0 - pi.pb;
         }
      }

    double Pnogy; //death prob of choosing the nog node at y
      size_t nnogs = x.nnogs();
      tree::tree_cp nxp = nx->getp();
      if(nxp==0) { //no parent, nx is the top and only node
         Pnogy=1.0;
      } else {
         //if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
         if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
            Pnogy = 1.0/nnogs;
         } else { //if parent is not a nog, y has one more nog.
           Pnogy = 1.0/(nnogs+1.0);
         }
      }

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
      getsuff(x,nx,v,c,xi,di,sl,sr);
       //--------------------------------------------------
      //compute alpha
		
      double alpha=0.0,alpha1=0.0,alpha2=0.0;
      double lill=0.0,lilr=0.0,lilt=0.0;
  
	  if((sl.n>=min_obs_node) && (sr.n>=min_obs_node)) { 
	  	lill = lil(sl.sum_y_div_v, sl.sum_v_inv, sl.n, pi, GAMMA_QUANT, TAU_SQ_QUANT);  
		lilr = lil(sr.sum_y_div_v, sr.sum_v_inv, sr.n, pi, GAMMA_QUANT, TAU_SQ_QUANT); 	
		lilt = lil(sl.sum_y_div_v + sr.sum_y_div_v, 
					sl.sum_v_inv + sr.sum_v_inv, 
					sl.n + sr.n, pi, GAMMA_QUANT, TAU_SQ_QUANT);  
		
         alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
         alpha2 = alpha1*exp(lill+lilr-lilt);
         alpha = std::min(1.0,alpha2);
		 
      } 
	  else {
         alpha=0.0; 
        }
	

		
        double mul,mur,theta_star_left, theta_star_right, sigma_star_left, 
		sigma_star_right;
         if(R::runif(0,1) < alpha) {
		 
        //draw mean and variance to try metropolis hastings

        theta_star_left = pi.sigma0*pi.sigma0*(sl.sum_y_div_v)/
						 (sl.sum_v_inv*pi.sigma0*pi.sigma0 + TAU_SQ_QUANT*pi.phi);
		sigma_star_left = pi.sigma0*sqrt(TAU_SQ_QUANT*pi.phi)/
						sqrt(sl.sum_v_inv*pi.sigma0*pi.sigma0 + TAU_SQ_QUANT*pi.phi);

		theta_star_right = pi.sigma0*pi.sigma0*(sr.sum_y_div_v)/
						 (sr.sum_v_inv*pi.sigma0*pi.sigma0 + TAU_SQ_QUANT*pi.phi);
		sigma_star_right = pi.sigma0*sqrt(TAU_SQ_QUANT*pi.phi)/
						sqrt(sr.sum_v_inv*pi.sigma0*pi.sigma0 + TAU_SQ_QUANT*pi.phi);
         //means for new bottom nodes, left and right
        mul = theta_star_left + R::rnorm(0,1)*sigma_star_left;
        mur = theta_star_right + R::rnorm(0,1)*sigma_star_right;
        x.birth(nx->nid(),v,c,mul,mur);
        varused.push_back(x.getv());
        return true;
         } else {
        return false;
      }

	  } else {
//--------------------------------------------------
      //draw proposal
      //draw nog node, any nog node is a possibility
      tree::npv nognds; //nog nodes
      x.getnogs(nognds);
      size_t ni = floor(R::runif(0,1)*nognds.size());
      tree::tree_p nx = nognds[ni]; //the nog node we might kill children at
      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double PGny; //prob the nog node grows
      size_t dny = nx->depth();
      PGny = pi.alpha/pow(1.0+dny,pi.betap);

      //better way to code these two?
      double PGlx = pgrow(nx->getl(),xi,pi);
      double PGrx = pgrow(nx->getr(),xi,pi);

      double PBy;  //prob of birth move at y
        if(!(nx->p)) { //is the nog node nx the top node
         PBy = 1.0;
      } else {
         PBy = pi.pb;
      }
          double Pboty;  //prob of choosing the nog as bot to split on when y
      int ngood = goodbots.size();
      if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
      if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
      ++ngood;  //know you can split at nx
      Pboty=1.0/ngood;

      double PDx = 1.0-PBx; //prob of a death step at x
      double Pnogx = 1.0/nognds.size();
       //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl,sr; 

        getsuff(x,nx->getl(),nx->getr(),xi,di,sl,sr);
          //compute alpha
		  
		double lill = lil(sl.sum_y_div_v, sl.sum_v_inv, sl.n, pi, GAMMA_QUANT, TAU_SQ_QUANT);  
		double lilr = lil(sr.sum_y_div_v, sr.sum_v_inv, sr.n, pi, GAMMA_QUANT, TAU_SQ_QUANT); 	
		double lilt = lil(sl.sum_y_div_v + sr.sum_y_div_v, 
					sl.sum_v_inv + sr.sum_v_inv, 
					sl.n + sr.n, pi, GAMMA_QUANT, TAU_SQ_QUANT);  

      double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
      double alpha2 = alpha1*exp(lilt - lill - lilr);
      double alpha = std::min(1.0,alpha2);

    //metrop
        double mu, theta_star_nog, sum_y_div_v_nog, sum_v_inv_nog, sigma_star_nog;
        size_t n;


         if(R::runif(0,1)<alpha) {
        n = sl.n + sr.n;
        sum_y_div_v_nog = sr.sum_y_div_v + sl.sum_y_div_v;
        sum_v_inv_nog = sr.sum_v_inv + sl.sum_v_inv;
		
		theta_star_nog = pi.sigma0*pi.sigma0*(sum_y_div_v_nog)/
						 (sum_v_inv_nog*pi.sigma0*pi.sigma0 + TAU_SQ_QUANT*pi.phi);
		sigma_star_nog = pi.sigma0*sqrt(TAU_SQ_QUANT*pi.phi)/
						sqrt(sum_v_inv_nog*pi.sigma0*pi.sigma0 + TAU_SQ_QUANT*pi.phi);

         //means for new bottom nodes, left and right
        mu = theta_star_nog + R::rnorm(0,1)*sigma_star_nog;


        //death
        x.death(nx->nid(),mu);
        //Rcpp::Rcout << "death"<< std::endl;
        return true;
        } else {
        return false;
        }

}
}

     


