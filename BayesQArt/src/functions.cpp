#include <cmath>
#include <R.h>
#include <Rmath.h>
#include "RcppArmadillo.h"
#include "functions.h"
#include <map>


//function to draw truncated normal
//above=1 means from above b=trunpt, a=-inf
//above=0 means from below a=trunpt, b= +inf   
double rtrun(double mu, double sigma,double trunpt, int above){

	double FA,FB,rnd,result,arg;
	if (above) {
		FA = 0.0; FB = R::pnorm(((trunpt-mu)/(sigma)),0.0,1.0,1,0);
	} else {
		FB = 1.0; FA = R::pnorm(((trunpt-mu)/(sigma)),0.0,1.0,1,0);
	}
	
  rnd = unif_rand(); 
	arg = rnd*(FB-FA)+FA;
	if(arg > .999999999) arg = .999999999;
	if(arg < .0000000001) arg = .0000000001;
	result = mu + sigma*R::qnorm(arg,0.0,1.0,1,0);

	return (result);
}

//--------------------------------------------------
//make xinfo = cutpoints
// [[Rcpp::depends(RcppArmadillo)]]
void makexinfo(arma::uword p, arma::uword n, arma::mat const& x, xinfo& xi, arma::uword nc)
{
	double xinc;
  arma::vec minx(p);
  arma::vec maxx(p);
  for(arma::uword j = 0 ; j< p; ++j){
	minx(j) =  min(x.col(j));
	maxx(j) =  max(x.col(j));
  }
	
 //make grid of nc cutpoints between min and max for each x.
  xi.resize(p);

for(size_t i=0;i<p;i++) {
      xinc = (maxx(i)-minx(i))/(nc+1.0);
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = minx(i) + (j+1)*xinc;
}
}
//--------------------------------------------------
//fit double array fv
void fit(tree& t, xinfo& xi, dinfo& di, arma::vec& fv)
{
   tree::tree_cp bn;
   for(size_t i=0;i<di.n;++i) {
       bn = t.bn(i,di,xi);
	  fv(i) = bn->getm();
   }
}

//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots, size_t maxdepth)
{
   double pb;  //prob of birth to be returned
   tree::npv bnv; //all the bottom nodes
   t.getbots(bnv);
   for(size_t i=0;i!=bnv.size();i++){
	tree::tree_p nptr = bnv[i];
	size_t dpth = nptr->depth();
	if(cansplit(bnv[i],xi) && dpth <= maxdepth) goodbots.push_back(bnv[i]);
   }
   if(goodbots.size()==0) { //are there any bottom nodes you can split on?
      pb=0.0;
   } else {
      if(t.treesize()==1) pb=1.0; //is there just one node?
      else pb=pi.pb;
   }
   
   
   return pb;
}


//--------------------------------------------------
//does this bottom node n have any variables it can split on.
bool cansplit(tree::tree_p n, xinfo& xi)
{
   int L,U;
   bool v_found = false; //have you found a variable you can split on
   size_t v=0;
   while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) v_found=true;
      v++;
   }
   return v_found;
}//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
   int L,U;
   for(size_t v=0;v!=xi.size();v++) {//try each variable
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) goodvars.push_back(v);
   }
}

//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
void getsuff(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
   sl.n = 0.0; sl.sum_v_inv = 0.0; sl.sum_y = 0.0; sl.sum_y_div_v= 0.0;
   sr.n = 0.0; sr.sum_v_inv = 0.0; sr.sum_y = 0.0; sr.sum_y_div_v= 0.0;

   for(size_t i=0;i<di.n;i++) {
      if(nx==x.bn(i, di,xi)) { //does the bottom node = xx's bottom node
	   if(di.x(i,v) < xi[v][c]) {
               sl.n++;
               sl.sum_y += di.y(i);
               sl.sum_v_inv += (double)(1.0/di.vv(i));
               sl.sum_y_div_v += (double)(di.y(i)/di.vv(i));
          } else {
               sr.n++;
               sr.sum_y += di.y(i);
               sr.sum_v_inv += (double)(1.0/di.vv(i));
               sr.sum_y_div_v += (double)(di.y(i)/di.vv(i));
          }
      }
   }
   

}

//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
void getsuff(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
   sl.n = 0.0; sl.sum_v_inv = 0.0; sl.sum_y = 0.0; sl.sum_y_div_v= 0.0;
   sr.n = 0.0; sr.sum_v_inv = 0.0; sr.sum_y = 0.0; sr.sum_y_div_v= 0.0;

   for(size_t i=0;i<di.n;i++) {
       tree::tree_cp bn = x.bn(i, di,xi);
        if(bn==nl) {
               sl.n++;
               sl.sum_y += di.y(i);
               sl.sum_v_inv += (double)(1.0/di.vv(i));
               sl.sum_y_div_v += (double)(di.y(i)/di.vv(i));
        }
        if(bn==nr) {
               sr.n++;
               sr.sum_y += di.y(i);
               sr.sum_v_inv += (double)(1.0/di.vv(i));
               sr.sum_y_div_v += (double)(di.y(i)/di.vv(i));
          }
      }
}

//--------------------------------------------------
// similar, but in case of perturb.
void getpertsuff(tree& x, tree::tree_p pertnode, tree::npv& bnv, size_t oldc, xinfo& xi, dinfo& di, std::vector<sinfo>& svold, std::vector<sinfo>& svnew)
{
	allsuff(x,xi,di,bnv,svnew);
	pertnode->setc(oldc);
	allsuff(x,xi,di,bnv,svold);		
}
void getpertsuff(tree& x, tree::tree_p pertnode, tree::npv& bnv, size_t oldc, size_t oldvar, bool didswaplr, xinfo& xi, dinfo& di, std::vector<sinfo>& svold, std::vector<sinfo>& svnew)
{
	allsuff(x,xi,di,bnv,svnew);
	pertnode->setv(oldvar);
	pertnode->setc(oldc);
	if(didswaplr) //undo swap
	{
		tree::tree_p tempt;
		tempt=pertnode->r;
		pertnode->r=pertnode->l;
		pertnode->l=tempt;
	}
	allsuff(x,xi,di,bnv,svold);		
}


//--------------------------------------------------
//log of the integrated likelihood
double lil(double sum_y_div_v, double sum_v_inv, size_t n, pinfo& pi,double GAMMA_QUANT,double TAU_SQ_QUANT)
{
   double val1, val2_num, val2_den, rv,sig2;
   sig2 = pi.sigma0*pi.sigma0;
   val1 = 0.5* (log(sig2*TAU_SQ_QUANT*pi.phi) - log(TAU_SQ_QUANT*pi.phi + sig2*sum_v_inv));
   val2_num = sig2*pow(sum_y_div_v, 2);
   val2_den = 2.0*TAU_SQ_QUANT*pi.phi*(TAU_SQ_QUANT*pi.phi + sum_v_inv*sig2);
   rv = val1 + val2_num/val2_den;
   
   return rv;
}

//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi)
{
   if(cansplit(n,xi)) {
      return pi.alpha/pow(1.0+n->depth(),pi.betap);
   } else {
      return 0.0;
   }
}

//--------------------------------------------------
// draw all the bottom node mu's
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double GAMMA_QUANT,
            double TAU_SQ_QUANT)
{
   tree::npv bnv;
   std::vector<sinfo> sv;
   allsuff(t,xi,di,bnv,sv);
   double theta_star, sigma_star, mu_star;
   for(tree::npv::size_type i=0;i!=bnv.size();i++) {


		theta_star = pi.sigma0*pi.sigma0*(sv[i].sum_y_div_v)/
						 (sv[i].sum_v_inv*pi.sigma0*pi.sigma0 + TAU_SQ_QUANT*pi.phi);
		sigma_star = pi.sigma0*sqrt(TAU_SQ_QUANT*pi.phi)/
						sqrt(sv[i].sum_v_inv*pi.sigma0*pi.sigma0 + TAU_SQ_QUANT*pi.phi);
	   
       mu_star = theta_star + R::rnorm(0,1)*sigma_star;
       bnv[i]->setm(mu_star);


   }
}

//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv)
{
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   bnv.clear();
   x.getbots(bnv);
   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   sv.resize(nb);
   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;
   for(size_t i=0;i<di.n;i++) {
      tbn = x.bn(i,di,xi);
      ni = bnmap[tbn];

      ++(sv[ni].n);
      sv[ni].sum_y += di.y(i);
      sv[ni].sum_v_inv += (double)(1.0/di.vv(i));
      sv[ni].sum_y_div_v += (double)(di.y(i)/di.vv(i));

   }

}


//--------------------------------------------------
// Get the L,U values for a node in the tree *given* the tree
// structure both above and below that node.
void getLU(tree::tree_p pertnode, xinfo& xi, int* L, int* U)
{
	tree::tree_p l,r;

	*L=0; *U = xi[pertnode->getv()].size()-1;
	l=pertnode->getl();
	r=pertnode->getr();

	bool usel,user;
	usel=l->nuse(pertnode->getv());
	user=r->nuse(pertnode->getv());
	if(usel && user)
	{
		l->rl(pertnode->getv(),L);
		r->ru(pertnode->getv(),U);
	}
	else if(usel)
	{
		pertnode->rg(pertnode->getv(),L,U);
		l->rl(pertnode->getv(),L);
	}
	else
	{
		pertnode->rg(pertnode->getv(),L,U);
		r->ru(pertnode->getv(),U);
	}
}

//--------------------------------------------------
// similar except we get it for a prescribed variable pertvar
void getpertLU(tree::tree_p pertnode, size_t pertvar, xinfo& xi, int* L, int* U)
{
		*L=0; *U = xi[pertvar].size()-1;

		bool usel,user;
		usel=pertnode->l->nuse(pertvar);
		user=pertnode->r->nuse(pertvar);
		if(usel && user)
		{
			pertnode->l->rl(pertvar,L);
			pertnode->r->ru(pertvar,U);
		}
		else if(usel)
		{
			pertnode->rg(pertvar,L,U);
			pertnode->l->rl(pertvar,L);
		}
		else
		{
			pertnode->rg(pertvar,L,U);
			pertnode->r->ru(pertvar,U);
		}
}

//--------------------------------------------------
// number of available cutpoints at node n for variable var
int getnumcuts(tree::tree_p n, xinfo& xi, size_t var)
{
	int L,U;

	getpertLU(n,var,xi,&L,&U);
	return std::max(0,U-L+1);
}

//--------------------------------------------------
// Find numbr of variables internal tree node n can split on
void getinternalvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
   int L,U;

   for(size_t v=0;v!=xi.size();v++) {//try each variable
      L=0; U = xi[v].size()-1;
      getpertLU(n,v,xi,&L,&U);
      if(U>=L) goodvars.push_back(v);
   }
}



//--------------------------------------------------
// Rotate a given rotatable node in the tree.
// node n must be the left child of n->p and also not a root or terminal leaf node.
void rotright(tree::tree_p n)
{
	tree::tree_cp ctstar;
	tree::tree_p tstar;
	tree::tree_p newt = new tree;
	size_t vprime;
	size_t cprime;

	tstar=n->p->r;
	ctstar=n->p->r;
	tree::tree_p newnr = new tree(*ctstar);
	tstar->p=0;

	vprime=n->v;
	cprime=n->c;	
	n->v=n->p->v;
	n->c=n->p->c;
	n->p->v=vprime;
	n->p->c=cprime;
	newt->v=n->v;
	newt->c=n->c;

	newt->p=n->p;
	newt->l=n->r;
	newt->l->p=newt;
	newt->r=tstar;
	tstar->p=newt;
	n->p->r=newt;
	
	n->r=newnr;
	newnr->p=n;
}

void rotleft(tree::tree_p n)
{
	tree::tree_cp ctstar;
	tree::tree_p tstar;
	tree::tree_p newt = new tree;
	size_t vprime;
	size_t cprime;

	tstar=n->p->l;
	ctstar=n->p->l;
	tree::tree_p newnl = new tree(*ctstar);
	tstar->p=0;

	vprime=n->v;
	cprime=n->c;
	n->v=n->p->v;
	n->c=n->p->c;
	n->p->v=vprime;
	n->p->c=cprime;
	newt->v=n->v;
	newt->c=n->c;

	newt->p=n->p;
	newt->r=n->l;
	newt->r->p=newt;
	newt->l=tstar;
	tstar->p=newt;
	n->p->l=newt;

	n->l=newnl;
	newnl->p=n;
}

// reduce the left sub-tree of the node that was rotated to the top
void reduceleft(tree::tree_p n, size_t v, size_t c)
{
	tree::tree_p temp;

	if(n->r->l && n->r->v==v) //right is not terminal and splits on v
		if(n->r->c >= c) { //then only keep left branch
			delete n->r->r;
			temp=n->r;
			n->r=temp->l;
			temp->l->p=n;
			temp->r=0;
			temp->l=0;
			temp->p=0;
			delete temp;
		}
	if(n->l->l && n->l->v==v) //left is not terminal and splits on v
		if(n->l->c >= c) { // then only keep left branch
			delete n->l->r;
			temp=n->l;
			n->l=temp->l;
			temp->l->p=n;
			temp->r=0;
			temp->l=0;
			temp->p=0;
			delete temp;
		}
}

void reduceright(tree::tree_p n, size_t v, size_t c)
{
	tree::tree_p temp;

	if(n->r->v==v && n->r->l) //right is not terminal and splits on v
		if(n->r->c <= c) { //then only keep right branch
			delete n->r->l;
			temp=n->r;
			n->r=temp->r;
			temp->r->p=n;
			temp->r=0;
			temp->l=0;
			temp->p=0;
			delete temp;
		}
	if(n->l->v==v && n->l->l) //left is not terminal and splits on v
		if(n->l->c <= c) { // then only keep right branch
			delete n->l->l;
			temp=n->l;
			n->l=temp->r;
			temp->r->p=n;
			temp->r=0;
			temp->l=0;
			temp->p=0;
			delete temp;
		}
}

void splitleft(tree::tree_p t, size_t v, size_t c)
{
	tree::tree_p temp;

	if(t->l) //not terminal node
	{
		if(t->v==v && t->c >= c)
		{
			temp=t->l;
			if(t->isleft())
			{
				t->p->l=temp;
				temp->p=t->p;
			}
			else //isright
			{
				t->p->r=temp;
				temp->p=t->p;
			}
			delete t->r;
			t->p=0;
			t->r=0;
			t->l=0;
			delete t;
			t=temp;
			splitleft(t,v,c);
		}
		else
		{
			splitleft(t->l,v,c);
			splitleft(t->r,v,c);
		}
	}
}

void splitright(tree::tree_p t, size_t v, size_t c)
{
	tree::tree_p temp;

	if(t->l) //not terminal node
	{
		if(t->v==v && t->c <= c)
		{
			temp=t->r;
			if(t->isleft())
			{
				t->p->l=temp;
				temp->p=t->p;
			}
			else //isright
			{
				t->p->r=temp;
				temp->p=t->p;
			}
			delete t->l;
			t->p=0;
			t->l=0;
			t->r=0;
			delete t;
			t=temp;
			splitright(t,v,c);
		}
		else
		{
			splitright(t->l,v,c);
			splitright(t->r,v,c);
		}
	}
}


bool merge2(tree::tree_p tl, tree::tree_p tr, tree::tree_p t, size_t v, size_t c)
{
	bool m1,m2;
	tree::tree_cp temptl,temptr;
	int tnwl=0,tnwr=0;
	double u;

	u=unif_rand();

	if(arenodesleafs(tl,tr)) {  //merging type 3
		if(u<0.5) {
			t->v=tl->v;
			t->c=tl->c;
			t->mu=tl->mu;  //doesn't matter actually it will be overwritten.
			t->l=0;
			t->r=0;
		}
		else 
		{
			temptl=tl;
			temptr=tr;
			t->v=v;
			t->c=c;
			t->l=new tree(*temptl);
			t->r=new tree(*temptr);
			t->l->p=t;
			t->r->p=t;
		}
		return true;
	}
	else if(arenodesequal(tl,tr) && !splitsonv(tl,tr,v)) {  //merging type 4
		m1=merge3(tl->l,tr->l,v,c,&tnwl);
		m2=merge3(tl->r,tr->r,v,c,&tnwr);
		if(u < (1.0/(tnwl+tnwr+1.0)) )
		{
			temptl=tl;
			temptr=tr;
			t->v=v;
			t->c=c;
			t->l=new tree(*temptl);
			t->r=new tree(*temptr);
			t->l->p=t;
			t->r->p=t;
		}
		else
		{
			t->v=tl->v;
			t->c=tl->c;
			t->l=new tree;
			t->r=new tree;
			t->l->p=t;
			t->r->p=t;
			tnwl=0;
			tnwr=0;
			m1=merge2(tl->l,tr->l,t->l,v,c);
			m2=merge2(tl->r,tr->r,t->r,v,c);
		}
		return (m1 & m2);
	}
	else if(splitsonv(tl,tr,v)) {  //merging type 7
		m1=merge3(tl->r,tr,v,c,&tnwr);
		m2=merge3(tl,tr->l,v,c,&tnwl);
		if(u < (1.0/(tnwr+tnwl+1.0)) )
		{
			temptl=tl;
			temptr=tr;
			t->v=v;
			t->c=c;
			t->l=new tree(*temptl);
			t->r=new tree(*temptr);
			t->l->p=t;
			t->r->p=t;
		}
		else if(u < ((1.0+tnwr)/(1.0+tnwr+tnwl)) )
		{
			temptl=tl->l;
			t->v=tl->v;
			t->c=tl->c;
			t->l=new tree(*temptl);
			t->l->p=t;
			t->r=new tree;
			t->r->p=t;
			m2=merge2(tl->r,tr,t->r,v,c);
		}
		else
		{
			temptr=tr->r;
			t->v=tr->v;
			t->c=tr->c;
			t->r=new tree(*temptr);
			t->r->p=t;
			t->l=new tree;
			t->l->p=t;
			m1=merge2(tl,tr->l,t->l,v,c);
		}
		if(!m1) Rcpp::Rcout << "doh7a" << std::endl;
		if(!m2) Rcpp::Rcout << "doh7b" << std::endl;	
		return (m1 & m2);
	}
	else if(splitsonv(tl,v) && isleaf(tr)) //merging type 1
	{
		m1=merge3(tl->r,tr,v,c,&tnwr);
		if(u < (1.0/(tnwr + 1.0)) )
		{
			temptl=tl;
			temptr=tr;
			t->v=v;
			t->c=c;
			t->l=new tree(*temptl);
			t->r=new tree(*temptr);
			t->l->p=t;
			t->r->p=t;
		}
		else
		{
			temptl=tl->l;
			t->v=tl->v;
			t->c=tl->c;
			t->l=new tree(*temptl);
			t->l->p=t;
			t->r=new tree;
			t->r->p=t;
			m1=merge2(tl->r,tr,t->r,v,c);
		}
		if(!m1) Rcpp::Rcout << "doh1(m1)" << std::endl;
		return m1;
	}
	else if(splitsonv(tr,v) && isleaf(tl)) //merging type 2
	{
		m2=merge3(tl,tr->l,v,c,&tnwl);
		if(u < (1.0/(tnwl+1.0)) )
		{
			temptl=tl;
			temptr=tr;
			t->v=v;
			t->c=c;
			t->l=new tree(*temptl);
			t->r=new tree(*temptr);
			t->l->p=t;
			t->r->p=t;			
		}
		else
		{
			temptr=tr->r;
			t->v=tr->v;
			t->c=tr->c;
			t->r=new tree(*temptr);
			t->r->p=t;
			t->l=new tree;
			t->l->p=t;
			m2=merge2(tl,tr->l,t->l,v,c);
		}
		if(!m2) Rcpp::Rcout << "doh2(m2)" << std::endl;
		return m2;
	}
	else if(!isleaf(tl) && !isleaf(tr) && splitsonv(tr,v)) { //merge type 6(i)
		m1=merge3(tl,tr->l,v,c,&tnwr);
		if(u < (1.0/(tnwr + 1.0)) )
		{
			temptl=tl;
			temptr=tr;
			t->v=v;
			t->c=c;
			t->l=new tree(*temptl);
			t->r=new tree(*temptr);
			t->l->p=t;
			t->r->p=t;			
		}
		else
		{
			temptr=tr->r;
			t->v=tr->v;
			t->c=tr->c;
			t->r=new tree(*temptr);
			t->r->p=t;
			t->l=new tree;
			t->l->p=t;
			m1=merge2(tl,tr->l,t->l,v,c);
		}
		if(!m1) Rcpp::Rcout << "doh6i(m1)" << std::endl;
		return m1;
	}
	else if(!isleaf(tl) && !isleaf(tr) && splitsonv(tl,v)) { //merge type 6(ii)
		m2=merge3(tl->r,tr,v,c,&tnwl);
		if(u < (1.0/(tnwl + 1.0)) )
		{
			temptl=tl;
			temptr=tr;
			t->v=v;
			t->c=c;
			t->l=new tree(*temptl);
			t->r=new tree(*temptr);
			t->l->p=t;
			t->r->p=t;
		}
		else
		{
			temptl=tl->l;
			t->v=tl->v;
			t->c=tl->c;
			t->l=new tree(*temptl);
			t->l->p=t;
			t->r=new tree;
			t->r->p=t;
			m2=merge2(tl->r,tr,t->r,v,c);
		}
		if(!m2) Rcpp::Rcout << "doh6ii(m2)" << std::endl;
		return m2;
	}
	else if(!splitsonv(tl,v) && isleaf(tr)) { //merge type 5(i)
		temptl=tl;
		temptr=tr;
		t->v=v;
		t->c=c;
		t->l=new tree(*temptl);
		t->r=new tree(*temptr);
		t->l->p=t;
		t->r->p=t;
		return true;
	}
	else if(!splitsonv(tr,v) && isleaf(tl)) { //merge type 5(ii)
		temptl=tl;
		temptr=tr;
		t->v=v;
		t->c=c;
		t->l=new tree(*temptl);
		t->r=new tree(*temptr);
		t->l->p=t;
		t->r->p=t;
		return true;
	}
	else // default type aka type 8
	{
		temptl=tl;
		temptr=tr;
		t->v=v;
		t->c=c;
		t->l=new tree(*temptl);
		t->r=new tree(*temptr);
		t->l->p=t;
		t->r->p=t;
		return true;
	}

	return false;
}

// only to get nways, not to actually do the merge.
bool merge3(tree::tree_p tl, tree::tree_p tr, size_t v, size_t c, int* nways)
{
	bool m1,m2;
	int tnwl=0,tnwr=0;

	if(arenodesleafs(tl,tr)) {  //merging type 3
		*nways += 2;
		return true;
	}
	else if(arenodesequal(tl,tr) && !splitsonv(tl,tr,v)) {  //merging type 4
		*nways += 1;
		m1=merge3(tl->l,tr->l,v,c,&tnwl);
		m2=merge3(tl->r,tr->r,v,c,&tnwr);
		*nways += (tnwr*tnwl);
		return (m1 & m2);
	}
	else if(splitsonv(tl,tr,v)) {  //merging type 7
		*nways += 1;  //for this one
		m1=merge3(tl->r,tr,v,c,&tnwr);
		m2=merge3(tl,tr->l,v,c,&tnwl);
		*nways+= (tnwr+tnwl);
		if(!m1) Rcpp::Rcout << "doh7a" << std::endl;
		if(!m2) Rcpp::Rcout << "doh7b" << std::endl;	
		return (m1 & m2);
	}
	else if(splitsonv(tl,v) && isleaf(tr)) //merging type 1
	{
		*nways += 1; //for this one
		m1=merge3(tl->r,tr,v,c,&tnwr);
		*nways += tnwr;
		if(!m1) Rcpp::Rcout << "doh1(m1)" << std::endl;
		return m1;
	}
	else if(splitsonv(tr,v) && isleaf(tl)) //merging type 2
	{
		*nways += 1; //for this one
		m2=merge3(tl,tr->l,v,c,&tnwl);
		*nways += tnwl;
		if(!m2) Rcpp::Rcout << "doh2(m2)" << std::endl;
		return m2;
	}
	else if(!isleaf(tl) && !isleaf(tr) && splitsonv(tr,v)) { //merge type 6(i)
		*nways += 1; //for this one
		m1=merge3(tl,tr->l,v,c,&tnwr);
		*nways += tnwr;
		if(!m1) Rcpp::Rcout << "doh6i(m1)" << std::endl;
		return m1;
	}
	else if(!isleaf(tl) && !isleaf(tr) && splitsonv(tl,v)) { //merge type 6(ii)
		*nways +=1 ; //for this one
		m2=merge3(tl->r,tr,v,c,&tnwl);
		*nways += tnwl;
		if(!m2) Rcpp::Rcout << "doh6ii(m2)" << std::endl;
		return m2;
	}
	else if(!splitsonv(tl,v) && isleaf(tr)) { //merge type 5(i)
		*nways += 1; //for this one
		return true;
	}
	else if(!splitsonv(tr,v) && isleaf(tl)) { //merge type 5(ii)
		*nways += 1; //for this one
		return true;
	}
	else // default type aka type 8
	{
		*nways += 1; //for this one
		return true;
	}

	return false;
}

//--------------------------------------------------
// These ones support the rotate code, but could be generally useful too.
bool hasvcsplit(tree::tree_p t, size_t v, size_t c)
{
   tree::npv tnodes;

   tnodes.clear();
   t->getnodesonvc(tnodes,v,c);
   if(tnodes.size())
      return true;
   return false;
}

bool splitsonv(tree::tree_p t, size_t v)
{
   if(!t->getl()) //its a leaf
      return false;
   if(t->getv()==v) return true;
   return false;
}

bool splitsonv(tree::tree_p nl, tree::tree_p nr, size_t v)
{
   if(!nl->getl() || !nr->getl()) //one is a leaf
      return false;
   if(nl->getv()==v && nr->getv()==v) return true;
   return false;
}

bool isleaf(tree::tree_p t)
{
   if(!t->getl()) return true;
   return false;
}

bool arenodesequal(tree::tree_p nl, tree::tree_p nr)
{
   if(!nl->getl() || !nr->getl())  //one is a leaf
      return false;
   if(nl->getv()==nr->getv() && nl->getc()==nr->getc()) //both split on v at c
      return true;
   return false;
}

bool arenodesleafs(tree::tree_p nl, tree::tree_p nr)
{
   if(!nl->getl() && !nr->getl()) //both are leaves
      return true;
   return false;
}

bool calirotp(tree::tree_p tnew, tree& x, xinfo& xi, dinfo& di, pinfo& pi, double GAMMA_QUANT, double TAU_SQ_QUANT, size_t minperbot)
{
	tree::tree_p rotp,temp;
	tree::tree_cp xp;
	tree::npv subtold, subtnew, nbold, nbnew;
	double Qold_to_new, Qnew_to_old;
	unsigned int rdx=0;
	bool twowaystoinvert=false;
	double prinew=1.0,priold=1.0;
	size_t rotid;
	bool hardreject=false;
	std::vector<size_t> goodvars; //variables an internal node can split on

	pi.tproposal++;

	// Get rot nodes
	tree::npv rnodes;
	tnew->getrotnodes(rnodes);
	if(rnodes.size()==0) return false;  //no rot nodes so that's a reject.

	rdx = (unsigned int)floor(unif_rand()*rnodes.size()); //which rotatable node will we rotate at?
	rotp = rnodes[rdx];
	rotid=rotp->nid();
	xp=x.getptr(rotid);

	int nwaysm1=0,nwaysm2=0,nwayss1=0,nwayss2=0;
	double pm1=1.0,pm2=1.0,ps1=1.0,ps2=1.0;
	if(rotp->isleft()) {
		if(rotp->v==rotp->p->v) //special case, faster to handle it direclty
		{
			rotright(rotp);
			rotp=tnew->getptr(rotid);
			delete rotp->r;
			temp=rotp->l;
			rotp->p->l=temp;
			temp->p=rotp->p;
			rotp->r=0;
			rotp->l=0;
			rotp->p=0;
			delete rotp;
			rotp=tnew->getptr(rotid);
			//pm1=pm2=ps1=ps2=1.0 in this case
		}
		else
		{
			rotright(rotp);
			rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
			reduceleft(rotp,rotp->p->v,rotp->p->c);
			rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
			reduceright(rotp->p->r,rotp->p->v,rotp->p->c);
			rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
			splitleft(rotp->r,rotp->p->v,rotp->p->c);
			splitright(rotp->p->r->r,rotp->p->v,rotp->p->c);

			merge3(rotp->r,rotp->p->r->r,rotp->p->v,rotp->p->c,&nwayss1);
			ps1=1.0/nwayss1;

			merge3(rotp->l,rotp->p->r->l,rotp->p->v,rotp->p->c,&nwayss2);
			ps2=1.0/nwayss2;

			tree::tree_p tmerged=new tree;
			tmerged->p=rotp->p;

			merge3(rotp->p->r->l,rotp->p->r->r,rotp->p->r->v,rotp->p->r->c,&nwaysm1);
			pm1=1.0/nwaysm1;
			merge2(rotp->p->r->l,rotp->p->r->r,tmerged,rotp->p->r->v,rotp->p->r->c);
			rotp->p->r->p=0;
			delete rotp->p->r;
			rotp->p->r=tmerged;

			tmerged=new tree;
			tmerged->p=rotp->p;

			merge3(rotp->l,rotp->r,rotp->v,rotp->c,&nwaysm2);
			pm2=1.0/nwaysm2;
			size_t v,c;
			v=rotp->v;
			c=rotp->c;
			merge2(rotp->l,rotp->r,tmerged,rotp->v,rotp->c);
			rotp->p->l=tmerged;
			rotp->p=0;
			delete rotp;
			rotp=tnew->getptr(rotid);

		//end of merge code if rotp isleft.
			if( (rotp->v!=v && rotp->c!=c) && (rotp->p->r->v != v && rotp->p->r->c != c))
				hardreject=true;
			if(rotp->p->r->v==rotp->v && rotp->p->r->c==rotp->c && !isleaf(rotp->p->r))
				twowaystoinvert=true;

		}
	}
	else { //isright
		if(rotp->v==rotp->p->v) //special case, faster to handle it directly
		{
			rotleft(rotp);
			rotp=tnew->getptr(rotid);
			delete rotp->l;
			temp=rotp->r;
			rotp->p->r=temp;
			temp->p=rotp->p;
			rotp->r=0;
			rotp->l=0;
			rotp->p=0;
			delete rotp;
			rotp=tnew->getptr(rotid);
			//pm1=pm2=ps1=ps2=1.0 in this case
		}
		else
		{
			rotleft(rotp);
			rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
			reduceleft(rotp->p->l,rotp->p->v,rotp->p->c);
			rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
			reduceright(rotp,rotp->p->v,rotp->p->c);
			rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
			splitleft(rotp->p->l->l,rotp->p->v,rotp->p->c);
			splitright(rotp->l,rotp->p->v,rotp->p->c);

			merge3(rotp->p->l->l,rotp->l,rotp->p->v,rotp->p->c,&nwayss1);
			ps1=1.0/nwayss1;

			merge3(rotp->p->l->r,rotp->r,rotp->p->v,rotp->p->c,&nwayss2);
			ps2=1.0/nwayss2;

			tree::tree_p tmerged=new tree;
			tmerged->p=rotp->p;

			merge3(rotp->p->l->l,rotp->p->l->r,rotp->p->l->v,rotp->p->l->c,&nwaysm1);
			pm1=1.0/nwaysm1;
			merge2(rotp->p->l->l,rotp->p->l->r,tmerged,rotp->p->l->v,rotp->p->l->c);
			rotp->p->l->p=0;
			delete rotp->p->l;
			rotp->p->l=tmerged;

			tmerged=new tree;
			tmerged->p=rotp->p;

			merge3(rotp->l,rotp->r,rotp->v,rotp->c,&nwaysm2);
			pm2=1.0/nwaysm2;
			size_t v,c;
			v=rotp->v;
			c=rotp->c;
			merge2(rotp->l,rotp->r,tmerged,rotp->v,rotp->c);
			rotp->p->r=tmerged;
			rotp->p=0;
			delete rotp;
			rotp=tnew->getptr(rotid);

		//end of merge code if rotp isright
			if( (rotp->v!=v && rotp->c!=c) && (rotp->p->l->v != v && rotp->p->l->c != c))
				hardreject=true;
			if(rotp->p->l->v==rotp->v && rotp->p->l->c==rotp->c && !isleaf(rotp->p->l))
				twowaystoinvert=true;
		}
	}

	// Calculate prior probabilities, we just need to use the subtree where the rotation occured of tnew and x.
	subtold.clear();
	subtnew.clear();
	xp->p->getnodes(subtold);
	rotp->p->getnodes(subtnew);

	for(size_t i=0;i<subtold.size();i++) {
		if(subtold[i]->l) { //interior node
			priold*=pi.alpha/pow(1.0 + subtold[i]->depth(),pi.betap);
			goodvars.clear();
			getinternalvars(subtold[i],xi,goodvars);
			priold*=1.0/((double)goodvars.size()); //prob split on v 
			priold*=1.0/((double)getnumcuts(subtold[i],xi,subtold[i]->v)); //prob split on v at c is 1/numcutpoints
		}
		else //terminal node
			priold*=(1.0-pi.alpha/pow(1.0 + subtold[i]->depth(),pi.betap));	
	}
	for(size_t i=0;i<subtnew.size();i++) {
		if(subtnew[i]->l) { //interior node
			prinew*=pi.alpha/pow(1.0 + subtnew[i]->depth(),pi.betap);
			goodvars.clear();
			getinternalvars(subtnew[i],xi,goodvars);
			prinew*=1.0/((double)goodvars.size()); //prob split on v
			prinew*=1.0/((double)getnumcuts(subtnew[i],xi,subtnew[i]->v)); //prob split on v at c is 1/numcutpoints
			if(getnumcuts(subtnew[i],xi,subtnew[i]->v)<1)
			{
				x.pr(true);
				tnew->pr(true);
			}
		}
		else //terminal node
			prinew*=(1.0-pi.alpha/pow(1.0 + subtnew[i]->depth(),pi.betap));	
	}

	Qold_to_new=1.0/((double)rnodes.size()); //proposal probability of rotating from x to tnew
	
	rnodes.clear();
	tnew->getrotnodes(rnodes);  //this is very inefficient, could make it much nicer later on.

	if(!twowaystoinvert)
		Qnew_to_old=1.0/((double)rnodes.size()); //proposal probability of rotating from tnew back to x
	else
		Qnew_to_old=2.0/((double)rnodes.size());

	// Calculate log integrated likelihoods for the subtree where the rotation occured of tnew and x.
	double lilold=0.0,lilnew=0.0;
	std::vector<sinfo> sold,snew;
	nbold.clear();
	nbnew.clear();
	x.getbots(nbold);
	tnew->getbots(nbnew);

	allsuff(x,xi,di,nbold,sold);
	allsuff(*tnew,xi,di,nbnew,snew);
	for(size_t i=0;i<nbold.size();i++) {
				lilold += lil(sold[i].sum_y_div_v, sold[i].sum_v_inv, sold[i].n, pi, GAMMA_QUANT, TAU_SQ_QUANT); 
	}
	for(size_t i=0;i<nbnew.size();i++) {
      	if( snew[i].n >= minperbot ) {
		
			lilnew += lil(snew[i].sum_y_div_v, snew[i].sum_v_inv, snew[i].n, pi, GAMMA_QUANT, TAU_SQ_QUANT); 
		}
		else hardreject=true;
	}

	double alpha1;
	alpha1=prinew*Qnew_to_old/priold/Qold_to_new/pm1/pm2*ps1*ps2;

    double alpha2 = alpha1*exp(lilnew-lilold);
	double alpha = std::min(1.0,alpha2);

	if(hardreject) alpha=0.0;

	if(unif_rand()<alpha) {
		pi.taccept++;
		x = *tnew;
		return true;
	}
	else {
		return false;
	}

	return false;  // we never actually get here.
}




void calichgv(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double GAMMA_QUANT, double TAU_SQ_QUANT, size_t minperbot)
{
	tree::tree_p pertnode;
	unsigned int pertdx;
	std::vector<size_t> tryvars;
	size_t oldvar,newvar,oldcut,newcut;
	std::vector<std::vector<double> > cov;
	bool didswap=false;


	if(x.treesize()==1)
		return;

	// Get nog nodes and propose new split value
	tree::npv nognds;
	x.getnobots(nognds);  //This is all internal nodes of the tree.
	for(size_t k=0;k<nognds.size();k++) {
		if(unif_rand()<0.5) //50% of the time we do change of variable, otherwise perturb.
		{
			pi.tproposal++;
			pertdx=k;
			pertnode = nognds[pertdx];

			int L,U;
			double oldcutprob=1.0, newcutprob=1.0;

			getLU(pertnode,xi,&L,&U);
			oldcutprob=((double)(U-L+1.0));

			cov=pi.cvpm; // a copy of the prior.  Start with the prior at every node of the tree...
			didswap=false;
			oldvar=pertnode->getv();
			oldcut=pertnode->getc();
			updatecormat(pertnode,xi,cov);
			normchgvrow(oldvar,cov);

			//choose new variable randomly
			newvar=getchgv(oldvar,cov);
			pertnode->setv(newvar);
			if(cov[oldvar][newvar]<0.0) {
				pertnode->swaplr();
				didswap=true;
			}

			//get L,U for the new variable, save it and set newcut.
			int Ln,Un;
			size_t newcut;
			getLU(pertnode,xi,&Ln,&Un);
			newcutprob=((double)(Un-Ln+1.0));
			newcut=Ln+(size_t)(floor(unif_rand()*(Un-Ln+1.0)));
			pertnode->setc(newcut);

	        //now we also need to update the row of chgv for newv->oldv to calc MH correctly
	        updatecormat(pertnode,xi,cov);
	        normchgvrow(newvar,cov);

	        //sanity check:
	        if(cov[newvar][oldvar]==0.0)
	           Rcpp::Rcout << "Proposal newv cannot return to oldv!  This is not possible!" << std::endl;
	        double alpha0=cov[newvar][oldvar]/cov[oldvar][newvar];  //proposal ratio for newvar->oldvar and oldvar->newvar


	        // get suff stats and do accept/reject
			std::vector<sinfo> svold;
			std::vector<sinfo> svnew;
			tree::npv bnv;
			getpertsuff(x,pertnode,bnv,oldcut,oldvar,didswap,xi,di,svold,svnew);

			typedef tree::npv::size_type bvsz;
			double lilold,lilnew;
			bool hardreject=false;
			lilold=0.0;
			for(bvsz j=0;j!=svold.size();j++) {

					lilold += lil(svold[j].sum_y_div_v, svold[j].sum_v_inv, svold[j].n, pi, GAMMA_QUANT, TAU_SQ_QUANT); 
			}

			lilnew=0.0;
			for(bvsz j=0;j!=svnew.size();j++) {
				if( (svnew[j].n) < minperbot)
					hardreject=true;
					else
					lilnew += lil(svnew[j].sum_y_div_v, svnew[j].sum_v_inv, svnew[j].n, pi, GAMMA_QUANT, TAU_SQ_QUANT); 
			}

			double alpha1 = oldcutprob/newcutprob; //anything from the prior?
			double alpha2=alpha0*alpha1*exp(lilnew-lilold);
			double alpha = std::min(1.0,alpha2);

			if(hardreject) alpha=0.0;  //change of variable led to an bottom node with <5 observations in it, we reject this.

			if(unif_rand()<alpha) {
				pi.taccept++;
				pertnode->setv(newvar); //because the call to getpertsuff changes it back to oldvar to calc old lil.
				pertnode->setc(newcut); //because the call to getpertsuff changes it back to oldcut to calc the old lil.
				if(didswap)
					pertnode->swaplr();
			}
			//else nothing, pertnode->c and pertnode->v is already reset to the old value and if a swap was done it is already unswapped.
		}
		else  //do a perturb proposal
		{
			didswap=false;
			pi.tproposal++;
			pertdx=k;
			pertnode = nognds[pertdx];
			int L,U;
			double oldcutprob=1.0, newcutprob=1.0;

			// Get allowable range for perturbing cv at pertnode
			getLU(pertnode,xi,&L,&U);
			newcutprob=1.0/(U-L+1.0);

			size_t oldc = pertnode->getc();
			size_t propc;
			int ai,bi,oldai,oldbi;
			ai=(int)(floor(oldc-pi.pertalpha*(U-L+1.0)/2.0));
			bi=(int)(floor(oldc+pi.pertalpha*(U-L+1.0)/2.0));
			ai=std::max(ai,L);
			bi=std::min(bi,U);
			propc = ai + (size_t)(floor(unif_rand()*(bi-ai+1.0)));
			pertnode->setc(propc);
			oldai=(int)(floor(propc-pi.pertalpha*(U-L+1.0)/2.0));
			oldbi=(int)(floor(propc+pi.pertalpha*(U-L+1.0)/2.0));
			oldai=std::max(oldai,L);
			oldbi=std::min(oldbi,U);
			newcutprob=1.0/(bi-ai+1.0);
			oldcutprob=1.0/(oldbi-oldai+1.0);

			std::vector<sinfo> svold;
			std::vector<sinfo> svnew;
			tree::npv bnv;
			getpertsuff(x,pertnode,bnv,oldc,pertnode->getv(),didswap,xi,di,svold,svnew);

			typedef tree::npv::size_type bvsz;
			double lilold,lilnew;
			bool hardreject=false;
			lilold=0.0;
			for(bvsz j=0;j!=svold.size();j++) {

					lilold += lil(svold[j].sum_y_div_v, svold[j].sum_v_inv, svold[j].n, pi, GAMMA_QUANT, TAU_SQ_QUANT); 
			}

			lilnew=0.0;
			for(bvsz j=0;j!=svnew.size();j++) {
				if( (svnew[j].n) < minperbot)
					hardreject=true;
					else
					lilnew += lil(svnew[j].sum_y_div_v, svnew[j].sum_v_inv, svnew[j].n, pi, GAMMA_QUANT, TAU_SQ_QUANT); 
			}

			double alpha1 = oldcutprob/newcutprob; //anything from the prior?
			double alpha2= alpha1*exp(lilnew-lilold);
			double alpha = std::min(1.0,alpha2);

			if(hardreject) alpha=0.0;  //perturb led to an bottom node with <5 observations in it, we reject this.

			if(unif_rand()<alpha) {
				pi.taccept++;
				pertnode->setc(propc); //because the call to getpertsuff changes it back to oldc to calc the old lil.
			}
			//else nothing, pertnode->c and pertnode->v is already reset to the old value and if a swap was done it is already unswapped.	
		}
	}
}



//--------------------------------------------------
// Functions to support change-of-variable proposal
//--------------------------------------------------
// update the correlation matrix for chgv move taking into account that not all
// variables may be eligible at pertnode.
void updatecormat(tree::tree_p pertnode, xinfo& xi, std::vector<std::vector<double> >& chgv)
{
   int Ln,Un; //L,U for the ``new'' variable
   size_t oldv=pertnode->getv();
   size_t p=chgv.size();

   for(size_t i=0;i<p;i++) {
      if(i!=oldv && std::abs(chgv[oldv][i])>0.0) {
         if(chgv[oldv][i]<0.0)  //swap left,right branches
            pertnode->swaplr();
         getpertLU(pertnode,i,xi,&Ln,&Un);
         if(chgv[oldv][i]<0.0)  //undo the swap
            pertnode->swaplr();
         if(Un<Ln) //we can't transition to variable i here according to the tree structure
            chgv[oldv][i]=0.0;
      }
   }
}
//--------------------------------------------------
// renormalize the correlation matrix so that the probability of row sums to 1.
void normchgvrow(size_t row, std::vector<std::vector<double> >& chgv)
{
   double tmp=0.0;
   size_t p=chgv.size();

   for(size_t i=0;i<p;i++)
      tmp+=std::abs(chgv[row][i]);
   for(size_t i=0;i<p;i++)
      chgv[row][i]/=tmp;
}
//--------------------------------------------------
// randomly choose a new variable to transition to from oldv
size_t getchgv(size_t oldv, std::vector<std::vector<double> >& chgv)
{
   double cp=unif_rand();
   size_t p=chgv.size();
   size_t newv=oldv;
   std::vector<double> cumprob;

   cumprob=chgv[oldv];
   cumprob[1]=std::abs(cumprob[1]);
   for(size_t i=1;i<p;i++)
      cumprob[i]=std::abs(cumprob[i])+cumprob[i-1];

   for(size_t i=0;i<p;i++) {
      if(cumprob[i] >= cp) { //propose transition to this variable
         newv=i;
         i=p;  //break out of loop early
      }        
   }
   return newv;
}
