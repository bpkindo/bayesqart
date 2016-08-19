#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <cmath>
#include <iostream>
#include "info.h"
#include "tree.h"
//function to draw truncated normal
//above=1 means from above b=trunpt, a=-inf
//above=0 means from below a=trunpt, b= +inf   

double rtrun(double mu, double sigma,double trunpt, int above);
void makexinfo(arma::uword p, arma::uword n, arma::mat const& x, xinfo& xi, arma::uword nc);
//fit double array fv
void fit(tree& t, xinfo& xi, dinfo& di, arma::vec& fv);

//--------------------------------------------------
//--------------------------------------------------
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots, size_t maxdepth);
//--------------------------------------------------
//does this bottom node n have any variables it can split on.
bool cansplit(tree::tree_p n, xinfo& xi);

//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi, std::vector<size_t>& goodvars);

//get sufficient stats for children (v,c) of node nx in tree x
void getsuff(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);
//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
void getsuff(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);

//--------------------------------------------------
// similar, but in case of perturb.
void getpertsuff(tree& x, tree::tree_p pertnode, tree::npv& bnv, size_t oldc, xinfo& xi, dinfo& di, std::vector<sinfo>& svold, std::vector<sinfo>& svnew);
void getpertsuff(tree& x, tree::tree_p pertnode, tree::npv& bnv, size_t oldc, size_t oldvar, bool didswaplr, xinfo& xi, dinfo& di, std::vector<sinfo>& svold, std::vector<sinfo>& svnew);
//--------------------------------------------------
//log of the integrated likelihood
double lil(double sum_y_div_v, double sum_v_inv, size_t n, pinfo& pi,double GAMMA_QUANT, double TAU_SQ_QUANT);
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else a/(1+d)^b
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi);
//--------------------------------------------------
// draw all the bottom node mu's
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double GAMMA_QUANT,
            double TAU_SQ_QUANT);
//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv);
//--------------------------------------------------
// Get the L,U values for a node in the tree *given* the tree
// structure both above and below that node.
void getLU(tree::tree_p pertnode, xinfo& xi, int* L, int* U);
//--------------------------------------------------
// similar except we get it for a prescribed variable pertvar
void getpertLU(tree::tree_p pertnode, size_t pertvar, xinfo& xi, int* L, int* U);
//--------------------------------------------------
// number of available cutpoints at node n for variable var
int getnumcuts(tree::tree_p n, xinfo& xi, size_t var);
//--------------------------------------------------
// Find numbr of variables internal tree node n can split on
void getinternalvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars);

//--------------------------------------------------
// number of available cutpoints at node n for variable var
int getnumcuts(tree::tree_p n, xinfo& xi, size_t var);
//--------------------------------------------------
// Rotate a given rotatable node in the tree.
// node n must be the left child of n->p and also not a root or terminal leaf node.
void rotright(tree::tree_p n);
void rotleft(tree::tree_p n);
void reduceleft(tree::tree_p n, size_t v, size_t c);
void reduceright(tree::tree_p n, size_t v, size_t c);
void splitleft(tree::tree_p t, size_t v, size_t c);
void splitright(tree::tree_p t, size_t v, size_t c);
bool merge2(tree::tree_p tl, tree::tree_p tr, tree::tree_p t, size_t v, size_t c);
// only to get nways, not to actually do the merge.
bool merge3(tree::tree_p tl, tree::tree_p tr, size_t v, size_t c, int* nways);
bool hasvcsplit(tree::tree_p t, size_t v, size_t c);
bool splitsonv(tree::tree_p t, size_t v);
bool splitsonv(tree::tree_p nl, tree::tree_p nr, size_t v);
bool isleaf(tree::tree_p t);
bool arenodesequal(tree::tree_p nl, tree::tree_p nr);
bool arenodesleafs(tree::tree_p nl, tree::tree_p nr);
bool calirotp(tree::tree_p tnew, tree& x, xinfo& xi, dinfo& di, pinfo& pi, double GAMMA_QUANT, double TAU_SQ_QUANT, size_t minperbot);

void getpertsuff(tree& x, tree::tree_p pertnode, tree::npv& bnv, size_t oldc, size_t oldv, bool didswaplr, xinfo& xi, dinfo& di, std::vector<sinfo>& svold, std::vector<sinfo>& svnew);
void getpertsuff(tree& x, tree::tree_p pertnode, tree::npv& bnv, size_t oldc, xinfo& xi, dinfo& di,
                     std::vector<sinfo>& svold, std::vector<sinfo>& svnew);
void calichgv(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double GAMMA_QUANT, double TAU_SQ_QUANT, size_t minperbot); 
// variables may be eligible at pertnode.
void updatecormat(tree::tree_p pertnode, xinfo& xi, std::vector<std::vector<double> >& chgv);
//--------------------------------------------------
// renormalize the correlation matrix so that the probability of row sums to 1.
void normchgvrow(size_t row, std::vector<std::vector<double> >& chgv);
//--------------------------------------------------
// randomly choose a new variable to transition to from oldv
size_t getchgv(size_t oldv, std::vector<std::vector<double> >& chgv);					 