#ifndef GUARD_tree_h
#define GUARD_tree_h

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstddef>
#include <vector>
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include "info.h"


/*
Class is a bit confusing because sometimes you
want to think about an instance as a node and
sometimes as the whole tree including the node
and all its children.
Fundamentally, we have a tree not a node.
*/

/*
three ways to access a node:
(i) node id: is the integer assigned by the node numbering system
     assuming they are all there
(ii) node ind: is the index into the array of
   node(tree) pointers returned by getnodes or getbots or getnogs
   which means you go left to right across the bottom of the tree
(iii) by its pointer (should only "give out" const pointers)
*/

//info contained in a node, used by input operator
struct node_info {
   std::size_t id; //node id
   std::size_t v;  //variable
   std::size_t c;  //cut point
   double m;       //mu
};

//xinfo xi, then xi[v][c] is the c^{th} cutpoint for variable v.
// left if x[v] < xi[v][c]
typedef std::vector<double> vec_d; //double vector
typedef std::vector<vec_d> xinfo; //vector of vectors, will be split rules

class tree {
public:
   
  //------------------------------
   //typedefs
   typedef tree* tree_p;
   typedef const tree* tree_cp;
   typedef std::vector<tree_p> npv; //Node Pointer Vector
   typedef std::vector<tree_cp> cnpv; //const Node Pointer Vector

   
   //------------------------------
   //friends
    friend bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi,double GAMMA_QUANT,
            double TAU_SQ_QUANT, size_t min_obs_node,std::vector<int>& varused, size_t maxdepth);
	friend void getLU(tree::tree_p pertnode, xinfo& xi, int* L, int* U);
    friend void getpertLU(tree::tree_p pertnode, size_t pertvar, xinfo& xi, int* L, int* U);
    friend void getpertsuff(tree& x, tree::tree_p pertnode, tree::npv& bnv, size_t oldc, size_t oldvar, bool didswaplr, xinfo& xi, dinfo& di, std::vector<sinfo>& svold, std::vector<sinfo>& svnew);
    friend void getinternalvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars);
    friend int getnumcuts(tree::tree_p n, xinfo& xi, size_t var);
	friend void rotright(tree::tree_p n);
	friend void rotleft(tree::tree_p n);
   friend void reduceleft(tree::tree_p n, size_t v, size_t c);
   friend void reduceright(tree::tree_p n, size_t v, size_t c);
   friend bool isequalvc(tree::tree_p t1, tree::tree_p t2);
   friend bool isequal(tree::tree_p t1, tree::tree_p t2);
    friend void splitleft(tree::tree_p t, size_t v, size_t c);
   friend void splitright(tree::tree_p t, size_t v, size_t c);
   friend bool splitsonv(tree::tree_p t, size_t v);
   friend bool splitsonv(tree::tree_p nl, tree::tree_p nr, size_t v);
   friend void nwaysonsplit(tree::tree_p tl, tree::tree_p tr, int *nways);
   friend bool hasvcsplit(tree::tree_p t, size_t v, size_t c);
   friend bool isleaf(tree::tree_p t);
   friend bool arenodesequal(tree::tree_p nl, tree::tree_p nr);
   friend bool arenodesleafs(tree::tree_p nl, tree::tree_p nr);
    friend bool canmerge(tree::tree_p tl, tree::tree_p tr, size_t v, size_t c, int* nways);
   friend bool merge2(tree::tree_p tl, tree::tree_p tr, tree::tree_p t, size_t v, size_t c);
   friend bool merge3(tree::tree_p tl, tree::tree_p tr, size_t v, size_t c, int* nways);
   friend bool merge2nt(tree::tree_p tl, tree::tree_p tr, tree::tree_p t, size_t v, size_t c, int* nways, int atway);
   friend bool merge3nt(tree::tree_p tl, tree::tree_p tr, size_t v, size_t c, int* nways);
   friend bool merge2triv(tree::tree_p tl, tree::tree_p tr, tree::tree_p t, size_t v, size_t c);
   
    
	friend bool calirotp(tree::tree_p tnew, tree& x, xinfo& xi, dinfo& di, pinfo& pi, double GAMMA_QUANT, double TAU_SQ_QUANT, size_t minperbot);
	friend void calichgv(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double GAMMA_QUANT, double TAU_SQ_QUANT, size_t minperbot);

   //------------------------------
   //tree constructors, destructors
   tree();
   tree(const tree&);
   tree(double);
   ~tree() {tonull();}

   //------------------------------
   //operators
   tree& operator=(const tree&);

   //------------------------------
   //node access
   //you are freely allowed to change mu, v, c
   //set----------
   void setm(double mu) {this->mu=mu;}
   void setv(size_t v) {this->v = v;}
   void setc(size_t c) {this->c = c;}
   //get----------
   double getm() const {return mu;}
   size_t getv() const {return v;}
   size_t getc() const {return c;}
   tree_p getp() const {return p;}  //should this be tree_cp?
   tree_p getl() const {return l;}
   tree_p getr() const {return r;}

   //------------------------------
   //tree functions
   //stuff about the tree----------
   size_t treesize() const; //number of nodes in tree

   size_t nnogs() const;    //number of nog nodes (no grandchildren nodes)
   size_t nbots() const;    //number of bottom nodes
   void pr(bool pc=true) const; //to screen, pc is "print children"
   //birth death using nid----------
   bool birth(size_t nid,size_t v, size_t c, double ml, double mr);
   bool death(size_t nid, double mu);
   //vectors of node pointers----------
   void getbots(npv& bv);         //get bottom nodes
   void getinternals(npv& internv); //get internal nodes
   void getnogs(npv& nv);         //get nog nodes (no granchildren)
   void getnodes(npv& v);         //get vector of all nodes
   void getnodes(cnpv& v) const;  //get all nodes
   
   
   void getintnodesnotonv(npv& nv, size_t var); //get alll nodes except either leafs or nodes splitting on var
   void getnodesonv(npv& v, size_t var); //get all nodes splitting on var
   void getnodesonvc(npv& v, size_t var, size_t cut); //get all nodes splitting on var at cut
   void getnobots(npv& v);            //get all nodes except leafs
   void getrotnodes(npv& v);          //get rot nodes
   void getrotelems(npv& v);
   
   //find node from x and region for var----------
   tree_cp bn(size_t i, dinfo& di,xinfo& xi); //find bottom node for x
   void rg(size_t v, int* L, int* U); //find region [L,U] for var v.
    void rl(size_t v, int *L);
   void ru(size_t v, int *U);
   bool isrightchildofvsplit(size_t v);  //is this node in a right subtree that splits on variable v?
   bool isleftchildofvsplit(size_t v);   //is this node in a left subtree that splits on variable v?
   // swap the left and right branches of this trees node.
	// does NOT check that this is not a bottom node!
	void swaplr();
   size_t nuse(size_t v); //how many times var v is used in a rule.
   //------------------------------
   //node functions
   size_t depth() const; //depth of a node
   size_t nid() const;   //node id
   char ntype() const;   //t:top;b:bottom;n:nog;i:interior, carefull a t can be bot
   tree_p getptr(size_t nid); //get node pointer from node id, 0 if not there.
   bool isnog() const;
   void tonull(); //like a "clear", null tree has just one node
   bool isleft() const;
   bool isright() const;

private:
   //------------------------------
   //parameter for node
   double mu;
   //------------------------------
   //rule: left if x[v] < xinfo[v][c]
   size_t v;
   size_t c;
   //------------------------------
   //tree structure
   tree_p p; //parent
   tree_p l; //left child
   tree_p r; //right child
   //------------------------------
   //utiity functions
   void cp(tree_p n,  tree_cp o); //copy tree o to n
   void birthp(tree_p np,size_t v, size_t c, double ml, double mr);
   void deathp(tree_p nb, double mu); //kill children of nog node nb
};


#endif
