#include <vector>
#include <map>
#include "tree.h"


//--------------------------------------------------
// constructors
tree::tree(): mu(0.0),v(0),c(0),p(0),l(0),r(0) {}
tree::tree(double m): mu(m),v(0),c(0),p(0),l(0),r(0) {}
tree::tree(const tree& n): mu(0.0),v(0),c(0),p(0),l(0),r(0) { cp(this,&n); }
//--------------------------------------------------
//operators
tree& tree::operator=(const tree& rhs)
{
   if(&rhs != this) {
      tonull(); //kill left hand side (this)
      cp(this,&rhs); //copy right hand side to left hand side
   }
   return *this;
}
//--------------------------------------------------
//public functions
// find bottom node pointer given x
//--------------------
tree::tree_cp tree::bn(size_t i, dinfo& di,xinfo& xi)
{
   if(l==0) return this; //bottom node
   if(di.x(i,v) < xi[v][c]) {
      return l->bn(i,di,xi);
   } else {
      return r->bn(i,di,xi);
   }
}
//--------------------
//find region for a given variable
void tree::rg(size_t v, int* L, int* U)
{ 
   if(p==0)  { //no parent
      return;
   }
   if(p->v == v) { //does my parent use v?
      if(this == p->l) { //am I left or right child
         if((int)(p->c) <= (*U)) *U = (p->c)-1; 
      } else {
         if((int)(p->c) >= *L) *L = (p->c)+1; 
      }
   } 
   p->rg(v,L,U);
}


void tree::rl(size_t v, int *L)
{
	if(l==0) { //no children
		return;
	}

	if(this->v==v && (int)(this->c) >= (*L)) {
		*L=(int)(this->c)+1;
		r->rl(v,L);
	}
	else {
		l->rl(v,L);
		r->rl(v,L);	
	}
	
}

void tree::ru(size_t v, int *U)
{
	if(l==0) { //no children
		return;
	}

	if(this->v==v && (int)(this->c) <= (*U)) {
		*U=(int)(this->c)-1;
		l->ru(v,U);
	}
	else {
		l->ru(v,U);
		r->ru(v,U);
	}

}

bool tree::isrightchildofvsplit(size_t v)
{
	if(p==0) //root node
		return false;
	else if(p->v!=v)
		return p->isrightchildofvsplit(v);
	else if(p->r==this)
		return true;
	return false;
}

bool tree::isleftchildofvsplit(size_t v)
{
	if(p==0) //root node
		return false;
	else if(p->v!=v)
		return p->isleftchildofvsplit(v);
	else if(p->l==this)
		return true;
	return false;
}

// swap the left and right branches of this trees node.
// does NOT check that this is not a bottom node!
void tree::swaplr()
{
	tree_p tempt;
	tempt=this->r;
	this->r=this->l;
	this->l=tempt;
}

//--------------------
//tree size
size_t tree::treesize() const
{
   if(l==0) return 1;  //if bottom node, tree size is 1
   else return (1+l->treesize()+r->treesize());
}


 
//--------------------
size_t tree::nnogs() const
{
   if(!l) return 0; //bottom node
   if(l->l || r->l) { //not a nog
      return (l->nnogs() + r->nnogs()); 
   } else { //is a nog
      return 1;
   }
}
size_t tree::nuse(size_t v)
{
   npv nds;
   this->getnodes(nds);
   size_t nu=0; //return value
   for(size_t i=0;i!=nds.size();i++) {
      if(nds[i]->l && nds[i]->v==v) nu+=1;
   }
   return nu;
}
//--------------------
size_t tree::nbots() const
{
   if(l==0) { //if a bottom node
      return 1;
   } else {
      return l->nbots() + r->nbots();
   }
}
//--------------------
//depth of node
size_t tree::depth() const
{
   if(!p) return 0; //no parents
   else return (1+p->depth());
}
//--------------------
// node id
size_t tree::nid() const
//recursion up the tree
{
   if(!p) return 1; //if you don't have a parent, you are the top
   if(this==p->l) return 2*(p->nid()); //if you are a left child
   else return 2*(p->nid())+1; //else you are a right child
}
//--------------------
//node type
char tree::ntype() const
{
   //t:top, b:bottom, n:no grandchildren, i:internal
   if(!p) return 't';
   if(!l) return 'b';
   if(!(l->l) && !(r->l)) return 'n';
   return 'i';
}
//--------------------
//get bottom nodes
//recursion down the tree
void tree::getbots(npv& bv) 
{
   if(l) { //have children
      l->getbots(bv);
      r->getbots(bv);
   } else {
      bv.push_back(this);
   }
}

//get internal nodes
//recursion down the tree
void tree::getinternals(npv& internv) 
{
   if(l) { //have children
      internv.push_back(this);
      l->getinternals(internv);
      r->getinternals(internv);
   } 
      
}
//--------------------
//get nog nodes
//recursion down the tree
void tree::getnogs(npv& nv) 
{
   if(l) { //have children
      if((l->l) || (r->l)) {  //have grandchildren
         if(l->l) l->getnogs(nv);
         if(r->l) r->getnogs(nv);
      } else {
         nv.push_back(this);
      }
   }
}



void tree::getintnodesnotonv(npv& nv, size_t var)
{
	if(l) { //have children
		if(this->v != var)      //as long as this doesn't split on var, it's eligible
			nv.push_back(this);
		l->getintnodesnotonv(nv,var);
		r->getintnodesnotonv(nv,var);
	}
}

void tree::getnodesonv(npv& v, size_t var)
{
	if(this->v==var)
		v.push_back(this);
	if(l) {
		l->getnodesonv(v,var);
		r->getnodesonv(v,var);
	}
}

void tree::getnodesonvc(npv& v, size_t var, size_t cut)
{
	if(this->v==var && this->c==cut)
		v.push_back(this);
	if(l) {
		l->getnodesonvc(v,var,cut);
		r->getnodesonvc(v,var,cut);
	}
}

//--------------------
//get all nodes
//recursion down the tree
void tree::getnodes(npv& v) 
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
void tree::getnodes(cnpv& v)  const
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
// Get nodes that are not leafs
void tree::getnobots(npv& v)
{
	if(this->l) //left node has children
	{
		v.push_back(this);
		this->l->getnobots(v);
		if(this->r->l)
			this->r->getnobots(v);
	}
}
//--------------------
// Get nodes of tree minus root and leafs
void tree::getrotnodes(npv& v)
{
	if(!this->p && this->l)  //this is the root node and it has children, so lets get the rot nodes
	{
		this->l->getnobots(v);
		this->r->getnobots(v);
	}
}

void tree::getrotelems(npv& v)
{
	if(this->l) //left node has children
	{
		if(this->v != this->p->v) v.push_back(this);
		this->l->getrotelems(v);
		if(this->r->l)
			this->r->getrotelems(v);
	}

}

//--------------------
//add children to  bot node nid
bool tree::birth(size_t nid,size_t v, size_t c, double ml, double mr)
{
   tree_p np = getptr(nid);
   if(np==0) {
      Rcpp::Rcout << "error in birth: bottom node not found\n";
      return false; //did not find note with that nid
   }
   if(np->l) {
      Rcpp::Rcout << "error in birth: found node has children\n";
      return false; //node is not a bottom node
   }

   //add children to bottom node np 
   tree_p l = new tree;
   l->mu=ml;
   tree_p r = new tree;
   r->mu=mr;
   np->l=l;
   np->r=r;
   np->v = v; np->c=c;
   l->p = np;
   r->p = np;

   return true;
}
//--------------------
//is the node a nog node
bool tree::isnog() const
{
   bool isnog=true;
   if(l) {
      if(l->l || r->l) isnog=false; //one of the children has children.
   } else {
      isnog=false; //no children
   }
   return isnog;
}
bool tree::isleft() const
{
	bool isleft=false;
	if(p && p->l==this)
		isleft=true;

	return isleft;
}

bool tree::isright() const
{
	bool isright=false;
	if(p && p->r==this)
		isright=true;

	return isright;
}

//--------------------
//kill children of  nog node nid
bool tree::death(size_t nid, double mu) 
{
   tree_p nb = getptr(nid);
   if(nb==0) {
      Rcpp::Rcout << "error in death, nid invalid\n";
      return false;
   }
   if(nb->isnog()) {
      delete nb->l;
      delete nb->r;
      nb->l=0;
      nb->r=0;
      nb->v=0;
      nb->c=0;
      nb->mu=mu;
      return true;
   } else {
      Rcpp::Rcout << "error in death, node is not a nog node\n";
      return false;
   }
}
//--------------------
//add children to bot node *np
void tree::birthp(tree_p np,size_t v, size_t c, double ml, double mr)
{
   tree_p l = new tree;
   l->mu=ml;
   tree_p r = new tree;
   r->mu=mr;
   np->l=l;
   np->r=r;
   np->v = v; np->c=c;
   l->p = np;
   r->p = np;
}
//--------------------
//kill children of  nog node *nb
void tree::deathp(tree_p nb, double mu) 
{
   delete nb->l;
   delete nb->r;
   nb->l=0;
   nb->r=0;
   nb->v=0;
   nb->c=0;
   nb->mu=mu;
}

//--------------------
// Print out tree(pc=true) or node(pc=false) information
// Uses recursion down
void tree::pr(bool pc) const
{
	size_t d = depth();
	size_t id = nid();

	size_t pid; 
	if(!p) pid=0; //parent of top node
	else pid = p->nid();

	std::string pad(2*d,' ');
	std::string sp(", ");
	if(pc && (ntype()=='t')) 
		Rcpp::Rcout << "tree size: " << treesize() << std::endl;
//	cout << pad << "(id,parent): " << id << sp << pid;
	Rcpp::Rcout << pad << "id: " << id;
	Rcpp::Rcout << sp << "(v,c): " << v << sp << c;
	Rcpp::Rcout << sp << "mu: " << mu; 
	Rcpp::Rcout << sp << "type: " << ntype();
	Rcpp::Rcout << sp << "depth: " << depth();
//	Rcpp::Rcout << sp << "pointer: " << this << std::endl;
//	Rcpp::Rcout << sp << "left: " << this->l << std::endl;
//	Rcpp::Rcout << sp << "right:" << this->r << std::endl;
	Rcpp::Rcout << std::endl;

	if(pc) {
		if(l) {
			l->pr(pc);
			r->pr(pc);
		}
	}
}
//--------------------------------------------------
//private functions
//--------------------
//copy tree o to tree n
void tree::cp(tree_p n, tree_cp o)
//assume n has no children (so we don't have to kill them)
//recursion down
{
   if(n->l) {
      Rcpp::Rcout << "cp:error node has children\n";
      return;
   }

   n->mu = o->mu;
   n->v = o->v;
   n->c = o->c;

   if(o->l) { //if o has children
      n->l = new tree;
      (n->l)->p = n;
      cp(n->l,o->l);
      n->r = new tree;
      (n->r)->p = n;
      cp(n->r,o->r);
   }
}
//--------------------
//cut back to one node
void tree::tonull()
{
   size_t ts = treesize();
   while(ts>1) { //if false ts=1
      npv nv;
      getnogs(nv);
      for(size_t i=0;i<nv.size();i++) {
         delete nv[i]->l;
         delete nv[i]->r;
         nv[i]->l=0;
         nv[i]->r=0;
      }
      ts = treesize(); 
   }
   mu=0.0;
   v=0;c=0;
   p=0;l=0;r=0;
}
//--------------------
// get pointer for node from its nid
tree::tree_p tree::getptr(size_t nid)
{
   if(this->nid() == nid) return this; //found it
   if(l==0) return 0; //no children, did not find it
   tree_p lp = l->getptr(nid);
   if(lp) return lp; //found on left
   tree_p rp = r->getptr(nid); 
   if(rp) return rp; //found on right
   return 0; //never found it
}
