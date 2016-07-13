#include <R.h>
#include <cmath>
#include <iostream>
#include "RcppArmadillo.h"
#include <Rcpp.h>
#include "rgig.h"

double _gig_mode(double lambda, double omega)
/*---------------------------------------------------------------------------*/
/* Compute mode of GIG distribution.                                         */
/*                                                                           */
/* Parameters:                                                               */
/*   lambda .. parameter for distribution                                    */
/*   omega ... parameter for distribution                                    */
/*                                                                           */
/* Return:                                                                   */
/*   mode                                                                    */
/*---------------------------------------------------------------------------*/
{
  if (lambda >= 1.0)
    /* mode of fgig(x) */
    return (sqrt((lambda-1.0)*(lambda-1.0) + omega*omega)+(lambda-1.0))/omega;
  else
    /* 0 <= lambda < 1: use mode of f(1/x) */
    return omega / (sqrt((1.0-lambda)*(1.0-lambda) + omega*omega)+(1.0-lambda));
} /* end of _gig_mode() */


double _rgig_ROU_noshift ( double lambda, double lambda_old, double omega, double alpha)
/*---------------------------------------------------------------------------*/
/* Tpye 1:                                                                   */
/* Ratio-of-uniforms without shift.                                          */
/*   Dagpunar (1988), Sect.~4.6.2                                            */
/*   Lehner (1989)                                                           */
/*---------------------------------------------------------------------------*/
{
    double rv;
  double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
  double ym, um;     /* location of maximum of x*sqrt(f(x)); umax of MBR */
  double s, t;       /* auxiliary variables */
  double U, V, X;    /* random variables */

  int count = 0;     /* counter for total number of iterations */

  /* -- Setup -------------------------------------------------------------- */

  /* shortcuts */
  t = 0.5 * (lambda-1.0);
  s = 0.25 * omega;

  /* mode = location of maximum of sqrt(f(x)) */
xm = _gig_mode(lambda, omega);

  /* normalization constant: c = log(sqrt(f(xm))) */
  nc = t*log(xm) - s*(xm + 1.0/xm);

  /* location of maximum of x*sqrt(f(x)):           */
  /* we need the positive root of                   */
  /*    omega/2*y^2 - (lambda+1)*y - omega/2 = 0    */
  ym = ((lambda+1.0) + sqrt((lambda+1.0)*(lambda+1.0) + omega*omega))/omega;

  /* boundaries of minmal bounding rectangle:                   */
  /* we us the "normalized" density f(x) / f(xm). hence         */
  /* upper boundary: vmax = 1.                                  */
  /* left hand boundary: umin = 0.                              */
  /* right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(xm)) */
  um = exp(0.5*(lambda+1.0)*log(ym) - s*(ym + 1.0/ym) - nc);

  /* -- Generate sample ---------------------------------------------------- */


    do {
      ++count;
      U = um * R::runif(0,1);        /* U(0,umax) */
      V = R::runif(0,1);             /* U(0,vmax) */
      X = U/V;
    }                              /* Acceptance/Rejection */
    while (((log(V)) > (t*log(X) - s*(X + 1.0/X) - nc)));

    /* store random point */
    rv = (lambda_old < 0.0) ? (alpha / X) : (alpha * X);

return rv;
  /* -- End ---------------------------------------------------------------- */

} /* end of _rgig_ROU_noshift() */





//_rgig_ROU_shift_alt
double _rgig_ROU_shift_alt ( double lambda, double lambda_old, double omega, double alpha)
/*---------------------------------------------------------------------------*/
/* Type 8:                                                                   */
/* Ratio-of-uniforms with shift by 'mode', alternative implementation.       */
/*   Dagpunar (1989)                                                         */
/*   Lehner (1989)                                                           */
/*---------------------------------------------------------------------------*/
{
double rv;
      //  Rcpp::Rcout << "_rgig_ROU_shift_alt" << endl;
  double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
  double s, t;       /* auxiliary variables */
  double U, V, X;    /* random variables */



  int count = 0;     /* counter for total number of iterations */

  double a, b, c;    /* coefficent of cubic */
  double p, q;       /* coefficents of depresressed cubic */
  double fi, fak;    /* auxiliary results for Cardano's rule */

  double y1, y2;     /* roots of (1/x)*sqrt(f((1/x)+m)) */

  double uplus, uminus;  /* maximum and minimum of x*sqrt(f(x+m)) */

  /* -- Setup -------------------------------------------------------------- */

  /* shortcuts */
  t = 0.5 * (lambda-1.0);
  s = 0.25 * omega;

  /* mode = location of maximum of sqrt(f(x)) */
  xm = _gig_mode(lambda, omega);

  /* normalization constant: c = log(sqrt(f(xm))) */
  nc = t*log(xm) - s*(xm + 1.0/xm);

  /* location of minimum and maximum of (1/x)*sqrt(f(1/x+m)):  */

  /* compute coeffients of cubic equation y^3+a*y^2+b*y+c=0 */
  a = -(2.0*(lambda+1.0)/omega + xm);       /* < 0 */
  b = (2.0*(lambda-1.0)*xm/omega - 1.0);
  c = xm;

  /* we need the roots in (0,xm) and (xm,inf) */

  /* substitute y=z-a/3 for depressed cubic equation z^3+p*z+q=0 */
  p = b - a*a/3.0;
  q = (2.0*a*a*a)/27.0 - (a*b)/3.0 + c;

  /* use Cardano's rule */
  fi = acos(-q/(2.0*sqrt(-(p*p*p)/27.0)));
  fak = 2.0*sqrt(-p/3.0);
  y1 = fak * cos(fi/3.0) - a/3.0;
  y2 = fak * cos(fi/3.0 + 4.0/3.0*M_PI) - a/3.0;

  /* boundaries of minmal bounding rectangle:                  */
  /* we us the "normalized" density f(x) / f(xm). hence        */
  /* upper boundary: vmax = 1.                                 */
  /* left hand boundary: uminus = (y2-xm) * sqrt(f(y2)) / sqrt(f(xm)) */
  /* right hand boundary: uplus = (y1-xm) * sqrt(f(y1)) / sqrt(f(xm)) */
  uplus  = (y1-xm) * exp(t*log(y1) - s*(y1 + 1.0/y1) - nc);
  uminus = (y2-xm) * exp(t*log(y2) - s*(y2 + 1.0/y2) - nc);

  /* -- Generate sample ---------------------------------------------------- */


    do {
      ++count;
      U = uminus + R::runif(0,1)  * (uplus - uminus);    /* U(u-,u+)  */
      V = R::runif(0,1);                                /* U(0,vmax) */
      X = U/V + xm;
    }                                         /* Acceptance/Rejection */
    while ((X <= 0.0) || ((log(V)) > (t*log(X) - s*(X + 1.0/X) - nc)));

    /* store random point */
    rv= (lambda_old < 0.0) ? (alpha / X) : (alpha * X);

return rv;

}



double  _rgig_newapproach1 ( double lambda, double lambda_old, double omega, double alpha)
/*---------------------------------------------------------------------------*/
/* Type 4:                                                                   */
/* New approach, constant hat in log-concave part.                           */
/* Draw sample from GIG distribution.                                        */
/*                                                                           */
/* Case: 0 < lambda < 1, 0 < omega < 1                                       */
/*                                                                           */
/* Parameters:                                                               */
/*   n ....... sample size (positive integer)                                */
/*   lambda .. parameter for distribution                                    */
/*   omega ... parameter for distribution                                    */
/*                                                                           */
/* Return:                                                                   */
/*   random sample of size 'n'                                               */
/*---------------------------------------------------------------------------*/
{
    double rv;
//    Rcpp::Rcout << "new method" << endl;
  /* parameters for hat function */
  double A[3], Atot;  /* area below hat */
  double k0;          /* maximum of PDF */
  double k1, k2;      /* multiplicative constant */

  double xm;          /* location of mode */
  double x0;          /* splitting point T-concave / T-convex */
  double a;           /* auxiliary variable */

  double U, V, X;     /* random numbers */
  double hx;          /* hat at X */


  int count = 0;      /* counter for total number of iterations */

  /* -- Check arguments ---------------------------------------------------- */

  if (lambda >= 1.0  || omega >1.0)
     {error( "invalid parameters\n");}

  /* -- Setup -------------------------------------------------------------- */

  /* mode = location of maximum of sqrt(f(x)) */
  xm = _gig_mode(lambda, omega);
  /* splitting point */
  x0 = omega/(1.0-lambda);

  /* domain [0, x_0] */
  k0 = exp((lambda-1.0)*log(xm) - 0.5*omega*(xm + 1.0/xm));     /* = f(xm) */
  A[0] = k0 * x0;

  /* domain [x_0, Infinity] */
  if (x0 >= 2.0/omega) {
    k1 = 0.0;
    A[1] = 0.0;
    k2 = pow(x0, lambda-1.0);
    A[2] = k2 * 2.0 * exp(-omega*x0/2.0)/omega;
  }

  else {
    /* domain [x_0, 2/omega] */
    k1 = exp(-omega);
    A[1] = (lambda == 0.0)
      ? k1 * log(2.0/(omega*omega))
      : k1 / lambda * ( pow(2.0/omega, lambda) - pow(x0, lambda) );

    /* domain [2/omega, Infinity] */
    k2 = pow(2/omega, lambda-1.0);
    A[2] = k2 * 2 * exp(-1.0)/omega;
  }

  /* total area */
  Atot = A[0] + A[1] + A[2];

  /* -- Generate sample ---------------------------------------------------- */


    do {
      ++count;

      /* get uniform random number */
      V = Atot * R::runif(0,1);

      do {

	/* domain [0, x_0] */
	if (V <= A[0]) {
	  X = x0 * V / A[0];
	  hx = k0;
	  break;
	}

	/* domain [x_0, 2/omega] */
	V -= A[0];
	if (V <= A[1]) {
	  if (lambda == 0.0) {
	    X = omega * exp(exp(omega)*V);
	    hx = k1 / X;
	  }
	  else {
	    X = pow(pow(x0, lambda) + (lambda / k1 * V), 1.0/lambda);
	    hx = k1 * pow(X, lambda-1.0);
	  }
	  break;
	}

	/* domain [max(x0,2/omega), Infinity] */
	V -= A[1];
	a = (x0 > 2.0/omega) ? x0 : 2.0/omega;
	X = -2.0/omega * log(exp(-omega/2.0 * a) - omega/(2.0*k2) * V);
	hx = k2 * exp(-omega/2.0 * X);
	break;

      } while(0);

      /* accept or reject */
      U = R::runif(0,1) * hx;

      if (log(U) <= (lambda-1.0) * log(X) - omega/2.0 * (X+1.0/X)) {
	/* store random point */
	rv = (lambda_old < 0.0) ? (alpha / X) : (alpha * X);
	break;
      }
    } while(1);


  /* -- End ---------------------------------------------------------------- */

return rv;
} /* end of _rgig_newapproach1() */


// [[Rcpp::depends(RcppArmadillo)]]

// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
double do_rgig(double lambda, double chi, double psi){
    double rv ;
    double omega, alpha;     /* parameters of standard distribution */
    /* check GIG parameters: */
  if ( !(std::isfinite(lambda) && std::isfinite(chi) && std::isfinite(psi)) ||
       (chi <  0.0 || psi < 0)      ||
       (chi == 0.0 && lambda <= 0.0) ||
       (psi == 0.0 && lambda >= 0.0) ) {
    {error(" invalid GIG arguments\n"); }
    }

  if (chi < ZTOL) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.00) {
      rv =  R::rgamma(lambda, 2.0/psi); 
    }
    else {
      rv = 1.0/R::rgamma(-lambda, 2.0/psi); 
    }
  }

    else if (psi < ZTOL) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      rv = 1.0/R::rgamma(lambda, 2.0/chi); 
    }
    else {
      rv = R::rgamma(-lambda, 2.0/chi); 
    }

  }

    else {
    double lambda_old = lambda;
    if (lambda < 0.0) lambda = -lambda;
    alpha = sqrt(chi/psi);
    omega = sqrt(psi*chi);

//Rcpp::Rcout <<"labmda, omega: "<< lambda << ", " << omega << endl;

    /* run generator */
    do {
      if (lambda > 2.0 || omega > 3.0) {
        /* Ratio-of-uniforms with shift by 'mode', alternative implementation */
      rv =   _rgig_ROU_shift_alt(lambda, lambda_old, omega, alpha);
        break;
      }

      if (lambda >= 1.00  - 2.25*omega*omega || omega > 0.2) {
        /* Ratio-of-uniforms without shift */
      //  Rcpp::Rcout <<  "_rgig_ROU_noshift" << endl;
        rv = _rgig_ROU_noshift(lambda, lambda_old, omega, alpha);

        break;
      }

      if (lambda >= 0.0 && omega > 0.0) {
        /* New approach, constant hat in log-concave part.0 */
       rv =  _rgig_newapproach1(lambda, lambda_old, omega, alpha);
        break;
      }

      /* else */
      error("parameters must satisfy lambda>=0 and omega>0.0 \n") ;

    } while (0);
  }

return rv;
}
