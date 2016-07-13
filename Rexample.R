ProjectLibrary <- "/..."
setwd(paste0(ProjectLibrary) )
library(devtools)
library(roxygen2)
library(RcppArmadillo)
library(rstudioapi)
library(Rcpp)
#RcppArmadillo.package.skeleton( "BayesQArt" )
#create("BayesQArt")
setwd(paste0(ProjectLibrary,"/BayesQArt/R"))
devtools::document
set.seed(9)
library("BayesQArt")

########################################### try friedman's data with only 5 predictors and mixture normal error#############
f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}

n = 100;
p = 10;
np = 100;
x=matrix(runif(n*p),n,p) 
xp=matrix(runif(np*p),np,p) 

fy = f(x)
fpy = f(xp)
pis = rbinom(n,1,.8);
y=fy+pis*rnorm(n) + (1-pis)*rnorm(n,mean = 3,sd = 3);

set.seed(99)
out = BayesQArt(y=y,X=x, Xtest = xp, quantile = 0.5, 
                burn = 3000, 
                nd = 5000,
                m = 200,
                min_obs_node= 5, 
                aa_parm = 2.0,
                bb_parm = 3.0,
                nc=50, 
                pbd=0.4, 
                pb=0.5, 
                alpha=0.95, 
                betap=2.0, 
                kappa = 2.0, 
                maxdept = 3)
            
wmad(fpy,out$pred_test,0.5)
