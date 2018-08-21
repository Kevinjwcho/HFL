rm(list=ls())
setwd("~/Desktop/simulation/datafunction")
source("betafun.R")
library(fda)
domain = c(0,1)
# generate n single sample path X(t) at discrete grids t1,t2,...,tm
s  = 200 #number of simulations
m  = 101 #number of observations
n  = 500 #sample size
nb = 3  #number of models for beta
M  = 100 #M+1:number of knots of B-spline basis
d  = 3 #degree of B-spline basis
tobs = seq(domain[1],domain[2],by=1/(m-1)) #observation grid

# define bsplines to generate X
norderx   = 4
nknotsx   = 64 
knotsx    = seq(0,1, length.out = nknotsx)
nbasisx   = nknotsx + norderx - 2 
basisx    = create.bspline.basis(knotsx,nbasisx,norderx)
basismatx  = eval.basis(tobs, basisx) # 101 66

# generate x
set.seed(20)
x = array(NA, c(m,n,s))
for(ss in 1:s)
{
  for(j in 1:n)
  {
    x[,j,ss] = rnorm(nbasisx, 0, 1)%*%t(basismatx)
  }
}

# beta
betav = array(NA, c(m, 3))
for(i in 1:3)
{
  betav[,i] = apply(as.matrix(tobs), 1, beta_fun, ii = i)
}

# y0 the signals
h = 1/(m-1)
cef = c(1, rep(c(4,2), (m-3)/2), 4, 1)
y0 = array(NA, c(n,s,nb))
for(ss in 1:s)
{
  y0[,ss,] =  h/3*t(x[,,ss])%*%diag(cef)%*%betav 
}

# eps0
eps0 = array(NA, c(s,nb))
for(i in 1:nb)
{
 eps0[,i] = apply(y0[,,i], 2, sd)
}

# y
y =  array(NA, c(n,s,nb))
for(i in 1:nb)
{
  for(ss in 1:s)
  {
    y[,ss,i] = y0[,ss,i] + rnorm(n,mean = 0, sd = eps0[ss,i]/sqrt(2))
  }
}

# define bsplines to approximate beta
norder   = d+1
nknots   = M+1
knots    = seq(0,1, length.out = nknots)
nbasis   = nknots + norder - 2 
basis    = create.bspline.basis(knots,nbasis,norder)
basismat  = eval.basis(tobs, basis) # 101 103
basismat2  = eval.basis(tobs, basis,2) # 101 103

# matrix u
u = array(NA, c(n,M+d,s))
for(ss in 1:s)
{
  u[,,ss] = h/3*t(x[,,ss])%*%diag(cef)%*%basismat
}

# matrix v
v = h/3*t(basismat2)%*%diag(cef)%*%basismat2

save(x, y, u, v, file = "dataSNR2.RData")
