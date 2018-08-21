rm(list = ls())
library(fda)
library(MASS)

betaind = 3 # change
samplesize = 1:500 # change
#nsim = 20
Qs = 2:9
Q = length(Qs)

# Q = no. of pc components to employ when this is fixed.
setwd("~/Desktop/simulation/beta3/tr") # change
source('~/Desktop/simulation/datafunction/ExpectedMSE7_ty.r')
source('~/Desktop/simulation/datafunction/betafun.r')
load("~/Desktop/simulation/datafunction/dataSNR2.RData")


# Covariates + examples

X = x[-1,samplesize,] # dim(x) 101 500 200 
dim(X) # 100 X 100 X 200

Y = y[samplesize,,betaind]  
dim(Y) # 100 X 200



# Parametric approximation -- how many terms to use
fbasis = create.fourier.basis(c(0,1),nbasis=25)
psi = eval.basis(seq(0.01,1,by = 0.01),fbasis)
nlins = 1:3
nlin = length(nlins)

Lin = psi[,1:nlins[nlin],drop=FALSE]


#psi = princomp(X)$scores

nobs = 500
nsim = 200

#bpsi = eval.basis((0:364)+0.5,bbasis)
#psi = eval.basis((0:364)+0.5,fbasis)
#X =  psi%*%diag(exp( -(0:24)/4))%*%matrix(rnorm(25*nobs),25,nobs)
#dim(X) # 365 X 100
#X = as.matrix(read.table('X.dat'))
#lin = cbind(rep(1,nrow(X)),1:nrow(X),(1:nrow(X))^2)
#lin = lin[,1:2,drop=FALSE]
#nlins = 1:3
#nlin = length(nlins)
#Lin = psi[,1:nlins[nlin],drop=FALSE]
#  Yhat = 40*apply((2-seq(1:180)/90)*X[1:180,],2,mean) 
#if(truemod == 1){
#  Yhat = 3e3*apply(psi[1:183,1]*X[1:183,],2,sum)/365*sqrt(26.25/29.2)
#}else if(truemod == 2){
#  Yhat = 3e3*apply((psi[1:183,2]-min(psi[1:183,2]))*X[1:183,],2,sum)/365
#}else if(truemod == 3){
#  Yhat = 2e3*apply((psi[1:183,3]-min(psi[1:183,3]))*X[1:183,],2,sum)/365
#}
#  Yhat = 40*apply(X[1:180,],2,mean) 
#psi = bpsi[,-c(1,ncol(bpsi))]
#Q = ncol(psi)
### So lets start with something long

theta.range = c(0.2,1) # try theta = 0.10, 0.11, 0.12,..., 0.99, 1.00
theta.range.int = c(20, 100)
#theta.range = c(362,365)
thetas = seq(0.2, 1, by = 0.01) # length(thetas) = 81
thetas.int = seq(20, 100, by = 1)

lambdas = exp( seq(-35,5,len=51) )
ms = seq(1,15,by=2)  # ?? temporary


BICr = rep(0,nsim)

qthetas    = rep(0,nsim)            #?? parametric theta??
thetabest0 = rep(0,nsim) 
thetabest1 = rep(0,nsim) 

betaest0 = array(0,c(nsim,nrow(X)))
betaest1 = array(0,c(nsim,nrow(X)))
qbetaest = array(0,c(nsim,nrow(X))) # ?? # parametric beta??
betaesta = array(0,c(nsim,nrow(X)))      # ??: this is the initial estimates of beta without truncation

mchoice0      = rep(0,nsim)
lambdachoice0 = rep(0,nsim)
lambdachoice1 = rep(0,nsim)
aest0 = rep(0,nsim)
aest1 = rep(0,nsim)

mod1ferr = array(0,c(nsim,diff(theta.range.int)+1) )  # ?? model 1 error???


psis = array(0,c(nrow(X),max(ms),length(thetas))) # nrow(X) = 100, ncol(X) = 500
designs = array(0,c(ncol(X),max(ms),length(thetas))) # similar to matrix u


for(ss in 1:nsim)
{
  for(i in 1:length(thetas)){
    tpsi = eigen(X[1:thetas.int[i],,ss]%*%t(X[1:thetas.int[i],,ss]))$vectors[,1:max(ms),drop=FALSE] # i = 1, dim(tpsi) = 20 X 15
    psis[,,i] = rbind(tpsi,matrix(0,nrow = theta.range.int[2] - thetas.int[i], ncol = max(ms)))  # i = 1, dim(psis[,,i]) = 100 X 15
    
    designs[,,i] = t(X[,,ss])%*%psis[,,i]/dim(psis)[1]    
  }
  
  psi = tpsi  # dim(tpsi)  100 X 15
  
  # Start with an initial estimate
  
  objq = rep(0,Q)
  for(q in 1:Q){
    modi = fPCR(Y[,ss],X[,,ss],psi[,1:Qs[q],drop=FALSE],nrow(X))
    
    objq[q] = length(Y[,ss])*log( sum( (Y[,ss] - modi$fit)^2 ) ) + log(length(Y[,ss]))*(modi$rank)
  }
  qT = Qs[which.min(objq)]
  
  mod = fPCR(Y[,ss],X[,,ss],psi[,1:qT],nrow(X))
  r = mod$resid
  sig = mean(r^2)
  aest = mod$coef[1]
  betaest = psi[,1:qT]%*%mod$coef[2:(qT+1)]
  
  betaesta[ss,] = betaest
  
  # And a parametric approximation
  qmod = fPCRsearch0(Y[,ss],X[,,ss],psi = Lin,theta.range = theta.range.int)
  qa = qmod$mod$coef[1]
  qbeta = Lin%*%qmod$mod$coef[1 + (1:ncol(Lin))]
  qtheta = qmod$thetabest
  qbeta[ min(qtheta+1,length(qbeta)):length(qbeta) ] = 0
  qbetaest[ss,] = qbeta
  
  qthetas[ss] = qtheta
  qY = qmod$mod$fit
  
  
  BICr[ss] = length(Y[,ss])*log(mean( qmod$mod$resid^2 )) +  log(length(Y[,ss]))*(qmod$mod$rank)
  
  # Mimic parametric fit:
  
  lams0 = rep(0,length(ms))
  th0 = rep(0,length(ms))
  obj0 = rep(0,length(ms))
  
  for(k in 1:length(ms)){
    LTchoicek = ChooseLambdas0(qY,qbeta,sig,designs[,1:ms[k],,drop=FALSE],psis[,1:ms[k],,drop=FALSE],thetas,lambdas)
    lams0[k] = LTchoicek$lbest
    th0[k] = LTchoicek$thbest
    
    mod0f = fPCRsearch(Y[,ss],designs[,1:ms[k],,drop=FALSE],psis[,1:ms[k],,drop=FALSE],thetas,lambdas[lams0[k]])
    obj0[k] = length(Y)*log( sum( (Y - mod0f$mod$fit)^2 ) ) + log(length(Y))*(mod0f$mod$rank+1)
  }
  
  mchoice0[ss] = which.min(obj0)
  lambdachoice0[ss] = lams0[mchoice0[ss]]
  
  LTchoice1 =  ChooseLambdas1(qY,qbeta,sig,X[,,ss],psi[,1:qT],seq(theta.range.int[1],theta.range.int[2]),lambdas)
  
  lambdachoice1[ss] = LTchoice1$lbest
  #  thetabest1[ss] =   LTchoice1$thbest
  
  
  ############ START !!!!!##############3
  mod0f = fPCRsearch(Y[,ss],designs[,1:ms[mchoice0[ss]],,drop=FALSE],psis[,1:ms[mchoice0[ss]],,drop=FALSE],thetas,lambdas[lambdachoice0[ss]])
  
  mod1f = fPCRsearch2(Y[,ss],X[,,ss],aest,as.vector(betaest),theta.range.int,lambdas[lambdachoice1[ss]])
  
  mod1ferr[ss,] = mod1f[[4]]
  
  thetabest0[ss] = mod0f$thetabest
  thetabest1[ss] = mod1f$thetabest
  
  aest0[ss] = mod0f$a
  aest1[ss] = mod1f$a
  
  betaest0[ss,] = mod0f$beta
  betaest1[ss,] = mod1f$beta
  print(paste("n=",ss, sep = ""))
}

# change
save(thetabest0, thetabest1, betaest0, betaest1, mchoice0,lambdachoice0,
     lambdachoice1, aest0, aest1, mod1ferr, psis, designs, file = "beta3_500_SNR2.RData")



# method A
plot(seq(0.01,1,by=0.01), betaest0[1,],type = "l", ylim = c(-1,2),ylab = expression(beta(t)), xlab = "t")
for(j in 1:90)
{
  lines(seq(0.01,1,by=0.01), betaest0[j,],col=j+1)
}
abline(v = 0.5)



# method B
plot(seq(0.01,1,by=0.01), betaest1[1,],type = "l", ylim = c(-1,2),ylab = expression(beta(t)), xlab = "t")
for(j in 1:90)
{
  lines(seq(0.01,1,by=0.01), betaest1[j,],col=j+1)
}
abline(v = 0.5)


mean(thetabest0)
median(thetabest0)
sd(thetabest0)


mean(thetabest1/100)
median(thetabest1/100)
sd(thetabest1/100)



