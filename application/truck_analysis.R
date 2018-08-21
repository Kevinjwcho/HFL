rm(list = ls())
# data Y, X
library(fda)
library(MASS)
library(psych)

# Read in data
setwd("~/Desktop/application/truck")
data = as.matrix(read.table('Truck1Run1.csv',sep=',',header=FALSE))


# Scale Y up by 1000
Y = data[,1]*1e3
Y = Y - 6.174 # centered Y
# Covariates are everything but the first column
Z = t(data[,-1])

# Basis expansion
bbasis = create.bspline.basis(c(0,90),norder=4,breaks = seq(0,90,by=3))


# Smooth for covariates, centered
zfd = smooth.basis((0:89)+0.5,Z,fdPar(bbasis,20))
zct = center.fd(zfd$fd)

# Obtain acceleration

tt = seq(0,90,by=0.5)
X = -eval.fd(tt,zct,1)
nsec  = 121

X = X[1:nsec,] 

# Plot accelrations
set.seed(22)
set1 = sample(1:108, 10, replace = F)
tt = seq(0, 60, by = 0.5)
par(mar = c(4, 4, 2, 0.5))
matplot(tt,X[,set1],type='l',lty=1,col=1,xlab='Second',ylab='Acceleration', main = "(a)")
plot(tt, X[,set1[1]],type = "l", ylim = c(-5, 4), yaxt='n',
     xlab='Second',ylab='Acceleration', main = "(a)")
axis(2, at=seq(-6,4,by=2), labels=seq(-6,4,by=2),las=2)
matlines(tt,X[,set1], lty = 1, col =1)



# plot all acceleration curves
matplot(tt,X,type='l',lty=1,col=1,xlab='Second',yaxt='n',
        ylim=c(-8,4),ylab='Acceleration', main = "(a)")
axis(2, at=seq(-8,4,by=2), labels=seq(-8,4,by=2),las=2)
arrows(0,-8, 60, -8, code=1, length = 0.1)
text(4,-6.8, "Now")
text(56,-6.8, "Past")

# calculate U
h   = 60/120
cef = c(1, rep(c(4,2), (120 - 2)/2), 4, 1)

# bspline for regression
Mb = 120
norderb   = 4
nknotsb   = Mb + 1 
knotsb    = seq(0, 60, length.out = nknotsb)
nbasisb   = nknotsb + norderb - 2 
basisb    = create.bspline.basis(knotsb, nbasisb, norderb)

basismatb = eval.basis(knotsb, basisb)
u   = h/3*t(X)%*%diag(cef)%*%basismatb # 108 X 53

# calculate V
basismatv = eval.basis(knotsb, basisb, 2) 
v = h/3*t(basismatv)%*%diag(cef)%*%basismatv # 103 X 53


#####################################
#        Smoothing Spline           #
#####################################
# source files that define functions
path = "~/Desktop/simulation/datafunction"
source(paste(path, "fun_b_smooth.R", sep = "/"))
source(paste(path, "fun_ystar_ustar.R", sep = "/"))
source(paste(path, "fun_beta_results.R", sep = "/"))
source(paste(path, "fun_AIC_BIC.R", sep = "/"))

# tunning parameter alpha for smooth method
Msim = Mb
dsim = norderb -1
nsample=108
alpha_s  = exp(seq(-5,13,len = 20))
#alpha_s  = 0
alpha_s_AIC = numeric(1)
alpha_s_BIC = numeric(1)

AIC_s      = rep(NA, length(alpha_s))
BIC_s      = rep(NA, length(alpha_s))
b_s_AIC    = rep(NA, Msim + dsim)
b_s_BIC    = rep(NA, Msim + dsim)
b_s_temp   = array(NA, c(Msim + dsim, length(alpha_s)))

for(iteralpha in 1:length(alpha_s))
{
  ptm = proc.time()
  Y_s  = Y
  U_s  = u

  b_s_temptemp = smoothb(Y = Y_s, U = U_s, V = v, n = length(Y_s), alpha = alpha_s[iteralpha], M = Msim, d = dsim)
  b_s_temp[, iteralpha] = b_s_temptemp
  
  AICBIC_s_temp = AICBIC.smooth.fun(Y = Y_s, b = b_s_temptemp, U = U_s, V = v, n = length(Y_s), alpha  = alpha_s[iteralpha])
  AIC_s[iteralpha]  = AICBIC_s_temp$AIC 
  BIC_s[iteralpha]  = AICBIC_s_temp$BIC 
  timeint = as.numeric((proc.time() - ptm)[3])
  print(paste("alpha =", alpha_s[iteralpha], "time", round(timeint,2)))
}

# alpha picked by AIC and BIC
plot(alpha_s, AIC_s)
plot(alpha_s, BIC_s)
alpha_s_AIC = alpha_s[which.min(AIC_s)]
alpha_s_BIC = alpha_s[which.min(BIC_s)]
b_s_AIC =  b_s_temp[, which.min(AIC_s)]
b_s_BIC =  b_s_temp[, which.min(BIC_s)]
beta_s_AIC = b_s_AIC%*%t(basismatb)
beta_s_BIC = b_s_BIC%*%t(basismatb)

# plot beta
plot(tt[1:121], beta_s_AIC, lty = 2, lwd = 2,
     xlab = "Second", ylab = expression(beta(t)))
abline(h = 0)

plot(tt[1:121], beta_s_BIC, type = "l", lwd = 2,ylim = c(-0.2, 0.8),yaxt = 'n',
     xlab = "Second", ylab = expression(beta(t)), main = "(b)")
abline(h = 0)
axis(2, at=seq(-0.2,0.8,by=0.2), labels=seq(-0.2,0.8,by=0.2),las=2)
arrows(0,-0.2, 60, -0.2, code=1, length = 0.1)
text(4,-0.1, "Now")
text(56,-0.1, "Past")


matplot(tt,X,type='l',lty=1,col=1,xlab='Second',yaxt='n',
        ylim=c(-8,4),ylab='Acceleration', main = "(a)")
axis(2, at=seq(-8,4,by=2), labels=seq(-8,4,by=2),las=2)
abline(h= 0)


#####################################
#          NGR Estimator           #
#####################################
# source files that define functions
library(glmnet)
path = "~/Desktop/simulation/datafunction"
source(paste(path, "fun_theta.R", sep = "/"))
source(paste(path, "fun_hg.R", sep = "/"))
source(paste(path, "fun_cj.R", sep = "/"))
source(paste(path, "fun_bhat_gbridge.R", sep = "/"))

Msim  = Mb
dsim  = 3
nsample  = length(Y)
nitersim = 100


# tunning parameter tau, gamma, alpha for nested group bridge method
gamma_gbr  = 0.5
alpha_gbr  = exp(seq(6,8.3,len = 8))
tau_gbr    = exp(seq(-15,-10,len = 8))
tau_gbr_AIC2   = numeric(1)
tau_gbr_BIC2   = numeric(1)
alpha_gbr_AIC2 = numeric(1)
alpha_gbr_BIC2 = numeric(1)
AIC2_gbr    = array(NA, c(length(alpha_gbr), length(tau_gbr)))
BIC2_gbr    = array(NA, c(length(alpha_gbr), length(tau_gbr)))

b_gbr_AIC2   = rep(NA, Msim + dsim)
b_gbr_BIC2   = rep(NA, Msim + dsim)
b_gbr_temp   = array(NA, c(Msim + dsim, length(alpha_gbr), length(tau_gbr)))

b0gbr  = b_s_BIC # use bsmooth based on BIC

for(iteralpha in 1:length(alpha_gbr))
{
  ptm = proc.time()
  for(itertau in 1:length(tau_gbr))
  {
    Ygbr  = Y
    Ugbr  = u
    b_gbr_temptemp = bhat_gbridge(b0 = b0gbr, Y = Ygbr, u = Ugbr, v = v, n = length(Ygbr), 
                                  alpha = alpha_gbr[iteralpha], tau = tau_gbr[itertau], gamma = gamma_gbr, 
                                  niter = nitersim, M = Msim, d = dsim) 
    b_gbr_temp[, iteralpha, itertau] = b_gbr_temptemp$bhatgbr
    Y_gbr_hat = Ugbr%*%b_gbr_temptemp$bhatgbr
    g_gbr     = b_gbr_temptemp$gk
    W_al      = diag(length(Ygbr)*g_gbr/abs(b_gbr_temptemp$bhatgbr))
    ABIC_gbr_temp = AICBIC.gbr.fun(Y = Ygbr, b = b_gbr_temptemp$bhatgbr, U = Ugbr, V = v, W = W_al, 
                                   n = length(Ygbr), alpha  = alpha_gbr[iteralpha])
 
    AIC2_gbr[iteralpha, itertau]  = ABIC_gbr_temp$AIC2 
    BIC2_gbr[iteralpha, itertau]  = ABIC_gbr_temp$BIC2
  }
  timeint = as.numeric((proc.time() - ptm)[3])
  print(paste("alpha =", alpha_gbr[iteralpha], "completed", round(timeint)))
}


# alpha and tau picked by AIC and BIC
idx_gbr_aic2 = which(AIC2_gbr == min(AIC2_gbr), arr.ind = TRUE)
idx_gbr_bic2 = which(BIC2_gbr == min(BIC2_gbr), arr.ind = TRUE)
idx_gbr_bic2

tau_gbr_AIC2   = tau_gbr[idx_gbr_aic2[2]]
tau_gbr_BIC2   = tau_gbr[idx_gbr_bic2[2]]
alpha_gbr_AIC2 = alpha_gbr[idx_gbr_aic2[1]]
alpha_gbr_BIC2 = alpha_gbr[idx_gbr_bic2[1]]

b_gbr_AIC2   = b_gbr_temp[, idx_gbr_aic2[1], idx_gbr_aic2[2]]
b_gbr_BIC2   = b_gbr_temp[, idx_gbr_bic2[1], idx_gbr_bic2[2]]
beta_gbr_AIC = b_gbr_AIC2%*%t(basismatb)
beta_gbr_BIC = b_gbr_BIC2%*%t(basismatb)

# plot beta
plot(tt[1:121], beta_gbr_AIC, type = "l", 
     xlab = "seconds", ylab = expression(beta(t)))

plot(tt[1:121], beta_gbr_BIC, type = "l", lwd = 2, main = "(c)", 
     yaxt = 'n', ylim = c(-0.2, 0.8), xlab = "Second", ylab = expression(beta(t)))
lines(tt[1:121], beta_s_BIC, lwd = 2, lty = 2, col = "grey")
abline(h = 0)
z = seq(-0.2, 0.8, by = 0.2)
axis(2, at=seq(-0.2, 0.8, by = 0.2), labels=seq(-0.2, 0.8, by = 0.2),las=2)

#save(b_gbr_BIC2 , beta_gbr_BIC2, file = "truck_ngr.RData")


#####################################
#          BOOTSTRAPPING            #
#####################################
Yhat = u%*%b_gbr_BIC2
err  = Y - Yhat

Bsim = 200

b_gbr_boot = array(NA, c(Msim + dsim, Bsim))
tau_gbr_BIC2_boot = rep(NA, Bsim)
alpha_gbr_BIC2_boot = rep(NA, Bsim)
set.seed(100)

nitersim = 500
alpha_s_boot = exp(seq(5,13,len = 20))
alpha_gbr_boot  = exp(seq(6,8,len = 8))
tau_gbr_boot    = exp(seq(-15,-10,len = 8))


BIC2_gbr_boot    = array(NA, c(length(alpha_gbr_boot), length(tau_gbr_boot)))
b_gbr_BIC2_boot  = rep(NA, Msim + dsim)
b_gbr_t_boot     = array(NA, c(Msim + dsim, length(alpha_gbr_boot), length(tau_gbr_boot)))

for(iboot in 1:Bsim)
{ ptm = proc.time()
  YstarY = Yhat + sample(err, replace = T)
  
  # obtain smoothing spline estimator
  b_s_t   = array(NA, c(Msim + dsim, length(alpha_s_boot)))
  AIC_s_boot = rep(NA, Bsim)
  BIC_s_boot = rep(NA, Bsim)
  
  for(iteralpha in 1:length(alpha_s_boot))
  {
    Y_s  = YstarY
    U_s  = u
    
    b_s_tt = smoothb(Y = Y_s, U = U_s, V = v, n = length(Y_s), alpha = alpha_s_boot[iteralpha], M = Msim, d = dsim)
    b_s_t[, iteralpha] = b_s_tt
    
    AICBIC_s_t = AICBIC.smooth.fun(Y = Y_s, b = b_s_t[,iteralpha], U = U_s, V = v, n = length(Y), alpha  = alpha_s_boot[iteralpha])
    
    AIC_s_boot[iteralpha]  = AICBIC_s_t$AIC 
    BIC_s_boot[iteralpha]  = AICBIC_s_t$BIC 
  }
  
  b_s_AIC_boot =  b_s_t[, which.min(AIC_s)]
  b_s_BIC_boot =  b_s_t[, which.min(BIC_s)]
  
  # ngr estimator
  nitersim = 500
  
  BIC2_gbr_boot    = array(NA, c(length(alpha_gbr_boot), length(tau_gbr_boot)))
  
  b_gbr_BIC2_boot  = rep(NA, Msim + dsim)
  b_gbr_t_boot     = array(NA, c(Msim + dsim, length(alpha_gbr_boot), length(tau_gbr_boot)))
  
  for(iteralpha in 1:length(alpha_gbr_boot))
  { 
  b0gbr = b_s_BIC_boot
  for(itertau in 1:length(tau_gbr_boot))
  {
    Ygbr  = YstarY
    Ugbr  = u
    b_gbr_tt = bhat_gbridge(b0 = b0gbr, Y = Ygbr, u = Ugbr, v = v, n = length(Ygbr), 
                            alpha = alpha_gbr_boot[iteralpha], tau = tau_gbr_boot[itertau], gamma = gamma_gbr, 
                            niter = nitersim, M = Msim, d = dsim) 
    b_gbr_t_boot[, iteralpha, itertau] = b_gbr_tt$bhatgbr
    
    if(sum(b_gbr_tt$bhatgbr == 0) == length(b_gbr_tt$bhatgbr)) {BIC2_gbr_boot[iteralpha, itertau] = 10^(20)} else{
      Y_gbr_hat = Ugbr%*%b_gbr_tt$bhatgbr
      g_gbr     = b_gbr_tt$gk
      W_al      = diag(length(Ygbr)*g_gbr/abs(b_gbr_tt$bhatgbr))
      ABIC_gbr_t = AICBIC.gbr.fun(Y = Ygbr, b = b_gbr_tt$bhatgbr, U = Ugbr, V = v, W = W_al, 
                                  n = length(Ygbr), alpha  = alpha_gbr_boot[iteralpha])
      
      BIC2_gbr_boot[iteralpha, itertau]  = ABIC_gbr_t$BIC2}
  }
  }
  
  # alpha and tau picked by AIC and BIC
  idx_gbr_bic2_boot = which(BIC2_gbr_boot == min(BIC2_gbr_boot), arr.ind = TRUE)
  
  tau_gbr_BIC2_boot[iboot]   = tau_gbr_boot[idx_gbr_bic2_boot[2]]
  alpha_gbr_BIC2_boot[iboot] = alpha_gbr_boot[idx_gbr_bic2_boot[1]]
  
  #plot(tau_gbr, BIC2_gbr, type = "b", xlab = expression(tau), ylab = "BIC", xaxt = "n") # so choose tau_gbr[5]
  
  b_gbr_boot[, iboot] = b_gbr_t_boot[, idx_gbr_bic2_boot[1], idx_gbr_bic2_boot[2]]
  timeint = as.numeric((proc.time() - ptm)[3])
  print(paste("n =", iboot, "completed", timeint))
}

save(b_gbr_boot,beta_gbr_BIC,basismatb, file = "truckboot.RData")
#load("truck_ngr.RData")
#load("truckboot.RData")
load("~/Desktop/Application/truck/truckboot.RData")

beta_gbr = basismatb%*%b_gbr_BIC2

# Pointwise standard deviation of betas
beta_gbr_boot = basismatb%*%b_gbr_boot  
beta_gbr_sd = apply(beta_gbr_boot,1,sd)     # ngr

# 
# for(i in 1:200)
#     {
#       plot(tt[1:121], beta_gbr_boot[,i], type = "l", 
#            main = i)
#     }

# Create and plot confidence intervals
# ngr ci piviotal interval
beta_gbr_ci = matrix(NA, nrow = length(beta_gbr), ncol = 3)
qa = 0.05
for(ici in 1:length(beta_gbr))
{
  qtt = quantile(beta_gbr_boot[ici,], probs = c(1-qa/2, qa/2))
  beta_gbr_ci[ici,] = c(2*beta_gbr[ici] - qtt[1], beta_gbr[ici],2*beta_gbr[ici] - qtt[2]) 
}

matplot(tt[1:121],beta_gbr_ci,type='l', 
        main = "Pivotal Interval", xlab = "seconds", ylab = expression(beta(t)))
matplot(tt[1:121],beta_gbr_ci,type='l', 
        xlab = "seconds", ylab = expression(beta(t)))
par(mar=c(4,4,2,0.5))
matplot(tt[1:121],beta_gbr_ci,type='l',lwd=c(1,2,1),lty=c(2,1,2),xlab='Second',
        yaxt='n', ylab=expression(beta(t)), col = "black", main = "(b)")
z = seq(0, 1, by = 0.2)
#axis(2, at = z, labels = z)
axis(2, at=z, labels=z,las=2)
abline(h = 0)


par(mar=c(4,4,2,0.5))
plot(ttu[1:(mu2 + 1)],beta_gbr_ci[,2],type='l',lwd = 2,xlab='seconds',ylim = c(-110,170),
     ylab=expression(beta(t)), col = "black", main = "(b)")
lines(ttu[1:(mu2 + 1)],beta_gbr_ci[,1],type='l',lwd = 1, lty = 2)
lines(ttu[1:(mu2 + 1)],beta_gbr_ci[,3],type='l',lwd = 1, lty = 2)
lines(ttu[1:(mu2 + 1)],beta_s_BIC,type='l',lwd = 1, lty = 1, col = 2)
lines(ttu[1:(mu2 + 1)],beta_gbr_ci[,2],type='l',lwd = 2)
abline(h = 0)
