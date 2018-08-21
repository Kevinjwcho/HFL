rm(list=ls())
setwd("~/Desktop/simulation/beta3/NGR") # change
library(psych)
library(fda)

betaind = 3 # change
samplesize = 1:500 # change

# source functions and load data
path = "~/Desktop/simulation/datafunction"
source(paste(path, "betafun.R", sep = "/"))  
source(paste(path, "fun_b_smooth.R", sep = "/"))
source(paste(path, "fun_ystar_ustar.R", sep = "/"))
source(paste(path, "fun_beta_results.R", sep = "/"))
source(paste(path, "fun_AIC_BIC.R", sep = "/"))
load("~/Desktop/simulation/datafunction/dataSNR2.RData")

# Ybeta1, ..., Ybeta5: matrix, 500 X 200
y100 = y[samplesize,,betaind] 
dim(y100)  # 500 200
u100 = u[samplesize,,]   
dim(u100)  # 500 103 200

# we generate nsim number of dataset
Msim  = 100
dsim  = 3
nsim = 200

# tunning parameter alpha for smooth method
alpha_s     = 10^(-(10:3))
alpha_s_AIC = rep(NA, nsim)
alpha_s_BIC = rep(NA, nsim)

AIC_s      = array(NA, c(length(alpha_s), nsim))
BIC_s      = array(NA, c(length(alpha_s), nsim))
b_s_AIC    = array(NA, c(Msim + dsim, nsim))
b_s_BIC    = array(NA, c(Msim + dsim, nsim))
b_s_temp   = array(NA, c(Msim + dsim, length(alpha_s), nsim))

for(itersim in 1:nsim)
{ptm = proc.time()
for(iteralpha in 1:length(alpha_s))
{
  Y_s  = y100[, itersim]
  U_s  = u100[,,itersim]
  b_s_temptemp = smoothb(Y = Y_s, U = U_s, V = v, n = length(Y_s), alpha = alpha_s[iteralpha], M = Msim, d = dsim)
  b_s_temp[, iteralpha, itersim] = b_s_temptemp
  Y_s_hat = U_s%*%b_s_temptemp
  AICBIC_s_temp = AICBIC.smooth.fun(Y = Y_s, b = b_s_temptemp, U = U_s, V = v, n = length(Y_s), alpha  = alpha_s[iteralpha])
  AIC_s[iteralpha, itersim]  = AICBIC_s_temp$AIC 
  BIC_s[iteralpha, itersim]  = AICBIC_s_temp$BIC 
}
timeint = as.numeric((proc.time() - ptm)[3])
print(paste("n =", itersim, "completed", timeint))
}

# alpha picked by AIC and BIC
idx_b_s_AIC = rep(NA, nsim)
idx_b_s_BIC = rep(NA, nsim)
for(alpi in 1:nsim)
{
  alpha_s_AIC[alpi] = alpha_s[which.min(AIC_s[,alpi])]
  alpha_s_BIC[alpi] = alpha_s[which.min(BIC_s[,alpi])]
  b_s_AIC[,alpi]    =  b_s_temp[, which.min(AIC_s[,alpi]), alpi]
  b_s_BIC[,alpi]    =  b_s_temp[, which.min(BIC_s[,alpi]), alpi]
}
hist(alpha_s_AIC)
hist(alpha_s_BIC)
table(alpha_s_AIC)
table(alpha_s_BIC)

save(alpha_s_BIC, b_s_BIC, file="bt3_ss_n500_SNR2.RData")   # change




# smooth method #
M        = 100
nknots   = M + 1
tobs     = seq(0,1,length.out=M+1)
norder   = 4
knots    = seq(0,1,length.out=nknots)
nbasis   = nknots + norder - 2 
basis    = create.bspline.basis(knots,nbasis,norder)
basismat = eval.basis(tobs, basis)
basismat.tobs  = basismat
# True #
beta_true  = apply(as.matrix(tobs), 1, beta_fun, ii = betaind)

# smooth method #
# bsmooth chosen by AIC
par(mfrow = c(1,2))
beta_hat_smooth_AIC  = apply(b_s_AIC, 2, FUN = betahat)
plot(tobs, beta_true, type="l", ylim = c(-1, 2.5), ylab = expression(beta[3](t)), xlab = "t", main = "Smooth AIC")   ### SUBJECT TO CHANGE
matlines(tobs, beta_hat_smooth_AIC)
lines(tobs, beta_true, lwd = 3)
mtext("100 samples  100 simulations")

# bsmooth chosen by BIC
beta_hat_smooth_BIC  = apply(b_s_BIC, 2, FUN = betahat)
plot(tobs, beta_true, type="l", ylim = c(-1, 2.5), ylab = expression(beta[3](t)), xlab = "t", main = "Smooth BIC")   ### SUBJECT TO CHANGE
matlines(tobs, beta_hat_smooth_BIC)
lines(tobs, beta_true, lwd = 3)
lines(tobs, beta_true, lwd = 3)
mtext("100 samples  100 simulations")