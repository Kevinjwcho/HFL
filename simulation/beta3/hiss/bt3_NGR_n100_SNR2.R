rm(list=ls())
setwd("~/Desktop/simulation/beta3/NGR") # change
library(psych)
library(fda)
library(glmnet)

betaind = 3 # change
samplesize = 1:100 # change
load("~/Desktop/simulation/beta3/NGR/bt3_ss_n100_SNR2.RData") # change

# source functions and load data
path = "~/Desktop/simulation/datafunction"
source(paste(path, "betafun.R", sep = "/"))      
source(paste(path, "fun_AIC_BIC.R", sep = "/"))
source(paste(path, "fun_theta.R", sep = "/"))
source(paste(path, "fun_hg.R", sep = "/"))
source(paste(path, "fun_cj.R", sep = "/"))
source(paste(path, "fun_bhat_gbridge.R", sep = "/"))
source(paste(path, "fun_ystar_ustar.R", sep = "/"))
source(paste(path, "fun_beta_results.R", sep = "/"))
load("~/Desktop/simulation/datafunction/dataSNR2.RData")


# Ybeta1, ..., Ybeta5: matrix, 500 X 200
y100 = y[samplesize,,betaind] 
dim(y100)  # 100 200
u100 = u[samplesize,,]   
dim(u100)  # 100 103 200

# we generate nsim number of dataset
Msim  = 100
dsim  = 3
nsim = 200
nitersim = 500

# optimal bsmooth and group bridge for nsim = 400 samples
# tunning parameter tau, gamma, alpha for group bridge method
gamma_gbr  = 0.5
alpha_gbr  = 10^(-(9:7))
tau_gbr = exp(seq(-38,-30,len=20))
tau_gbr_AIC2   = rep(NA, nsim)
tau_gbr_BIC2   = rep(NA, nsim)
alpha_gbr_AIC2   = rep(NA, nsim)
alpha_gbr_BIC2   = rep(NA, nsim)

AIC2_gbr    = array(NA, c(length(alpha_gbr), length(tau_gbr), nsim))
BIC2_gbr    = array(NA, c(length(alpha_gbr), length(tau_gbr), nsim))

b_gbr_AIC2   = array(NA, c(Msim + dsim, nsim))
b_gbr_BIC2   = array(NA, c(Msim + dsim, nsim))
b_gbr_temp   = array(NA, c(Msim + dsim, length(alpha_gbr), length(tau_gbr), nsim))

for(itersim in 1:nsim)
{
  ptm = proc.time()
  b0gbr    = b_s_BIC[, itersim] # use bsmooth based on BIC
  for(iteralpha in 1:length(alpha_gbr))
  {
    for(itertau in 1:length(tau_gbr))
    {
        Ygbr  = y100[, itersim]
        Ugbr  = u100[,,itersim]
        b_gbr_temptemp = bhat_gbridge(b0 = b0gbr, Y = Ygbr, u = Ugbr, v = v, n = length(Ygbr), 
                                     alpha = alpha_gbr[iteralpha], tau = tau_gbr[itertau], gamma = gamma_gbr, 
                                     niter = nitersim, M = Msim, d = dsim) 
        b_gbr_temp[, iteralpha, itertau, itersim] = b_gbr_temptemp$bhatgbr
        Y_gbr_hat = Ugbr%*%b_gbr_temptemp$bhatgbr
        g_gbr     = b_gbr_temptemp$gk
        W_al      = diag(length(Ygbr)*g_gbr/abs(b_gbr_temptemp$bhatgbr))
        ABIC_gbr_temp = AICBIC.gbr.fun(Y = Ygbr, b = b_gbr_temptemp$bhatgbr, U = Ugbr, V = v, W = W_al, 
                                       n = length(Ygbr), alpha  = alpha_gbr[iteralpha])
        #AIC1_gbr[iteralpha, itertau, itersim]  = ABIC_gbr_temp$AIC1 
        AIC2_gbr[iteralpha, itertau, itersim]  = ABIC_gbr_temp$AIC2 
        #BIC1_gbr[iteralpha, itertau, itersim]  = ABIC_gbr_temp$BIC1
        BIC2_gbr[iteralpha, itertau, itersim]  = ABIC_gbr_temp$BIC2
    }
  }
  timeint = as.numeric((proc.time() - ptm)[3])
  print(paste("n =", itersim, "completed", timeint))
}



# # alpha and tau picked by AIC and BIC
for(ali in 1:nsim)
{
  idx_gbr_aic2 = which(AIC2_gbr[,,ali] == min(AIC2_gbr[,,ali]), arr.ind = TRUE)
  idx_gbr_bic2 = which(BIC2_gbr[,,ali] == min(BIC2_gbr[,,ali]), arr.ind = TRUE)
  
  tau_gbr_AIC2[ali]   = tau_gbr[idx_gbr_aic2[2]]
  tau_gbr_BIC2[ali]   = tau_gbr[idx_gbr_bic2[2]]
  alpha_gbr_AIC2[ali] = alpha_gbr[idx_gbr_aic2[1]]
  alpha_gbr_BIC2[ali] = alpha_gbr[idx_gbr_bic2[1]]
  
  b_gbr_AIC2[,ali]   = b_gbr_temp[, idx_gbr_aic2[1], idx_gbr_aic2[2], ali]
  b_gbr_BIC2[,ali]   = b_gbr_temp[, idx_gbr_bic2[1], idx_gbr_bic2[2], ali]
}

table(alpha_gbr_AIC2)
table(alpha_gbr_BIC2)
table(tau_gbr_AIC2)
table(tau_gbr_BIC2)

# change
save(alpha_gbr_AIC2, tau_gbr_AIC2, alpha_gbr_BIC2, tau_gbr_BIC2, 
     b_gbr_AIC2, b_gbr_BIC2, file="bt3_NGR_n100_SNR2.RData")     

# plot result #
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


beta_hat_gbr_T1 = apply(b_gbr_BIC2, 2, FUN = betahat)
T1hat = apply(beta_hat_gbr_T1, 2, T1.hat)
median(T1hat)
mean(T1hat)
sd(T1hat)
par(mfrow = c(1,1))
hist(T1hat, probability = TRUE, xlab=expression(T[1]), main ="")
abline(v = 0.5, col = "red")


# estimate of delta (T1 in first version)
par(mfrow = c(1,1))
beta_hat_gbr  = apply(b_gbr_BIC2, 2, FUN = betahat)
plot(tobs, beta_true, type="l", ylim = c(-0.5, 2.5), ylab = expression(beta(t)), 
     xlab = "t", main = "BIC2")
lines(seq(0,1,length.out=101), beta_hat_gbr[,1], col="red")
matlines(tobs, beta_hat_gbr)
lines(tobs, beta_true, lwd = 3)