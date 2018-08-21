# define function for the estimator given alpha, tau, gamma and niter
# INPUT: 1,   alpha : tunning parameter for roughness penalty term
#        2,   tau   : tunning parameter related to group bridge penalty term
#        3,4, gamma : tunning parameter ralated to group bridge penalty term 
#        5,6, xuse, yuse : training U^** and Y^* 
#        7,8, xvalid, yvalid : validation U and Y 

# OUTPUT: 1, bhat: estimator
#         2, RSS : RSS of the validation sets
bhat_gbridge = function(b0, Y, u, v, n, alpha, tau, gamma, niter, M, d) 
{
  biter = matrix(NA, nrow = length(b0), ncol = niter)
  for(iterbg in 1:niter)
  {
    if(iterbg == 1) {biter[,iterbg] = b0}
    ##################################
    #    Step 1: compute theta_js    #
    ##################################
    else {
      theta = rep(0, M)
      h     = rep(0, M)
      for(j in 1:M)
      {
        theta[j] = theta_js(b_sminus1 = biter[,(iterbg-1)], cj = cj(M = M, d = d, gamma = gamma, j = j, bsmooth = b0), n = n, gamma = gamma, tau = tau, M = M, d = d, j = j)
        h[j]     = h_js(theta[j], cj(M = M, d = d, gamma = gamma, j = j, bsmooth = b0), gamma)
      }
      
      
      ##################################
      #    Step 2: compute g_ks        #
      ##################################  
      g = rep(0, (M + d)) 
      for(k in 1:(M + d))
      {
        g[k] = g_ks(h_js = h, k = k, M = M, d = d)
      }
      
      ##################################
      #    Step 3: compute bs          #
      ##################################  
      Ustarbhatgbr = Ustar(u = u, v = v, n = n, alpha = alpha, M = M, d = d)
      Ystarbhatgbr = Ystar(Y = Y, alpha = alpha, M = M, d = d)
      Ustarstargbr = Ustarbhatgbr%*%diag(1/(n*g))
      lassomodel   = glmnet(x = Ustarstargbr, y = Ystarbhatgbr, standardize = FALSE, alpha = 1, lambda = 0.5/n, family = "gaussian", intercept = FALSE)
      biter[,iterbg] =  coef(lassomodel)[2:length(coef(lassomodel)),]/(n*g)
      
      ################################################
      # decide whether to stop the iteration         #
      ################################################
     
      difratio = rep(0, length(biter[,iterbg]))
      if(iterbg >= 3){
        idx0 = which((biter[,iterbg]-biter[,(iterbg-1)]) == 0)
        if(length(idx0) == 0){difratio = (biter[,iterbg] - biter[, (iterbg-1)])/biter[,(iterbg-1)]}
        else {difratio[-idx0] = (biter[-idx0,iterbg] - biter[-idx0, (iterbg-1)])/biter[-idx0,(iterbg-1)]}
        if(max(difratio) < 10^(-4)) break}
      }
  }
  
  bhatgbr = biter[,iterbg] 
  return(list(bhatgbr = bhatgbr, gk = g))
}



# define function for the estimator given alpha = 0, tau, gamma and niter
# INPUT: 1. nothing
#        2,   tau   : tunning parameter related to group bridge penalty term
#        3,4, gamma : tunning parameter ralated to group bridge penalty term 
#        5,6, xuse, yuse : training U^** and Y^* 
#        7,8, xvalid, yvalid : validation U and Y 

# OUTPUT: 1, bhat: estimator
#         2, RSS : RSS of the validation sets
bhat_gbridge2 = function(b0, Y, u, n, tau, gamma, niter, M, d) 
{
  biter = matrix(NA, nrow = length(b0), ncol = niter)
  for(iterbg in 1:niter)
  {
    if(iterbg == 1) {biter[,iterbg] = b0}
    ##################################
    #    Step 1: compute theta_js    #
    ##################################
    else {
      theta = rep(0, M)
      h     = rep(0, M)
      for(j in 1:M)
      {
        theta[j] = theta_js(b_sminus1 = biter[,(iterbg-1)], cj = cj(M = M, d = d, gamma = gamma, j = j, bsmooth = b0), n = n, gamma = gamma, tau = tau, M = M, d = d, j = j)
        h[j]     = h_js(theta[j], cj(M = M, d = d, gamma = gamma, j = j, bsmooth = b0), gamma)
      }
      
      
      ##################################
      #    Step 2: compute g_ks        #
      ##################################  
      g = rep(0, (M + d)) 
      for(k in 1:(M + d))
      {
        g[k] = g_ks(h_js = h, k = k, M = M, d = d)
      }
      
      ##################################
      #    Step 3: compute bs          #
      ##################################  
      Ustarstargbr = u%*%diag(1/(n*g))
      lassomodel   = glmnet(x = Ustarstargbr, y = Y, standardize = FALSE, alpha = 1, lambda = 0.5/n, family = "gaussian", intercept = FALSE)
      biter[,iterbg] =  coef(lassomodel)[2:length(coef(lassomodel)),]/(n*g)
      
      ################################################
      # decide whether to stop the iteration         #
      ################################################
      
      difratio = rep(0, length(biter[,iterbg]))
      if(iterbg >= 3){
        idx0 = which((biter[,iterbg]-biter[,(iterbg-1)]) == 0)
        if(length(idx0) == 0){difratio = (biter[,iterbg] - biter[, (iterbg-1)])/biter[,(iterbg-1)]}
        else {difratio[-idx0] = (biter[-idx0,iterbg] - biter[-idx0, (iterbg-1)])/biter[-idx0,(iterbg-1)]}
        if(max(difratio) < 10^(-4)) break}
    }
  }
  
  bhatgbr = biter[,iterbg] 
  return(list(bhatgbr = bhatgbr, gk = g))
}