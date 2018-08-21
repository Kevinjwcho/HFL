# Method A
ChooseLambdas0 = function(qY,qbeta,sig,designs,psis,thetas,lambdas){
  
  # qY -- predictions from parametric representation for beta
  # qbeta -- parametric representation of beta
  # sig -- residual standard error
  # designs -- array of design matrices corresponding to different tehta
  # psis -- array of eigen decompositions for different thetas (should be zero after truncation point)
  # thetas -- thetas to use
  # lambdas -- set of lambda values to consider
  
  # vector to store expected MSE's for each theta
  beta0MSE = rep(0,length(thetas))
  Yhat0MSE = rep(0,length(thetas))
  
  beta0surf = matrix(0,nrow(psis),length(thetas))
  
  # Loop over values of theta
  for(i in 1:length(thetas)){
    
    # add intercept to design
    idesign = cbind( rep(1,dim(designs)[1]), designs[,,i])
    Psi0 = cbind( rep(1,dim(psis)[1]), psis[,,i])
    
    # regress on design matrix + extract coefficients
    beta0lm = lm(qY ~ designs[,,i])
    beta0coef = beta0lm$coef
    
    
    # Estimate beta
    beta0est =  (as.matrix(psis[,,i])%*%beta0coef[-1])
    
    # Expected mean squared error for beta
    beta0MSE[i] =  mean(  (beta0est - qbeta)^2 ) +
      sig*sum(diag( (ginv(t(idesign)%*%idesign)%*% 
                       t(Psi0)%*%Psi0)))/dim(psis)[1]
    
    # Expected predictive squared error 
    Yhat0MSE[i] = mean( ( beta0lm$fit - qY)^2 ) + sig*beta0lm$rank/ncol(X)
    
    # values of beta
    beta0surf[,i] = beta0est
  }
  
  
  # Now store errors for lambda
  b0MSE = rep(0,length(lambdas))
  b0MSE2 = rep(0,length(lambdas))
  thbest0 = rep(0,length(lambdas))
  thetavar0 = rep(0,length(lambdas))
  
  thindex = 1:length(thetas)
  
  # Look over lambdas
  for(j in 1:length(lambdas)){
    
    # Which is the best theta at this lambda
    thbest0[j] = which.min( Yhat0MSE + lambdas[j]*thetas^2 )
    
    
    # Tale secpmd differemces (assuming unit spacing of theta)
    
    if(thbest0[j] == 1){
      d0MSE =  (Yhat0MSE[thbest0[j]] + Yhat0MSE[thbest0[j]+2] - 2*Yhat0MSE[thbest0[j]+1])
    }else if(thbest0[j] == length(thetas)){
      d0MSE =  (Yhat0MSE[thbest0[j]-2] + Yhat0MSE[thbest0[j]] - 2*Yhat0MSE[thbest0[j]-1])
    }else{
      d0MSE =  (Yhat0MSE[thbest0[j]-1] + Yhat0MSE[thbest0[j]+1] - 2*Yhat0MSE[thbest0[j]])
    }
    
    # Laplace approximation to distribution of theta
    thetavar0[j] = (sig/ncol(X))*d0MSE/(d0MSE+lambdas[j])^2
    
    W0 = exp( - (thetas-thetas[thbest0[j]])^2/(2*thetavar0[j]) )
    
    
    # Take a weighted mean of MSE for estimating beta over distribution of thetas. 
    b0MSE[j] = weighted.mean(beta0MSE,W0,na.rm=TRUE)
  }
  
  # Which is the best lambda to use
  lbest0 = which.min(b0MSE)
  
  # And the corresponding theta
  thetabest0 = thbest0[lbest0]
  
  
  # return these, along with objective function, expected squared error for beta over lambdas, EMSE for Y,
  # and EMSE for beta of theta. 
  return( list(lbest = lbest0, thbest = thetabest0, obj = b0MSE[lbest0],b0MSE,Yhat0MSE,beta0MSE) )
}





ChooseLambdas1 = function(qY,qbeta,sig,X,psi1,thetas,lambdas){
  
  # qY -- predictions from parametric representation for beta
  # qbeta -- parametric representation of beta
  # sig -- residual standard error
  # X -- values of covariate
  # psi1 -- eigen decomposition on whole range (this will be truncated)
  # thetas -- thetas to use
  # lambdas -- set of lambda values to consider
  
  
  # Store some quantities over values of theta. 
  beta1MSE = rep(0,length(thetas))
  beta1bias = rep(0,length(thetas))
  beta1var = rep(0,length(thetas))
  Yhat1MSE = rep(0,length(thetas))
  
  # Designs (X*Psi1) for each level of truncation
  designs1 = array(0,c(ncol(X),ncol(psi1)+1,length(thetas)))
  for(i in 1:length(thetas)){
    designs1[,,i] = cbind( rep(1,ncol(X)), t(X[1:thetas[i],])%*%psi1[1:thetas[i],]/nrow(psi1))
  }
  
  # Length of theta vector
  nt = length(thetas)
  
  # Add either an intercept of explicitly not an intercept to Psi
  Psi1 = cbind(rep(1,nrow(psi1)),psi1)
  Psi01 = cbind(rep(0,nrow(psi1)),psi1)
  
  # Estiamte on whole range
  beta1lm = lm(qY ~ designs1[,,nt]-1)
  
  # Coefficients
  beta1coef = beta1lm$coef
  
  # Equivalent beta
  beta1est = psi1%*%beta1coef[-1]
  
  # Fitted values
  beta1pred = beta1lm$fit
  
  # Equivalent of X'*X for the design matrix
  D1 =  t(designs1[,,nt])%*%designs1[,,nt]
  
  times = 1:max(thetas)
  
  # Loop over thetas
  for(i in 1:length(thetas)){
    
    # How much error do I get for truncating?
    beta1bias[i] =  mean( ( beta1est*(times <= thetas[i]) - qbeta)^2 )
    
    # Variance at each level of truncation 
    beta1var[i] = sig*sum(diag( solve(D1,  t(Psi01[1:thetas[i],])%*%Psi01[1:thetas[i],])))/nrow(psi1)
    
    # Mean squared error
    beta1MSE[i] = mean( ( beta1est*(times <= thetas[i]) - qbeta)^2 ) +
              sig*sum(diag( solve(D1,  t(Psi01[1:thetas[i],])%*%Psi01[1:thetas[i],])))/nrow(psi1)
    
    # Predicted values at each level of truncation, but adjust the mean
    ibeta1pred = designs1[,,i]%*%beta1coef
    ibeta1pred = ibeta1pred - mean(ibeta1pred) + mean(beta1pred)
    
    # Expected predicted squared error
    Yhat1MSE[i] = mean( (ibeta1pred  - qY)^2 ) +
      sig*sum(diag( solve(D1, t(designs1[,,i])%*%designs1[,,i])))/ncol(X) + 
      sig*sum(diag( solve(D1, t(designs1[,,nt]-designs1[,,i])%*%matrix(1/ncol(X),ncol(X),ncol(X))%*%(designs1[,,nt]-designs1[,,i]))))/ncol(X)
  }
  
  
  # Now we'll look over values of lambda
  
  b1MSE = rep(0,length(lambdas))
  thbest1 = rep(0,length(lambdas))
  thetavar1 = rep(0,length(lambdas))
  
  thindex = 1:length(thetas)
  
  for(j in 1:length(lambdas)){
    
    # First difference of Yhat1MSE
    dYhat1MSE = c(diff(Yhat1MSE + lambdas[j]*thetas^2),1)
    
    # Find the first minimum
    thbest1[j] = thindex[dYhat1MSE > 0][1]
    
    # Now look at curvature  
    if(thbest1[j] == 1){
      d1MSE =  (Yhat1MSE[thbest1[j]] + Yhat1MSE[thbest1[j]+2] - 2*Yhat1MSE[thbest1[j]+1])
    }else if(thbest1[j] == length(thetas)){
      d1MSE =  (Yhat1MSE[thbest1[j]-2] + Yhat1MSE[thbest1[j]] - 2*Yhat1MSE[thbest1[j]-1])
    }else{
      d1MSE =  (Yhat1MSE[thbest1[j]-1] + Yhat1MSE[thbest1[j]+1] - 2*Yhat1MSE[thbest1[j]])
    }
    
    # Laplace approximation to distribution of theta at this lambda
    thetavar1[j] = (sig/ncol(X))*d1MSE/(d1MSE+lambdas[j])^2
    W1 = exp( - (thetas-thetas[thbest1[j]])^2/(2*thetavar1[j]) )
    
    # Expectation over this laplace approximation
    b1MSE[j] = weighted.mean(beta1MSE,W1,na.rm=TRUE)
  }
  
  # Which is the best lambda
  b1MSE[which(is.na(b1MSE))]=10^5
  lbest1 = which.min(b1MSE)
  
  # Which is the best theta
  thetabest1 = thbest1[lbest1]
  
  
  return( list(lbest = lbest1, thetabest = thbest1, b1MSE,Yhat1MSE,beta1MSE,beta1bias,beta1var,thbest1,thetavar1)  )
}





# fPCR search functions

fPCR = function(Y,X,psi,theta){
  # Conducts principal components regression truncating the eigenfunctions 
  # truncated at the value theta
  XX = t(X[1:theta,])%*%psi[1:theta,]/nrow(psi)
  mod = lm(Y~XX)
  return(mod)
}

fPCRsearch0 = function(Y,X,psi,theta.range = 1:nrow(X),lambda=0)
{
  # theta.range = range of theta values to examine
  # lambda =  penalty parameter
  
  # Searches over penalized fPCR to find first minimum in  theta for each lambdda
  
  pensse = rep(0,diff(theta.range)+1)
  for(th in theta.range[1]:theta.range[2]){
    mod = fPCR(Y,X,psi,th)
    pensse[th-theta.range[1]+1] = mean(mod$resid^2) + lambda*th^2
  }
  dpensse = c(diff(pensse),1)
  thetabest = (theta.range[1]:theta.range[2])[dpensse > 0][1] 
  # thetabest = which.min(pensse) + theta.range[1] - 1
  
  mod = fPCR(Y,X,psi,thetabest)
  return( list(mod = mod, thetabest = thetabest,pensse-lambda*th^2) )
  
}


fPCRsearch = function(Y,designs,psis,thetas,lambda=0)
{
  # theta.range = range of theta values to examine
  # lambda =  penalty parameter
  
  # Searches over penalized fPCR to find minimizing theta for each lambdda
  pensse = rep(0,length(thetas))
  for(i in 1:length(thetas)){
    mod = lm(Y~designs[,,i])
    pensse[i] = mean(mod$resid^2) + lambda*thetas[i]^2
  }
  
  thbest = which.min(pensse) 
  
  mod = lm(Y~designs[,,thbest])
  a = mod$coef[1]
  beta = as.matrix(psis[,,thbest])%*%mod$coef[-1]
  return( list(mod = mod, thetabest = thetas[thbest],a = a,beta=beta,pensse-lambda*thetas^2) )
  
}


fPCRsearch2 = function(Y,X,a,beta,theta.range = 1:nrow(X),lambda = 0){
  #  Takes an already estimated value of beta and searches over 
  #  possible truncation levels for it, returning the first minimum
  #  in the penalized sum of squares 
  
  # Note that a is an intercept, which we considered modifying for each 
  # truncation level, but found that the unmodified version is more 
  # numerically stable. 
  
  meanpred = mean(t(X)%*%beta)
  
  pensse = rep(0,diff(theta.range)+1)
  for(th in theta.range[1]:theta.range[2]){
    pred = t(X[1:th,])%*%beta[1:th]/length(beta)
    
    #   new.a =  a + mean(pred) - meanpred
    new.a = a
    allpred = pred + new.a
    
    pensse[th-theta.range[1]+1] = mean( (Y - allpred)^2) + lambda*th^2
  }
  dpensse = c(diff(pensse),1)
  thetabest = (theta.range[1]:theta.range[2])[dpensse > 0][1] 
  #  thetabest = which.min(pensse) + theta.range[1] -1
  beta[ thetabest:length(beta) ] = 0
  
  new.a = a + mean(t(X)%*%beta) - meanpred
  
  return( list(thetabest=thetabest,a=new.a,beta=beta,pensse-lambda*th^2) )
  
}
