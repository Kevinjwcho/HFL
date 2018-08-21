# smooth method AIC and BIC
library(psych)
AICBIC.smooth.fun = function(Y, b, U, V, n, alpha)
{
  hat.s  = U%*%solve(t(U)%*%U + n*alpha*V)%*%t(U)
  df.s   = tr(hat.s)
  Yhat.s = U%*%b
  RSS.s  = t(Y - Yhat.s)%*%(Y - Yhat.s)
  AICtemp.s = n*log(RSS.s/n) + 2*df.s
  BICtemp.s = n*log(RSS.s/n) + log(n)*df.s
  return(list(AIC = AICtemp.s, BIC = BICtemp.s))
}



# grouop bridge AIC and BIC
AICBIC.gbr.fun = function(Y, b, U, V, W, n, alpha)
{
  sparse.idx   = which(b == 0)
  if(length(sparse.idx) == 0)
  {
    ula = U
    vla = v
    wla = W
  }
  else{
  ula  = U[, -sparse.idx]
  vla  = V[-sparse.idx, -sparse.idx]
  wla  = W[-sparse.idx, -sparse.idx]
  }
  #hat1 = ula%*%solve(t(ula)%*%ula + n*alpha*vla + 0.5*wla)%*%t(ula)
  hat2 = ula%*%solve(t(ula)%*%ula + n*alpha*vla)%*%t(ula)
  #df1  = tr(hat1)
  df2  = tr(hat2)
  Yhat = U%*%b
  RSS.gbric  = t(Y - Yhat)%*%(Y - Yhat)
  #AIC1temp.gbr = n*log(RSS.gbric/n) + 2*df1
  AIC2temp.gbr = n*log(RSS.gbric/n) + 2*df2
  #BIC1temp.gbr = n*log(RSS.gbric/n) + log(n)*df1
  BIC2temp.gbr = n*log(RSS.gbric/n) + log(n)*df2
  return(list(AIC2 = AIC2temp.gbr, BIC2 = BIC2temp.gbr))
}
