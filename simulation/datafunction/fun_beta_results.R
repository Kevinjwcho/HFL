# Beta hat, obtained from b #
betahat =  function(b, basismat = basismat.tobs)
{
  temp = t(as.matrix(b))%*%t(basismat)
  return(temp)
}

# T1 estimate
# input "beta" is a vector containing beta(t_i) values
T1.hat = function(beta)
{
  T1n = length(beta)
  t1hatfuntemp = 1
  for(t1i in 1:(T1n-1))
  {
    if(beta[t1i]!=0 && sum(abs(beta)[(t1i+1):T1n])==0)
      {t1hatfuntemp = t1i/(T1n-1)}
  }
return(t1hatfuntemp)
}




# bias of beta: betahat - beta #
pointwise.beta.bias =  function(beta.hat, beta.true = beta_true, n.data = ndataset)
{
  beta.true.matrix = matrix(rep(beta.true, n.data), ncol = n.data)
  bias.matrix      = beta.hat - beta.true.matrix
  temp             = apply(bias.matrix, 1, mean)
  return(temp)
}


# difference of y and yhat #
dif.y =  function(y.hat, y.true)
{
  square.dif.matrix = (y.hat - y.true)^2
  temp              = apply(square.dif.matrix, 2, mean)
  return(temp)
}


# ISE
ISE.beta.all = function(bias)
{
  temp   = h/3*sum(cef*bias)
  return(temp)
}

# ISE01: for null region and non-null region seperated
ISE.beta.01 = function(bias)
{
  temp   = h/3*sum(cef01*bias)
  return(temp)
}

# ISE0: for null region and non-null region seperated
ISE.beta.0 = function(bias)
{
  temp   = h*sum(cef0*bias)
  return(temp)
}


# ISE1: for null region and non-null region seperated
ISE.beta.1 = function(bias)
{
  temp   = h*sum(cef1*bias)
  return(temp)
}
