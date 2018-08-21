# function for smooth estimator (with roughness penalty)
smoothb = function(Y, U, V, n, alpha, M, d)
{ 
  Usmoothb  = Ustar(u = U, v = v, n = n, alpha = alpha, M = M, d = d)
  Ysmoothb  = Ystar(Y = Y, alpha = alpha, M = M, d = d)
  smoothfit = lm(Ysmoothb ~  0 + Usmoothb)
  return(bmoothhat = smoothfit$coefficients)
}

