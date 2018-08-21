# Y^*, U^*
# Y^*
Ystar = function(Y, alpha, M, d)
{if(alpha == 0){tempystar = Y}else{tempystar = as.matrix(rbind(as.matrix(Y, ncol = 1), matrix(0, nrow = M + d, ncol = 1)))}
  return(tempystar)
}


# U^*
Ustar = function(u, v, n, alpha, M, d)
{ 
  if(alpha==0){tempustar = u}else{eig = eigen(v) 
  eig$values[eig$values < 0] = 0 
  w    = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
  tempustar = rbind(u, sqrt(n*alpha)*w)}
  return(tempustar)
}