# theta_j^(s) for iteratively update theta
# c_j(d, gamma, j), d is the degree = order -1, gamma is the tunning parameter
# j is the from 1:M (M is the number of subintervals divided by knots)

theta_js = function(b_sminus1, cj, n, gamma, tau, M, d, j)
{ 
  b2bargamma = (sum(abs(b_sminus1)[j:(d+M)]))^gamma 
  tempthetajs = cj*((1/gamma-1)/tau)^gamma*b2bargamma
  return(tempthetajs)
}