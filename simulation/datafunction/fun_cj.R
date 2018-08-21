# theta_j^(s) for iteratively update theta
# A_j, j =1, ..., 100
# bsmooth is a vector with the number of entries = nbasis (b_1, b_2, ... b_(M+d))

cj = function(M, d, gamma, j, bsmooth)
{
  tempcj  = (d + M + 1 - j)^(1-gamma)
  bsm_Aj_norm_gamma  = sqrt(sum((bsmooth[j:(M+d)])^2))^gamma
  tempcjfinal = tempcj/bsm_Aj_norm_gamma
  return(tempcjfinal)
}

