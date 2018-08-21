# function h and g
h_js = function(theta_js, cj, gamma)
{
  temphjs = theta_js^(1-1/gamma)*cj^(1/gamma)
  return(temphjs)
}
  
g_ks = function(h_js, k, M, d)
{
  if(k <= M) {tempgks = sum(h_js[1:k])}
  else {tempgks = sum(h_js[1:M])}
}