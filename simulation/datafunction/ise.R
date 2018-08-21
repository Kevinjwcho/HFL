# ISE
ise = function(b,b0,domain) # input b and b0 is a vector, domain is c(start, end) and is c(0, 1) here
{
  lgb = length(b) 
  difb2=(b-b0)^2
  if(lgb%%2 == 1)
  {h = (domain[2]-domain[1])/(lgb-1)
  cef = c(1, rep(c(4,2), (lgb-3)/2), 4, 1)
  isetp = h/3*cef%*%difb2 
  }
  else
  {
    h = (domain[2]-domain[1])/lgb
    isetp = h*sum(difb2) }
  return(isetp)
}