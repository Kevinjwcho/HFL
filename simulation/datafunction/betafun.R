beta_fun = function(t, ii)
{
  if(ii == 1) {bf = ifelse(t<=0.5 && t>=0,1,0)}
  else if(ii == 2){bf = sin(2*pi*t)*ifelse(t<=0.5 && t>=0,1,0)}
  else if(ii == 3){bf = (cos(2*pi*t)+1)*ifelse(t<=0.5 && t>=0,1,0)}
  else if(ii == 4){bf = (-100*(t-0.5)^3 - 200*(t - 0.5)^4)*ifelse(t<=0.5 && t>=0,1,0)}
  else{print("model does not exit")}
  return(bf)
}

# # plot betas
domain = c(0,1)
m=101
tobs <- seq(domain[1],domain[2],by=1/(m-1)) #observation grid
 for(i in 1:3)
 {
  betatobs =  apply(as.matrix(tobs), 1, beta_fun, ii = i)
  plot(tobs, betatobs, main = paste("Model", i), type = "l")
 }



t   = seq(0,1,by=0.0001)
# plot beta 1
i=1
beta_t =  apply(as.matrix(t), 1, beta_fun, ii = i)
par(mar=c(4, 4, 2, 0.5))
plot(t,beta_t, type="l", ylim = c(-0.2, 1), col = "white", 
     yaxt = 'n', ylab=expression(beta(t)), main="Model 1")
lines(t[1:5001],beta_t[1:5001])
lines(t[5002:10001],beta_t[5002:10001])
axis(2, at=seq(-0.2,1,by=0.2), labels=seq(-0.2,1,by=0.2),las=2)
abline(v=0.5,lty=2,col="red") 
 
# plot beta 2
i=2
beta_t =  apply(as.matrix(t), 1, beta_fun, ii = i)
par(mar=c(4, 4, 2, 0.5))
plot(t,beta_t, type="l", ylim = c(-0.2, 1),  
     yaxt = 'n', ylab=expression(beta(t)), main="Model 2")
axis(2, at=seq(-0.2,1,by=0.2), labels=seq(-0.2,1,by=0.2),las=2)
abline(v=0.5,lty=2,col="red") 

# plot beta 3
i=3
beta_t =  apply(as.matrix(t), 1, beta_fun, ii = i)
par(mar=c(4, 4, 2, 0.5))
plot(t,beta_t, type="l", ylim = c(-0.5, 2),  
     yaxt = 'n', ylab=expression(beta(t)), main="Model 3")
axis(2, at=seq(-0.4,2,by=0.4), labels=seq(-0.4,2,by=0.4),las=2)
abline(v=0.5,lty=2,col="red") 
