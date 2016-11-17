kateVOL.k <- function(ret, k1=50, k2=1250, w=0.8, vol.min=0.001) {
  vol = pmax((w*volEWMA(ret,k1)+(1+w)*volEWMA(ret,k2)),vol.min)
  vol= as.numeric(vol)
  return(vol)
}

wsVOL.k <-function(ret,k1=50, k2=1250,vol.min=0.0001) {
  vol = pmax(sqrt(0.8*volEWMA(ret,k1)**2 + 0.2*volEWMA(ret,k2)**2),vol.min)
  vol=as.numeric(vol)
  return(vol)
}

volEWMA <- function(x, N=1000, init=var(x ,na.rm=TRUE)) {
  y= x^2
  y[is.na(x)] = var(x, na.rm=TRUE)
  lambda = 1-1/N
  z = filter((1-lambda)*y,filter=(lambda),method="recursive", init=init)
  zz=sqrt(z)
  return(zz)
}

EWMA <-function(x, days=22, init=mean(x[which(!is.na(x))[1:5]], na.rm=TRUE)) {
  y=x
  y[is.na(x)]=mean(x,na.rm=T)
  lambda=1-1/days
  z=filter((1-lambda)*y,filter=(lambda),method="recursive", init=init)
  z=as.numeric(z)
  return(z)
}
  
OscN <- function(x,N) {
  oN = EWMA(x,2**N)-EWMA(x,3*2**N)
  oN = as.numeric(oN)
  return(oN)
}
  
xN4 <-function(x) {
  xN4=x*exp(-(x**2)/4)
  return(xN4)
}
  
xN4_var <-function(x) {
  var=exp((x**2)/4)
  return(var)
}
  

