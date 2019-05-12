###################################################
### code chunk number 34: eco2.rnw:1100-1130
###################################################
per<-function(x,plot_T)
{
  len<-length(x)
  per<-0:(len/2)
  DFT<-per
  
  for (k in 0:(len/2))
  {
    cexp <- complex(arg=-(1:len)*2*pi*k/len)
    DFT[k+1]<-sum(cexp*x*sqrt(1/(2*pi*len)))
  }
# Weighths wk: if length of data sample is even then DFT in frequency pi is scaled by 1/sqrt(2) (Periodogram in pi is weighted by 1/2)  
  if (abs(as.integer(len/2)-len/2)<0.1)
    DFT[k+1]<-DFT[k+1]/sqrt(2)
  per<-abs(DFT)^2
  if (plot_T)
  {
    par(mfrow=c(2,1))
    plot(per,type="l",axes=F,xlab="Frequency",ylab="Periodogram",
         main="Periodogram")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
                                    "4pi/6","5pi/6","pi"))
    axis(2)
    box()
    plot(log(per),type="l",axes=F,xlab="Frequency",ylab="Log-periodogram",
         main="Log-periodogram")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
                                    "4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }
  return(list(DFT=DFT,per=per))
}








###################################################
### code chunk number 72: eco2.rnw:2846-2879
###################################################
# This function computes mean-square DFA-solutions
# L is the length of the MA filter,
# periodogram is the frequency weighting function in the DFA
# Gamma is the transferfunction of the symmetric filter (target) and
# Lag is the lag-parameter: Lag=0 implies real-time filtering, Lag=L/2
#     implies symmetric filter
# The function returns optimal coefficients as well as the transfer function of the
#     optimized real-time filter
#
#
dfa_ms<-function(L,periodogram,Lag,Gamma)
{
  
  K<-length(periodogram)-1
  X<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)*sqrt(periodogram)
  X_y<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)
  for (l in 2:L)          #l<-L<-21
  {
    X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                  1.i*sin((l-1-Lag)*pi*(0:(K))/(K)))*sqrt(periodogram))
    X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                      1.i*sin((l-1-Lag)*pi*(0:(K))/(K))))
  }
  xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
  # MA-Filtercoefficients
  b<-as.vector(solve(xtx)%*%(t(Re(X_y))%*%(Gamma*periodogram)))
  # Transferfunction
  trffkt<-1:(K+1)
  trffkt[1]<-sum(b)
  for (k in 1:(K))#k<-1
  {
    trffkt[k+1]<-(b%*%exp(1.i*k*(0:(length(b)-1))*pi/(K)))
  }
  return(list(b=b,trffkt=trffkt))
}






###################################################
### code chunk number 101: eco2.rnw:4148-4218
###################################################
# This function computes analytic DFA-solutions
# L is the length of the MA-filter,
# weight_func is the periodogram,
# lambda emphasizes phase artifacts in the customized criterion,
# eta emphasizes noise-suppression/smoothness
# Gamma is the transferfunction of the symmetric filter (target) and
# Lag is the lag-parameter: Lag=0 implies real-time filtering, Lag=L/2
#     implies symmetric filter
# i1 and i2 allow for filter restrictions in frequency zero
# The function returns the weights of the MA-Filter as well as its transferfunction
#
#
dfa_analytic<-function(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2)
{
  K<-length(weight_func)-1
  # Define the amplitude weighting function weight_h (W(omega_k))
  omega_Gamma<-as.integer(cutoff*K/pi)
  if ((K-omega_Gamma+1)>0)
  {
    weight_h<-weight_func*(c(rep(1,omega_Gamma),(1:(K-omega_Gamma+1))^(eta)))
  } else
  {
    weight_h<-weight_func*rep(1,K+1)
  }
  # First order filter restriction: assigning a `large' weight to frequency zero
  if (i1)
    weight_h[1]<-max(1.e+10,weight_h[1])
  
  X<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)*sqrt(weight_h)
  X_y<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)
  if (i2)
  {
    # Second order restriction: time shift in frequency zero vanishes
    for (l in 2:(L-1))
    {
      X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
                    cos((L-1-Lag)*pi*(0:(K))/(K))+
                    sqrt(1+Gamma*lambda)*1.i*(sin((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
                                                sin((L-1-Lag)*pi*(0:(K))/(K))))*sqrt(weight_h))
      X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
                        cos((L-1-Lag)*pi*(0:(K))/(K))+
                        1.i*(sin((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*sin((L-1-Lag)*pi*(0:(K))/(K)))))
    }
    xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
    # MA-Filterweights
    b<-as.vector(solve(xtx)%*%(t(Re(X_y))%*%(Gamma*weight_h)))
    # the last weight is a function of the previous ones through the second order restriction
    b<-c(b,-sum(b*(0:(length(b)-1)))/(length(b)))
  } else
  {
    for (l in 2:L)
    {
      X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                    sqrt(1+Gamma*lambda)*1.i*sin((l-1-Lag)*pi*(0:(K))/(K)))*sqrt(weight_h))
      X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                        1.i*sin((l-1-Lag)*pi*(0:(K))/(K))))
    }
    xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
    # MA-Filterweights
    b<-as.vector(solve(xtx)%*%(t(Re(X_y))%*%(Gamma*weight_h)))
  }
  # Transferfunction
  trffkt<-1:(K+1)
  trffkt[1]<-sum(b)
  for (k in 1:(K))#k<-1
  {
    trffkt[k+1]<-(b%*%exp(1.i*k*(0:(length(b)-1))*pi/(K)))
  }
  return(list(b=b,trffkt=trffkt))
}
