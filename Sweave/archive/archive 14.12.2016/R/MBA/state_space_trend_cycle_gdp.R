
ss_model_i1<-function(lgdp)
{
  # I(1)-model
  
  # set up SS model: (Morley, Nelson, and Zivot, 2003): http://research.economics.unsw.edu.au/jmorley/mnz03.pdf
  #   Does not work if latest data is accounted for (great recession): AR is non-stationary and cycle is completely bogus
  #   Seems to work if Data up to Dec 2007 is used
  #   However cycle length is impossible: 7 quarters...
  # set up SS model
  ssm2 <- function(parm){
    dlm <- dlmModPoly(2,dV=1e-7,dW=c(parm[4]^2,0)) + 
      dlmModARMA(ar=c(parm[1],parm[2]), ma=NULL, sigma2=parm[3]^2)
    # get distribution variance of initial state
    tmp0 <- matrix(c(parm[1],parm[2],1,0),nr=2)
    tmp1 <- matrix(c(parm[3]^2,0,0,0),nc=1)
    tmp <- solve(diag(4)-tmp0%x%tmp0)%*%tmp1
    dlm$C0[3:4,3:4] <- matrix(tmp,nr=2)
    return( dlm )
  }
  
  # estimate parameters
  # use estimates in Morley, Nelson, Zivot (2003) as initializations
  fit2 <- dlmMLE(y=lgdp,parm=c(1.53,-.61,.62,.81),build=ssm2,hessian=T)
  mat_est<-c(fit2$value,abs(polyroot(c(-fit2$par[2:1],1))),abs((2*pi/Arg(polyroot(c(-fit2$par[2:1],1))))/4)[1],fit2$par[3:4])
  mat_parm<-fit2$par
  # use proxy in Fossati (2013): https://www.ualberta.ca/~sfossati/e509/files/slides/lec5.r
  fit2 <- dlmMLE(y=lgdp,parm=c(1.5,-.6,.4,.3),build=ssm2,hessian=T)
  mat_est<-rbind(mat_est,c(fit2$value,abs(polyroot(c(-fit2$par[2:1],1))),abs((2*pi/Arg(polyroot(c(-fit2$par[2:1],1))))/4)[1],fit2$par[3:4]))
  mat_parm<-rbind(mat_parm,fit2$par)
  # use estimates of best model below as initialization
  fit2 <- dlmMLE(y=lgdp,parm=c(1.3328158121,-0.403844596,0.8889332691,0.0001191306),build=ssm2,hessian=T)
  mat_est<-rbind(mat_est,c(fit2$value,abs(polyroot(c(-fit2$par[2:1],1))),abs((2*pi/Arg(polyroot(c(-fit2$par[2:1],1))))/4)[1],fit2$par[3:4]))
  mat_parm<-rbind(mat_parm,fit2$par)
  
  dimnames(mat_est)[[1]]<-c("Morley,Nelson,Zivot","Fossati","Own")
  dimnames(mat_est)[[2]]<-c("Criterion","rho 1" ,"rho 2","Duration (in years)","Sigma cycle","Sigma Trend")
  
  # Likelihood
  fit2$value
  parm<-fit2$par
  
  # Persistence and periodicity 
  abs(polyroot(c(-fit2$par[2:1],1)))
  (2*pi/Arg(polyroot(c(-fit2$par[2:1],1))))/4
  (2*pi/Arg(polyroot(c(1,-fit2$par[1:2]))))/4
  
  return(list(mat_est=mat_est,mat_parm=mat_parm))
}









ss_model_i2<-function(lgdp)
{
  
  # I(2)-model: it estimates a variable drift but no level innovation 
  ssm2 <- function(parm){
    dlm <- dlmModPoly(2,dV=1e-7,dW=c(0,parm[4]^2)) + 
      dlmModARMA(ar=c(parm[1],parm[2]), ma=NULL, sigma2=parm[3]^2)
    # get distribution variance of initial state
    tmp0 <- matrix(c(parm[1],parm[2],1,0),nr=2)
    tmp1 <- matrix(c(parm[3]^2,0,0,0),nc=1)
    tmp <- solve(diag(4)-tmp0%x%tmp0)%*%tmp1
    dlm$C0[3:4,3:4] <- matrix(tmp,nr=2)
    return( dlm )
  }
  
  # Estimate parameters: use simple initialization of parameters
  a_1<-2*0.89963*cos(2*pi/18.74201)
  a_2<--(0.89963^2)
  abs(polyroot(c(-a_2,-a_1,1)))
  2*pi/Arg(polyroot(c(-a_2,-a_1,1)))
  
  #   Likelihood is larger than I(1)-model, Innovations of trend+drift ~0: cycle takes it all (it is non-stationary); cycles looks bad;
  #   roots are real i.e. frequency=infty; 
  fit2 <- dlmMLE(y=lgdp,parm=c(1.5,-.6,sqrt(.4),sqrt(.3),sqrt(0.01)),build=ssm2,hessian=T)#,method="SANN")
  
#  sigma_eta<-0.01
#  fit2 <- dlmMLE(y=lgdp,parm=c(a_1,a_2,.5,sigma_eta),build=ssm2,hessian=T)
  mat_est<-c(fit2$value,abs(polyroot(c(-fit2$par[2:1],1))),abs((2*pi/Arg(polyroot(c(-fit2$par[2:1],1))))/4),fit2$par[3],0,fit2$par[4])
  mat_parm<-c(fit2$par[1:3],0,fit2$par[4])
  parm_best<-fit2$par
  
#------------------------------------  
# I(2)-model with imposed cycle, no level innovation 
  ssm2 <- function(parm){
    dlm <- dlmModPoly(2,dV=1e-7,dW=c(0,parm[3]^2)) + 
      dlmModARMA(ar=c(2*parm[1]*cos(2*pi/(6*4)),-(parm[1])^2), ma=NULL, sigma2=parm[2]^2)
    # get distribution variance of initial state
    tmp0 <- matrix(c(c(2*parm[1]*cos(2*pi/(6*4)),-(parm[1])^2),1,0),nr=2)
    tmp1 <- matrix(c(parm[2]^2,0,0,0),nc=1)
    tmp <- solve(diag(4)-tmp0%x%tmp0)%*%tmp1
    dlm$C0[3:4,3:4] <- matrix(tmp,nr=2)
    return( dlm )
  }
  
# Estimate parameters: use simple initialization of parameters
  a_1<-0.89963

#   Likelihood is larger than I(1)-model, Innovations of trend+drift ~0: cycle takes it all (it is non-stationary); cycles looks bad;
#   roots are real i.e. frequency=infty; 
  fit2 <- dlmMLE(y=lgdp,parm=c(a_1,sqrt(.4),sqrt(0.01)),build=ssm2,hessian=T)#,method="SANN")

#  sigma_eta<-0.01
#  fit2 <- dlmMLE(y=lgdp,parm=c(a_1,0,.5,sigma_eta),build=ssm2,hessian=T)
  mat_est<-rbind(mat_est,c(fit2$value,rep(fit2$par[1],2),rep(6,2),fit2$par[2],0,fit2$par[3]))
  mat_parm<-rbind(mat_parm,c(2*cos(2*pi/(4*6))*fit2$par[1],-(fit2$par[1]^2),fit2$par[2],0,fit2$par[3]))
#  parm_best<-fit2$par

#------------------------------------  

# I(1) + I(2) model
  
  ssm2 <- function(parm){
    dlm <- dlmModPoly(2,dV=1e-7,dW=c(parm[4]^2,parm[5]^2)) + 
      dlmModARMA(ar=c(parm[1],parm[2]), ma=NULL, sigma2=parm[3]^2)
    # get distribution variance of initial state
    tmp0 <- matrix(c(parm[1],parm[2],1,0),nr=2)
    tmp1 <- matrix(c(parm[3]^2,0,0,0),nc=1)
    tmp <- solve(diag(4)-tmp0%x%tmp0)%*%tmp1
    dlm$C0[3:4,3:4] <- matrix(tmp,nr=2)
    return( dlm )
  }
  
# Initialization is previous I(2)-model. Therefore the resulting estimate must be `at least as good'  
  fit2 <- dlmMLE(y=lgdp,parm=c(1.5,-.6,sqrt(.4),sqrt(.3),sqrt(0.01)),build=ssm2,hessian=T)#,method="SANN")
# fit2 <- dlmMLE(y=lgdp,parm=c(parm_best[1:3],0.,parm_best[4]),build=ssm2,hessian=T)
  mat_est<-rbind(mat_est,c(fit2$value,abs(polyroot(c(-fit2$par[2:1],1))),abs((2*pi/Arg(polyroot(c(-fit2$par[2:1],1))))/4),fit2$par[3:5]))
  mat_parm<-rbind(mat_parm,fit2$par)

#------------------------------------
# I(1) + I(2) model with imposed cycle-frequency

  ssm2 <- function(parm){
  dlm <- dlmModPoly(2,dV=1e-7,dW=c(parm[3]^2,parm[4]^2)) + 
    dlmModARMA(ar=c(2*parm[1]*cos(2*pi/(4*6)),-(parm[1])^2), ma=NULL, sigma2=parm[2]^2)
  # get distribution variance of initial state
  tmp0 <- matrix(c(c(2*parm[1]*cos(2*pi/(4*6)),-(parm[1])^2),1,0),nr=2)
  tmp1 <- matrix(c(parm[2]^2,0,0,0),nc=1)
  tmp <- solve(diag(4)-tmp0%x%tmp0)%*%tmp1
  dlm$C0[3:4,3:4] <- matrix(tmp,nr=2)
  return( dlm )
  }

  fit2 <- dlmMLE(y=lgdp,parm=c(0.89,sqrt(.4),sqrt(.3),sqrt(0.01)),build=ssm2,hessian=T)#,method="SANN")
#  fit2 <- dlmMLE(y=lgdp,parm=c(parm_best[1:3],0.3,0.6),build=ssm2,hessian=T)
  mat_est<-rbind(mat_est,c(fit2$value,rep(fit2$par[1],2),rep(6,2),fit2$par[2:4]))
  mat_parm<-rbind(mat_parm,c(2*cos(2*pi/(4*6))*fit2$par[1],-(fit2$par[1]^2),fit2$par[2:4]))

  
  dimnames(mat_est)[[1]]<-c("s_11=0","s_11=0, length=6","Unconstrained","Length=6")
  dimnames(mat_est)[[2]]<-c("Criterion","rho","rho","Duration","Duration","Sigma cycle","Sigma trend","Sigma drift")
  dimnames(mat_parm)[[1]]<-dimnames(mat_est)[[1]]

  
  return(list(mat_est=mat_est,mat_parm=mat_parm))
}






plot_ss_i1<-function(parm,lgdp,start_year)
{
  nobs<-length(lgdp)
  ssm2 <- function(parm){
    dlm <- dlmModPoly(2,dV=1e-7,dW=c(parm[4],0)) + 
      dlmModARMA(ar=c(parm[1],parm[2]), ma=NULL, sigma2=parm[3])
    # get distribution variance of initial state
    tmp0 <- matrix(c(parm[1],parm[2],1,0),nr=2)
    tmp1 <- matrix(c(parm[3],0,0,0),nc=1)
    tmp <- solve(diag(4)-tmp0%x%tmp0)%*%tmp1
    dlm$C0[3:4,3:4] <- matrix(tmp,nr=2)
    return( dlm )
  }
  
  mod2 <- ssm2(parm) 
  mod2f <- dlmFilter(lgdp,mod2)
  
  
  # get parameter estimates 
  drift <- mod2f$m[nobs+1,2]
  covar <- dlmSvd2var(mod2f$U.C[[nobs+1]],mod2f$D.C[nobs+1,])
  coef.se <- sqrt(covar[2,2])
  
  
  
  # filtered values: remove first initializations
  xtfilt <- ts(mod2f$m[-1,1],start=start_year,frequency=4)
  ctfilt <- ts(mod2f$m[-1,c(2,3)],start=start_year,frequency=4)
  # Remove first values because of initialization issues
  ctfilt[1:3,]<-NA
  # plot: 
  # Upper graph: data, filtered state (trend and cycle), 
  # Bottom graph: cycle and drift
  par(mfrow=c(1,2))
  mplot<-cbind(lgdp,xtfilt)
  ymin<-min(apply(mplot[-1,],1,min))
  ymax<-max(apply(mplot[-1,],1,max))#dim(mplot)
  plot(mplot,ylim=c(ymin,ymax),xlim=c(start_year,end_year),
       plot.type='s',col=c("black","blue"),ylab="",main="Log Real US GDP ")
  nberShade()
  mtext("GDP", side = 3, line = -1,at=mean(c(start_year,end_year)),col="black")
  mtext("Cycle-adjusted Component", side = 3, line = -2,at=mean(c(start_year,end_year)),col="blue")
  lines(lgdp)
  lines(xtfilt,col="blue")
  
  plot(ctfilt[,1],plot.type='s',ylim=c(min(ctfilt,na.rm=T),max(ctfilt,na.rm=T)),ylab="",main="Drift and Cycle I(1)-Model")
  nberShade()
  lines(ctfilt[,1])
  lines(ctfilt[,2],col="blue")
  mtext("drift", side = 3, line = -1,at=mean(c(start_year,end_year)),col="black")
  mtext("Cycle", side = 3, line = -2,at=mean(c(start_year,end_year)),col="blue")
  abline(h=0)
  
  
  
}




plot_ss_i2<-function(parm,lgdp,start_year,plot_which,title_main)
{
  nobs<-length(lgdp)
  ssm2 <- function(parm){
    dlm <- dlmModPoly(2,dV=1e-7,dW=c(parm[4],parm[5])) + 
      dlmModARMA(ar=c(parm[1],parm[2]), ma=NULL, sigma2=parm[3])
    # get distribution variance of initial state
    tmp0 <- matrix(c(parm[1],parm[2],1,0),nr=2)
    tmp1 <- matrix(c(parm[3],0,0,0),nc=1)
    tmp <- solve(diag(4)-tmp0%x%tmp0)%*%tmp1
    dlm$C0[3:4,3:4] <- matrix(tmp,nr=2)
    return( dlm )
  }
  
  mod2 <- ssm2(parm) 
  mod2f <- dlmFilter(lgdp,mod2)
  
  
  # get parameter estimates 
  drift <- mod2f$m[nobs+1,2]
  covar <- dlmSvd2var(mod2f$U.C[[nobs+1]],mod2f$D.C[nobs+1,])
  coef.se <- sqrt(covar[2,2])
  
  
  
  # filtered values: remove first initializations
  xtfilt <- ts(mod2f$m[-1,1],start=start_year,frequency=4)
  ctfilt <- ts(mod2f$m[-1,c(2,3)],start=start_year,frequency=4)
  # Remove first values because of initialization issues
  ctfilt[1:3,]<-NA
  # plot: 
  # Upper graph: data, filtered state (trend and cycle), 
  # Bottom graph: cycle and drift
  if (sum(plot_which)==2)
  {
    par(mfrow=c(1,2))
    mplot<-cbind(lgdp,xtfilt)
    ymin<-min(apply(mplot[-1,],1,min))
    ymax<-max(apply(mplot[-1,],1,max))#dim(mplot)
    plot(mplot,ylim=c(ymin,ymax),xlim=c(start_year,end_year),
         plot.type='s',col=c("black","blue"),ylab="",main="Log Real US GDP ")
    nberShade()
    mtext("GDP", side = 3, line = -1,at=mean(c(start_year,end_year)),col="black")
    mtext("Cycle-adjusted Component", side = 3, line = -2,at=mean(c(start_year,end_year)),col="blue")
    lines(lgdp)
    lines(xtfilt,col="blue")
    
    plot(ctfilt[,1],plot.type='s',ylim=c(min(ctfilt,na.rm=T),max(ctfilt,na.rm=T)),ylab="",main="Drift and Cycle I(2)-model")
    nberShade()
    lines(ctfilt[,1])
    lines(ctfilt[,2],col="blue")
    mtext("drift", side = 3, line = -1,at=mean(c(start_year,end_year)),col="black")
    mtext("Cycle", side = 3, line = -2,at=mean(c(start_year,end_year)),col="blue")
    abline(h=0)
  }  else
  if (plot_which[1]==T&plot_which[2]==F)
  {
    mplot<-cbind(lgdp,xtfilt)
    ymin<-min(apply(mplot[-1,],1,min))
    ymax<-max(apply(mplot[-1,],1,max))#dim(mplot)
    plot(mplot,ylim=c(ymin,ymax),xlim=c(start_year,end_year),
         plot.type='s',col=c("black","blue"),ylab="",main="Log Real US GDP ")
    nberShade()
    mtext("GDP", side = 3, line = -1,at=mean(c(start_year,end_year)),col="black")
    mtext("Cycle-adjusted Component", side = 3, line = -2,at=mean(c(start_year,end_year)),col="blue")
    lines(lgdp)
    lines(xtfilt,col="blue")    
  }  
  if (plot_which[1]==F&plot_which[2]==T)
  {
    
    plot(ctfilt[,1],plot.type='s',ylim=c(min(ctfilt,na.rm=T),max(ctfilt,na.rm=T)),ylab="",main=title_main)
    nberShade()
    lines(ctfilt[,1])
    lines(ctfilt[,2],col="blue")
    mtext("drift", side = 3, line = -1,at=mean(c(start_year,end_year)),col="black")
    mtext("Cycle", side = 3, line = -2,at=mean(c(start_year,end_year)),col="blue")
    abline(h=0)
  }
  return(list(ctfilt=ctfilt))  
}

