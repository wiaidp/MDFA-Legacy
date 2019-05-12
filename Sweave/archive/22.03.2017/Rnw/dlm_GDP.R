#https://www.ualberta.ca/~sfossati/e509/files/slides/lec5.r

# UC-0 (Morley, Nelson, and Zivot, 2003)
library(dlm)
# load numDeriv package 
library(numDeriv)
# recession shading
library(tis)  
# Quandl
# Load required libraries
library(RCurl)    # For getURL() and curl handler / cookie / google login
library(stringr)  # For str_trim() to trip whitespace from strings
library(Quandl)
library(xts)
require (quantmod)
require (Quandl)

#????????????????????????
Quandl.auth("yWbLG4aMd_UZAzCPyXe5")





source(file=paste(path_MBA.pgm,"state_space_trend_cycle_gdp.r",sep=""))



#-----------------------------------------------------------------------------
# Read data
# Last data point
end_date<-format(Sys.time(), "%Y-%m-%d")

# Data up to great recession
end_date<-"2007-12-07"

start_year<-1947
end_year<-as.double(substr(end_date,1,4))
start_date=paste(start_year,'-01-01',sep="")

#mydata<-Quandl(c('FRED/GDP'),start_date=start_date,end_date=end_date,type='xts')
mydata<-Quandl(c('FRED/GDPC96'),start_date=start_date,end_date=end_date,type='xts')

# Annualized sharpe of GDP series
sharpe_GDP<-sqrt(4)*mean(diff(lgdp))/sqrt(var(diff(lgdp)))
sharpe_GDP

mydata<-na.omit(mydata)
lgdp <- ts(100*log(mydata),start=start_year,frequency=4)
nobs <- length(lgdp)

ss_obj_i1<-ss_model_i1(lgdp)

mat_est_i1<-ss_obj_i1$mat_est
mat_parm_i1<-ss_obj_i1$mat_parm

# All three plots are almost indistinguishable (although cycles-lengths differ)
plot_ss_i1(mat_parm_i1[1,],lgdp,start_year)


ss_obj_i2<-ss_model_i2(lgdp)

mat_est_i2<-ss_obj_i2$mat_est
mat_parm_i2<-ss_obj_i2$mat_parm

# First and second models are identical i.e. level innovation is not used

# Cycle of first (or second) model is poor. 
par(mfrow=c(2,2))
for (i in 1:nrow(mat_parm_i2))
{
  plot_which<-c(F,T)
  title_main<-dimnames(mat_parm_i2)[[1]][i]
  plot_ss_i2(mat_parm_i2[i,],lgdp,start_year,plot_which,title_main)
}



#-------------------------------------------------------------------------

end_date<-format(Sys.time(), "%Y-%m-%d")


start_year<-1947
end_year<-as.double(substr(end_date,1,4))
start_date=paste(start_year,'-01-01',sep="")

#mydata<-Quandl(c('FRED/GDP'),start_date=start_date,end_date=end_date,type='xts')
mydata<-Quandl(c('FRED/GDPC96'),start_date=start_date,end_date=end_date,type='xts')

# Annualized sharpe of GDP series
sharpe_GDP<-sqrt(4)*mean(diff(lgdp))/sqrt(var(diff(lgdp)))
sharpe_GDP

mydata<-na.omit(mydata)
lgdp <- ts(100*log(mydata),start=start_year,frequency=4)
nobs <- length(lgdp)

ss_obj_i1<-ss_model_i1(lgdp)

mat_est_i1<-ss_obj_i1$mat_est
mat_parm_i1<-ss_obj_i1$mat_parm

# All three models are almost identical and cycle estimate is very poor (non-stationary cycle)
plot_ss_i1(mat_parm_i1[3,],lgdp,start_year)


ss_obj_i2<-ss_model_i2(lgdp)

mat_est_i2<-ss_obj_i2$mat_est
mat_parm_i2<-ss_obj_i2$mat_parm

# First and second models are different: level innovation is large. Criterion value gets a big boost.

# Second and fourth models are identical: they have the highest criterion value
# First and third models differ mainly in sigma of drift and persistence of cycle: criterion values almost identical
par(mfrow=c(2,2))
for (i in 1:nrow(mat_parm_i2))
{
  plot_which<-c(F,T)
  title_main<-dimnames(mat_parm_i2)[[1]][i]
  plot_ss_i2(mat_parm_i2[i,],lgdp,start_year,plot_which,title_main)
}


#-----------------------------------------------------------------------------





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
fit2 <- dlmMLE(y=lgdp,parm=c(1.5,-.6,.4,.3),build=ssm2,hessian=T)
mod2 <- ssm2(fit2$par); mod2f <- dlmFilter(lgdp,mod2)





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

# estimate parameters
fit2 <- dlmMLE(y=lgdp,parm=c(1.5,-.6,sqrt(.4),sqrt(.3),sqrt(0.01)),build=ssm2,hessian=T)#,method="SANN")
mod2 <- ssm2(fit2$par); mod2f <- dlmFilter(lgdp,mod2)

(2*pi/Arg(polyroot(c(1,-fit2$par[1:2]))))[1]/4


# get estimates for ARMA(2,0) part
coef <- fit2$par
var <- solve(fit2$hessian)
## print results
coef; sqrt(diag(var))

# get parameter estimates 
drift <- mod2f$m[nobs+1,2]
covar <- dlmSvd2var(mod2f$U.C[[nobs+1]],mod2f$D.C[nobs+1,])
coef.se <- sqrt(covar[2,2])
drift; coef.se

# filtered values
xtfilt <- ts(mod2f$m[-1,1],start=1947,frequency=4)
ctfilt <- ts(mod2f$m[-1,3],start=1947,frequency=4)

# plot fitered state (trend and cycle)
plot(ctfilt,ylim=c(-7,7),xlim=c(1947,2010))
nberShade()
lines(ctfilt)
abline(h=0)
plot(cbind(lgdp,xtfilt),ylim=c(740,960),xlim=c(1947,2010),
     plot.type='s',col=c("black","blue"))
nberShade()
lines(lgdp)
lines(xtfilt)
