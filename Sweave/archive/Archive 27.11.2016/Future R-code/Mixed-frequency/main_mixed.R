
path.main<-paste(disk_id,":\\wia_desktop\\Projekte\\2014\\MDFA-Legacy\\Sweave\\",sep="")
# Source mixed-frequency functions
source(paste(path.main,"Future R-code\\Mixed-frequency\\mixed_frequency_functions.r",sep=""))
# Source regularization matrices: this will source a new version of reg_mat_func 
source(paste(path.main,"Future R-code\\Regularization\\regularization_matrices.r",sep=""))

# Parallel computing
library(doParallel)


#--------------------------------------------------------------
# generate artificial data




# Example 1:
# White noise

len<-1200
set.seed(1)
target_unobserved_diff<-rnorm(len)
x1<-target_unobserved_diff+0.*rnorm(len)
x2<-0.1*rnorm(len)+target_unobserved_diff
x3<-rnorm(len)+target_unobserved_diff
x2[-(5*(1:(len/5)))]<-NA
x3[-(20*(1:(len/20)))]<-NA
target_diff<-target_unobserved_diff
target_diff[-(20*(1:(len/20))+12)]<-NA
  
x_diff<-cbind(target_diff,x1,x2,x3)
xhh_diff<-x_diff
ts.plot(xhh_diff)
head(xhh_diff,100)



#----------------------
# Example 2:
# Random-walks
#   1. The data gets standardized in the estimation file (to ensure consistency with zero-shrinkage/regularization)
#   2. Therefore explaining data should not be contaminated too strongly with noise
#       Otherwise scaling will destroy level information and one will observe seasonal effects in mixed-ferquency estimate
#       because levels of series are different.
#   3. Possible solution: track levels of series

# This setting must be tried with and without regularization!!!
#  with_reg==T and high_freq_only==F (low_freq_only<-F):
#MSE mixed-frequency MSE target MSE target shifted MSE white noise
#[1,]            7.799083   93.27092           109.8302        34.70185
# with_reg==F and high_freq_only==F (low_freq_only<-F):
#MSE mixed-frequency MSE target MSE target shifted MSE white noise
#[1,]            81.15432   93.27092           109.8302        34.70185
# The difference is huge because estimates are entirely out-of-sample!!!!

# with_reg==T and high_freq_only==T (low_freq_only<-F):
#MSE mixed-frequency MSE target MSE target shifted MSE white noise
#[1,]            8.628317   93.27092           109.8302        34.70185
> 
set.seed(3)
target_unobserved_level<-cumsum(rnorm(len))+3*sqrt(20)*rnorm(len)
ts.plot(target_unobserved_level)
target_level<-target_unobserved_level
target_level[-(20*(1:(len/20))+12)]<-NA
# Explaining low-frequency: no noise but delayed
x1<-c(target[1:20],target[1:(length(target)-20)])
# high-frequency explaining and low-frequency target are cointegrated
x2<-target_unobserved_level+sqrt(20)*rnorm(len)
# The lower the frequency the smaller the cointegration term
x3<-x3_unobserved<-sqrt(5)*rnorm(len)+target_unobserved_level 
x4<-x4_unobserved<-rnorm(len)+target_unobserved_level
x3[-(5*(1:(len/5)))]<-NA
x4[-(20*(1:(len/20)))]<-NA
low_freq_only<-F
# Only high-frequency data as explaining
high_freq_only<-T
if (high_freq_only)
{
  x_level<-cbind(target_level,x2)
} else
{
  # Low frequency only
  if (low_freq_only)
  {
    x_level<-cbind(target_level,x1,x3,x4)
  }
}xhh_level<-x_level
ts.plot(xhh_level)

cor(target_unobserved_level,x2)


#----------------------------
# Example 3: The explaining data is strongly contamintade by noise: x2<-target_unobserved_level+10*sqrt(20)*rnorm(len)
#   Interesting positive and negative example depending on regularization settings
# 3.1 Negative/bad example: do not impose shrinkage (with_reg<-F)
#       Effect: seasonality!!!!!!!!!!!!
#       The (high frequency) explaining is strongly contaminated by noise
#       As a consequence  scaling (used for shrinkage!!!) will strongly affect level of (high frequency) series
#       Therefore filter output jumps back to low-frequency series when data is released: seasonal effects can be detected
# 3.2 Positive/good example: impose shrinkage (with_reg<-T)
#       Textbook behaviour: amplitude of low-frequency data shrinks monotonically up to next release time point
#       Weekly is quite strong and monthly less so
#       Shrinkage applies mostly to monthly, then weekly and to a lesser extent to daily.
#       Performance multivariate improves strongly (because high-frequency data is noisy)

# High frequency data only is worse than low-frequency only (high-frequency is very noisy)
# Both are worse than combining everything
# Without regularization out-of-sample performance is awful and seasonality becomes apparent

# with_reg==T and high_freq_only==F (low_freq_only<-F):
#MSE mixed-frequency MSE target MSE target shifted MSE white noise
#[1,]            13.36485   93.27092           109.8302        34.70185

# with_reg==F and high_freq_only==F (low_freq_only<-F):
# MSE mixed-frequency MSE target MSE target shifted MSE white noise
# [1,]            149.8608   93.27092           109.8302        34.70185

# with_reg==T and high_freq_only==T (low_freq_only<-F):
# MSE mixed-frequency MSE target MSE target shifted MSE white noise
# [1,]            26.05531   93.27092           109.8302        34.70185

# with_reg==F and high_freq_only==T (low_freq_only<-F):
# MSE mixed-frequency MSE target MSE target shifted MSE white noise
# [1,]            25.00792   93.27092           109.8302        34.70185

# with_reg==T and low_freq_only<-T (high_freq_only==F):
# MSE mixed-frequency MSE target MSE target shifted MSE white noise
# [1,]            19.03946   93.27092           109.8302        34.70185

# with_reg==F and low_freq_only<-T (high_freq_only==F):
#MSE mixed-frequency MSE target MSE target shifted MSE white noise
#[1,]            185.1382   93.27092           109.8302        34.70185
set.seed(3)
target_unobserved_level<-cumsum(rnorm(len))+3*sqrt(20)*rnorm(len)
ts.plot(target_unobserved_level)
target_level<-target_unobserved_level
target_level[-(20*(1:(len/20))+12)]<-NA
# Explaining low-frequency: no noise but delayed
x1<-c(target[1:20],target[1:(length(target)-20)])
# high-frequency explaining and low-frequency target are cointegrated
x2<-target_unobserved_level+10*sqrt(20)*rnorm(len)
# The lower the frequency the smaller the cointegration term
x3<-x3_unobserved<-sqrt(5)*rnorm(len)+target_unobserved_level
x4<-x4_unobserved<-rnorm(len)+target_unobserved_level
x3[-(5*(1:(len/5)))]<-NA
x4[-(20*(1:(len/20)))]<-NA

x_level<-cbind(target_level,x1,x2,x3,x4)

# Only high-frequency data as explaining
high_freq_only<-F
low_freq_only<-T
if (high_freq_only)
{
  x_level<-cbind(target_level,x2)
} else
{
# Low frequency only
  if (low_freq_only)
  {
    x_level<-cbind(target_level,x1,x3,x4)
  }
}

xhh_level<-x_level
ts.plot(xhh_level)

cor(target_unobserved_level,x2)


#___________________________________________________________________________________------------------------------------------------------------------
#___________________________________________________________________________________------------------------------------------------------------------
# Estimation


# filter settings
with_reg<-T
L<-30
d<-0
lambda<-0
Lag<-0
eta<-0
lin_eta<-F
i1<-F
i2<-F
weight_constraint<-rep(1/(ncol(weight_func)-1),ncol(weight_func)-1)
lambda_cross<-lambda_smooth<-0
lambda_decay<-c(0,0)
if (with_reg)
{
  lambda_cross<-0.
  lambda_smooth<-0.98#0.3
  lambda_decay<-c(0.2,0.98)#c(0.5,0.8)
} else
{
  lambda_cross<-0.
  lambda_smooth<-0
  lambda_decay<-c(0.,0.)  
}
lin_expweight<-F
shift_constraint<-rep(0,ncol(weight_func)-1)
grand_mean<-F
b0_H0<-NULL
c_eta<-F
weights_only<-F
weight_structure<-c(0,0)
white_noise<-F
synchronicity<-F
cutoff<-pi/30
estimation_length<-660
plots<-T

# Re-estimation intervall for DFT: filter is estimated daily but DFT is re-estimated in intervals of length re_estimate_DFT
re_estimate_DFT<-20


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Estimation: non-parallelized algorithm
# Estimation: the NA-structure is changing each day: re-estimate is required each day
# The DFT does not have to be re-estimated each day, though.

# 1. Differenced data

plots<-T
system.time(
est_obj<-estim_mixed_freq(L,estimation_length,len,xhh_diff,d,cutoff,Lag,lambda,eta,i1,i2,weight_constraint,lambda_cross,lambda_decay,lambda_smooth,plots,lin_expweight,shift_constraint,grand_mean,b0_H0,c_eta,weights_only,weight_structure, white_noise,synchronicity,re_estimate_DFT,lin_eta)
)
xf_diff<-est_obj$xf


ord<-200
un_obj<-unobserved_high_freq_signal(target_unobserved_diff,ord,L,estimation_length,len)
filt_target_unobserved_diff<-un_obj$filt_target_unobserved

standardize<-T
plot_obj<-plot_mixed_func(xf_diff,filt_target_unobserved_diff,target_diff,estimation_length,standardize)

plot_obj$perf_mat


# 2. Data in levels

# Estimation: the NA-structure is changing each day: re-estimate is required each day
# The DFT does not have to be re-estimated each day, though.
plots<-T
system.time(
est_obj<-estim_mixed_freq(L,estimation_length,len,xhh_level,d,cutoff,Lag,lambda,eta,i1,i2,weight_constraint,lambda_cross,lambda_decay,lambda_smooth,plots,lin_expweight,shift_constraint,grand_mean,b0_H0,c_eta,weights_only,weight_structure, white_noise,synchronicity,re_estimate_DFT,lin_eta)
)
xf_level<-est_obj$xf


ord<-200
un_obj<-unobserved_high_freq_signal(target_unobserved_level,ord,L,estimation_length,len)
filt_target_unobserved_level<-un_obj$filt_target_unobserved

#ts.plot(cbind(filt_target_unobserved_level,target_unobserved_level))

standardize<-F
plot_obj<-plot_mixed_func(xf_level,filt_target_unobserved_level,target_level,estimation_length,standardize)

plot_obj$perf_mat

acf(na.exclude(xf_level))
acf(na.exclude(diff(xf_level)))

#ts.plot(cbind(xf_level,x4_unobserved)[1100:1200,],lty=1:3)
#mplot<-plot_obj$mplot
#ts.plot(cbind(mplot[,1],target_unobserved_level))




#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
# Parallel computing
# Multi-core setting

# 1. Data in differences
plots<-F
anz_cores<-detectCores() - 1

registerDoParallel(cl=ifelse(plots, 1, anz_cores),cores = ifelse(plots, 1, anz_cores))

system.time(
  est_obj<-estim_mixed_freq_parallel_compute(L,estimation_length,len,xhh_diff,d,cutoff,Lag,lambda,eta,i1,i2,weight_constraint,lambda_cross,lambda_decay,lambda_smooth,plots,lin_expweight,shift_constraint,grand_mean,b0_H0,c_eta,weights_only,weight_structure, white_noise,synchronicity,re_estimate_DFT,lin_eta)
  
)  

xf_diff<-est_obj$xf


ord<-200
un_obj<-unobserved_high_freq_signal(target_unobserved_diff,ord,L,estimation_length,len)
filt_target_unobserved_diff<-un_obj$filt_target_unobserved


standardize<-T
plot_obj<-plot_mixed_func(xf_diff,filt_target_unobserved_diff,target_diff,estimation_length,standardize)

plot_obj$perf_mat

acf(na.exclude(xf_diff))
acf(na.exclude(diff(xf_diff)))


#--------------------------------
# 2. Data in levels

plots<-F
anz_cores<-detectCores() - 1

registerDoParallel(cl=ifelse(plots, 1, anz_cores),cores = ifelse(plots, 1, anz_cores))

system.time(
  est_obj<-estim_mixed_freq_parallel_compute(L,estimation_length,len,xhh_level,d,cutoff,Lag,lambda,eta,i1,i2,weight_constraint,lambda_cross,lambda_decay,lambda_smooth,plots,lin_expweight,shift_constraint,grand_mean,b0_H0,c_eta,weights_only,weight_structure, white_noise,synchronicity,re_estimate_DFT,lin_eta)
    
)  

xf_level<-est_obj$xf


ord<-200
un_obj<-unobserved_high_freq_signal(target_unobserved_level,ord,L,estimation_length,len)
filt_target_unobserved_level<-un_obj$filt_target_unobserved

#ts.plot(cbind(filt_target_unobserved_level,target_unobserved_level))

standardize<-F
plot_obj<-plot_mixed_func(xf_level,filt_target_unobserved_level,target_level,estimation_length,standardize)

plot_obj$perf_mat

acf(na.exclude(xf_level))
acf(na.exclude(diff(xf_level)))


