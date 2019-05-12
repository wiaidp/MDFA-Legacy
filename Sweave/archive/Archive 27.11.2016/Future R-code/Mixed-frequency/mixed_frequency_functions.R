


# Function for determining lag_mat

#x<-isData
compute_lags_mixed_frequency<-function(x,L,ind)
{
  lag_mat<-matrix(ncol=ncol(x)-1,nrow=L)
  for (j in 2:ncol(x))#j<-3
  {
    z<-x[,j]
    lag_mat[,j-1]<-which(!is.na(z[ind:1]))[1:L]-1
    
  }
  return(list(lag_mat=lag_mat))
}


# x<-isData
# insamp<-nrow(x)
spec_comp <- function(insamp, x, d) 
{

#target    
  weight_func <- general_spec(x[1 : insamp, 1], 0, insamp)$fourtrans
# explaining variables
  if(length(weight_func) > 1) {
    for(j in 2 : ncol(x)) #j<-2
    {
      DFT <- general_spec(x[1 : insamp, j], 0, insamp)$fourtrans
      weight_func <- cbind(weight_func, DFT)
    }
  }
  colnames(weight_func) <- colnames(x)

# Re-scale the amplitude functions of the low sampling frequency:
#   Sum of spectrum must be identical with variance
#   Folding of the DFT changes the scale!!!!
  for (i in 1:ncol(weight_func))#i<-2
  {
    sigma<-var(x[,i],na.rm=T)
#    print(sigma)
    weight_func[,i]<-weight_func[,i]*sigma/sqrt(2*pi*mean(abs(weight_func[,i])^2))
#    print(2*pi*mean(abs(weight_func[,i])^2))
  }    
# return the spectral estimate
  return(list(weight_func = weight_func))
}



#  ts.plot(((Im(weight_func))[1:30,c(1,2)]),lty=1:4)   ts.plot(((Re(weight_func))[1:30,c(1,2)]),lty=1:4)
#  ts.plot(((abs(weight_func))[1:30,c(1,3)]),lty=1:4)






# This function computes the DFT in the case of general mixed_frequency data
#   1. It accounts for NA's
#   2. It folds the DFT of the low-frequency data

#y<-x[,1]
general_spec<-function(y, d, insamp)
{
# generate even length (if length is odd then oldest observation is removed)
  if (!(as.integer(length(y)/2)==length(y)/2))
    y<-y[-1]
# Remove NA's (if present)
  z<-y[!is.na(y)]
  # we could compute the lag until the first observation is available and rotate the DFT accordingly
  # But a better/cleaner solution is to postpone these relative lags of the data to the computations 
  # of the MDFA spectral matrix
#  which(y==y[!is.na(y)][1])[1]-1
  # compute the DFT of z: undo scaling by the number of observations!!!!  
#  dft<-periodogram_bp(z,d, length(z))$fourtrans#*sqrt(pi*length(z))   sum(x[,1],na.rm=T)/sqrt(sum(!is.na(x[,1])))  sum(x[,2],na.rm=T)/sqrt(sum(!is.na(x[,2])))

# Use standard per function in I-DFA folder of MDFA-Legacy project
  dft<-per(z-mean(z),F)$DFT

# Extending the high-frequency part
#   Assumption: AR(1)  

  if (length(z)<length(y))
  {
    ar_obj<-arima(z-mean(z),order=c(1,0,0),include.mean=F)
    coef<-ar_obj$coef
    sigmas<-ar_obj$sigma
# align value for pi (on low-frequency scale) to last frequency on hihj-ferquency scale
    value_pi<-sqrt(sigmas)/(sqrt(2*pi)*(1-coef*exp(1.i*(pi/(length(dft)-1))*((length(dft)-1)))))
# Extension of frequencies    
    freq_ext<-(length(dft)):(as.integer(length(y)/2))
# Extension of DFT based on AR(1)-model
    dft_ext<-sqrt(sigmas)/(sqrt(2*pi)*(1-coef*exp(-1.i*(pi/(length(y)/2))*freq_ext)))
# Append folded DFT and DFT-extension
    fourtrans<-c(dft*abs(dft_ext[1]/value_pi),dft_ext)   #ts.plot(Re(fourtrans))  #ts.plot(Im(fourtrans))  #ts.plot(abs(fourtrans))
  } else
  {
    fourtrans<-dft
  }
  return(list(fourtrans=fourtrans))
}




#set.seed(1)
#dft<-per((rnorm(100)),F)$DFT
#ts.plot(Re(dft))
#ts.plot(Im(dft))
#ts.plot(Arg(dft))
#ts.plot(Arg(dft)/(pi*(0:(length(dft)-1))/(length(dft)-1)))


IMDFA_comp <- function(Lag, K, L, lambda, weight_func, Gamma, eta, cutoff,
                       i1, i2, weight_constraint, lambda_cross, lambda_decay,
                       lambda_smooth, x, plots,  shift_constraint,
                       grand_mean, b0_H0, c_eta, weights_only,weight_structure,
                       white_noise,synchronicity,lag_mat,lin_eta) {
  # call the mdfa
  #  i_mdfa <- mdfa_analytic(K, L, lambda, weight_func, Lag, Gamma, expweight,
  #                              cutoff, i1, i2, weight_constraint, lambda_cross,
  #                              lambda_decay, lambda_smooth,  
  #                              shift_constraint, grand_mean, b0_H0,
  #                              chris_expweight, weights_only,weight_structure,white_noise,
  #                              synchronicity,lag_mat)
  i_mdfa<-mdfa_analytic(K,L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                        lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                        b0_H0,c_eta,weights_only=F,weight_structure,white_noise,
                        synchronicity,lag_mat)
  
  # print the number of freezed degrees and the criterion value
  print(paste("Freezed degrees =", i_mdfa$freezed_degrees))
  print(paste("Criterion value =", i_mdfa$rever))
  # plot diagnostic charts
  if(plots) {
    par(mfrow = c(3, 1))
    # amplitude functions
    mplot <- abs(i_mdfa$trffkt)
    # x-axis
    freq_axe <- rep(NA, K + 1)
    freq_axe[1] <- 0
    freq_axe[1 + (1 : 6) * K / 6] <- c(paste0(c("", 2 : 5), "pi/6"), "pi")
    ax <- freq_axe
    # colors, title and additional titles
    insamp <- 1.e+90
    colo <- NULL
    plot_title <- "Amplitude Functions"
    title_more <- colnames(x[, -1])
    mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
    
    # time-shift
    mplot <- Arg(t(sign(apply(i_mdfa$b, 2, sum)) * t(i_mdfa$trffkt))) /
      ((0 : (nrow(i_mdfa$trffkt) - 1)) * pi / (nrow(i_mdfa$trffkt) - 1))
    mplot[1, ] <- apply(i_mdfa$b * ((0 : (L - 1))), 2, sum) / 
      apply(i_mdfa$b, 2, sum)
    plot_title <- "Time-Shift"
    mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
    
    # filter coefficients
    mplot <- i_mdfa$b
    ax <- Lag + 0 : (L-1)
    plot_title <- "Filter Coefficients"
    mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
  }
  # first order constraint
  print(paste("First order restriction:", apply(i_mdfa$b, 2, sum)))
  # 02.08.2013: correct formula for time-shift: a warning is issued when the 
  # shift depends non-linearly on filter coefficients
  # second order constraint
  if(sum(abs(shift_constraint)) > 0 & i2 & !i1) {
    print(rep("!", 100))
    print("Warning: shift constraint is different from zero and i2 <- TRUE,
          i1 == FALSE. Then the problem cannot be solved analytically, see the 
          elements-paper, section 8.1.2. You may have to try some iterative 
          search procedure instead")
    print(rep("!", 100))
  }
  if(i1 & any(weight_constraint == 0)) {
    print(rep("!", 100))
    print("Warning: filter vanishes in frequency zero and therefore the 
          time-shift is not defined")
    print(rep("!", 100))
  }
  print(paste("Second order constraint:", 
              apply(i_mdfa$b * (0 : (L - 1)), 2, sum) / apply(i_mdfa$b, 2, sum)))
  # 02.08.2013: end correction formula for time-shift
  
  # filter the series
  xff <- NULL
  xff <- matrix(nrow = nrow(as.matrix(x[, -1])), 
                ncol = ncol(as.matrix(x[, -1])))
  # bn <- scale(i_mdfa$b, center = TRUE, scale = FALSE)
  bn <- as.matrix(i_mdfa$b)
  for(i in L : nrow(x)) {
    xff[i, ] <- 0
    for(j in 2 : ncol(x)) {
      xff[i, j - 1] <- xff[i, j - 1] + bn[, j - 1] %*% x[i : (i - L + 1), j]
    }
  }
  return(list(xff = xff, i_mdfa = i_mdfa))
}


#xhh<-xhh_level
#xhh<-xhh_diff
#lambda_cross<-0.9
#lambda_smooth<-0.9
#lambda_decay<-c(0.1,0.8)
#plots<-T
# cor((Re(weight_func))[1:16,c(1,3)])    cor((Im(weight_func))[1:16,c(1,3)])  cor((abs(weight_func))[1:16,c(1,3)])

estim_mixed_freq<-function(L,estimation_length,len,xhh,d,cutoff,Lag,lambda,eta,i1,i2,weight_constraint,
                           lambda_cross,lambda_decay,lambda_smooth,plots,shift_constraint,
                           grand_mean,b0_H0,c_eta,weights_only,weight_structure, white_noise,synchronicity,re_estimate_DFT,lin_eta)
{
  xf_series<-as.matrix(xhh[,-1])#head(xhh)
  xf<-rep(NA,nrow(xhh))
  for (i in max(L,estimation_length):len) #i<-estimation_length #i<-1196
  {
    print(i)
    isData<-xhh[(i-estimation_length+1):i,]  #var(isData[,2],na.rm=T)
# Data is scaled: to be in accordance with zero-shrinkage!!!
    for (j in 1:ncol(isData))
      isData[,j]<-isData[,j]/sqrt(var(xhh[(i-estimation_length+1):(i-1),j],na.rm=T))
# Re-estimate DFT once per month: MDFA will be full out-of-sample because 
#   1. Reestimation is only once in an interval of length re_estimate_DFT interval (avoid unnecessary revisions)
#   2. Once re-estimated, last observation of isData is removed
    if (i==max(L,estimation_length)|as.integer(i/re_estimate_DFT)==i/re_estimate_DFT)
      weight_func <- spec_comp(estimation_length-1, isData[-nrow(isData),], d)$weight_func#head(weight_func,100)
      
    #  ts.plot(((Re(weight_func))[1:30,c(1,3)]),lty=1:4)  2*pi*apply(abs(weight_func)^2,2,mean)
    #  ts.plot(((Im(weight_func))[1:30,c(1,3)]),lty=1:4)  2*pi*apply(abs(weight_func)^2,2,mean)
    #  ts.plot(((abs(weight_func))[,c(1,3)]),lty=1:4)  #ts.plot(xhh[,3])
    # compute the ideal symmetric filter
    K <- nrow(weight_func) - 1
    Gamma <- (0 : K) < (cutoff * K / pi)
    
    # Compute the effective lags: here the daily lag-structure is relevant
    
    lag_mat<-compute_lags_mixed_frequency(isData,L,nrow(isData))$lag_mat
    
    # Estimate filter coefficients:
    #  1. They are based on a daily lag-structure which changes each day in case of mixed-frequencies (or NA's) 
    #  2. The DFT is kept fixed over a time span determined by re_estimate_dft_Index 
    #  3. The DFTs are not refreshed on a daily basis in order to avoid revision noise
    
    res <- IMDFA_comp(Lag, K, L, lambda, weight_func, Gamma, eta, cutoff,i1, i2, weight_constraint, lambda_cross, lambda_decay, lambda_smooth, isData, plots,  shift_constraint, grand_mean, b0_H0, c_eta, weights_only,weight_structure, white_noise,synchronicity,lag_mat,lin_eta) 
    
    # extract the coefficients
    bn <- res$i_mdfa$b
    # filter the data
    for (j in 1:ncol(bn))#j<-3
    {
      filter_data<-as.matrix(isData[,-1])#ts.plot(isData)
# Remove NA's from low-frequency data    
      exp_data<-na.exclude(filter_data[,j])
      xf_series[i,j]<-bn[,j] %*% exp_data[length(exp_data):(length(exp_data)-L+1)]#xf_series[i,]
    }
    xf[i]<-sum(xf_series[i,])*sqrt(var(xhh[(i-estimation_length+1):(i-1),1],na.rm=T))   
  }
# Compute aggregate series and re-scale back to original scale of target
  return(list(xf=xf))
}








estim_mixed_freq_parallel_compute<-function(L,estimation_length,len,xhh,d,cutoff,Lag,lambda,eta,i1,i2,weight_constraint,
                           lambda_cross,lambda_decay,lambda_smooth,plots,shift_constraint,
                           grand_mean,b0_H0,c_eta,weights_only,weight_structure, white_noise,synchronicity,re_estimate_DFT,lin_eta)
{
  
  obj_par <- foreach(i = max(L,estimation_length):nrow(xhh), .combine = rbind, .verbose = TRUE,.export = 
c("estim_mixed_freq_parallel_comp","spec_comp","per","IMDFA_comp","mdfa_analytic","spec_mat_comp","mat_func","mplot_func",
  "general_spec","compute_lags_mixed_frequency","centraldev_original_func","reg_mat_func","w_eight_func","des_mat_func","structure_func","white_noise_synchronicity","Q_reg_func")) %dopar%
  estim_mixed_freq_parallel_comp(i,L,estimation_length,nrow(xhh),xhh,d,cutoff,Lag,lambda,eta,i1,i2,weight_constraint,lambda_cross,lambda_decay,lambda_smooth,plots,shift_constraint,grand_mean,b0_H0,c_eta,weights_only,weight_structure, white_noise,synchronicity,re_estimate_DFT,lin_eta)
# assign the coefficients and the in- and out-of-sample filtered series to the  
# local environment
  xf<-c(rep(NA,(max(L,estimation_length)-1)),unlist(obj_par))  
  return(list(xf=xf))
}  
















estim_mixed_freq_parallel_comp<-function(i,L,estimation_length,len,xhh,d,cutoff,Lag,lambda,eta,i1,i2,weight_constraint,
                           lambda_cross,lambda_decay,lambda_smooth,plots,shift_constraint,
                           grand_mean,b0_H0,c_eta,weights_only,weight_structure, white_noise,synchronicity,re_estimate_DFT,lin_eta)
{

  isData<-xhh[(i-estimation_length+1):i,] #i<-1001
# DFT is not re-estimated each day but in intervals of length re_estimate_DFT
# In parallel computing the DFT must be re-estimated for each i because the previous DFT is not available (as is the case when not running code in parallel)
# MDFA is full out-of-sample because 
#   1. Reestimation is only once in an interval of length re_estimate_DFT interval (avoid unnecessary revisions)
#   2. Once re-estimated, last observation is removed i.e. the last observation in isData (for which MDFA is run) will always be out-of-sample
  DFT_index<-as.integer(i/re_estimate_DFT)*re_estimate_DFT-1
  DFT_data<-xhh[max(1,DFT_index-estimation_length+2):DFT_index,]
# Data is scaled: to be in accordance with zero-shrinkage!!!
  for (j in 1:ncol(isData))
  {
    isData[,j]<-isData[,j]/sqrt(var(xhh[(i-estimation_length+1):(i-1),j],na.rm=T))
    DFT_data[,j]<-DFT_data[,j]/sqrt(var(xhh[(DFT_index-estimation_length+2):DFT_index,j],na.rm=T))
  }
# Re-estimate DFT for each i (must be done when running code in parallel because DFT is not shared amongts threads)
  weight_func <- spec_comp(nrow(DFT_data), DFT_data, d)$weight_func#head(weight_func,100)

    #  ts.plot(((Re(weight_func))[1:30,c(1,3)]),lty=1:4)  2*pi*apply(abs(weight_func)^2,2,mean)
    #  ts.plot(((Im(weight_func))[1:30,c(1,3)]),lty=1:4)  2*pi*apply(abs(weight_func)^2,2,mean)
    #  ts.plot(((abs(weight_func))[,c(1,3)]),lty=1:4)  #ts.plot(xhh[,3])
    # compute the ideal symmetric filter
  K <- nrow(weight_func) - 1
  Gamma <- (0 : K) < (cutoff * K / pi)
    
    # Compute the effective lags: here the daily lag-structure is relevant
    
  lag_mat<-compute_lags_mixed_frequency(isData,L,nrow(isData))$lag_mat
    
    # Estimate filter coefficients:
    #  1. They are based on a daily lag-structure which changes each day in case of mixed-frequencies (or NA's) 
    #  2. The DFT is kept fixed over a time span determined by re_estimate_dft_Index 
    #  3. The DFTs are not refreshed on a daily basis in order to avoid revision noise
    
  res <- IMDFA_comp(Lag, K, L, lambda, weight_func, Gamma, eta, cutoff,i1, i2, weight_constraint, lambda_cross, lambda_decay, lambda_smooth, isData, plots,  shift_constraint, grand_mean, b0_H0, c_eta, weights_only,weight_structure, white_noise,synchronicity,lag_mat,lin_eta) 
    
    # extract the coefficients
  bn <- res$i_mdfa$b
# filter the data
  xf_series<-rep(NA,ncol(bn))
  for (j in 1:ncol(bn))#j<-3
  {
    filter_data<-as.matrix(isData[,-1])#ts.plot(isData)
# Remove NA's from low-frequency data    
    exp_data<-na.exclude(filter_data[,j])
    xf_series[j]<-bn[,j] %*% exp_data[length(exp_data):(length(exp_data)-L+1)]#xf_series[i,]
  }

  # Compute aggregate series and re-scale back to original scale of target
  xf<-sum(xf_series)*sqrt(var(xhh[(i-estimation_length+1):(i-1),1],na.rm=T))
  
  return(list(xf=xf))
}









# Makes a step-function from a time series with NA's (random-walk extrapolation)
na.complete_func<-function(x)
{
  x<-as.matrix(x)
  x_step<-x
  for (i in 1:ncol(x))#i<-1
  {
    na_ex<-na.exclude(x[,i])
    na_ind<-which(!is.na(x[,i]))
    if (na_ind[1]>1)
      x_step[1:na_ind[1],i]<-na_ex[1]
    for (j in 1:(length(na_ex)-1))
    {
      
      x_step[na_ind[j]:(na_ind[j+1]-1),i]<-na_ex[j]
    }
    if (na_ind[length(na_ind)]<nrow(x))
      x_step[na_ind[length(na_ind)]:nrow(x),i]<-na_ex[length(na_ex)]
    
  }
  return(list(x_step=x_step))
}





# Computes symmetric filter output
unobserved_high_freq_signal<-function(target_unobserved,ord,L,estimation_length,len)
{
  # Symmetric filter applied to unobserved high-frequency target
  gamma<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))
  sum(gamma)+sum(gamma[2:ord])
  filt_target_unobserved<-rep(NA,length(target_unobserved))
  for (i in max(L,estimation_length):(len-ord)) #i<-estimation_length
  {
    filt_target_unobserved[i]<-gamma%*%target_unobserved[i:(i-ord)]+gamma[2:ord]%*%target_unobserved[(i+1):(i+ord-1)]
  }
  return(list(filt_target_unobserved=filt_target_unobserved))
}







#xf<-xf_level
#filt_target_unobserved<-filt_target_unobserved_level
#target<-target_level
plot_mixed_func<-function(xf,filt_target_unobserved,target,estimation_length,standardize)
{
  
  target_step<-target
  na_vec<-which(!is.na(target))
  for (i in 1:(length(na_vec)-1))
  {
    target_step[(na_vec[i]+1):na_vec[i+1]]<-target[na_vec[i]]
  }
  #  ts.plot(target_step)
  if (standardize)
  {
    mplot<-scale(cbind(target_step,filt_target_unobserved,xf,target)[estimation_length:length(xf),])
  } else
  {
    mplot<-cbind(target_step,filt_target_unobserved,xf,target)[estimation_length:length(xf),]
  }
  dimnames(mplot)[[2]]<-c("Monthly target","Ideal daily target","Real-time estimate","Target")
  perf_mat<-matrix(c(mean(abs(mplot[,3]-mplot[,2])^2,na.rm=T),mean(abs(mplot[,1]-mplot[,2])^2,na.rm=T),mean(abs(c(mplot[1:20,1],mplot[1:(nrow(mplot)-20),1])-mplot[,2])^2,na.rm=T),mean(abs(mplot[,2])^2,na.rm=T)),nrow=1)
  dimnames(perf_mat)[[2]]<-c("MSE mixed-frequency","MSE target","MSE target shifted","MSE white noise")
# remove low-frequency target
  mplot<-mplot[,-ncol(mplot)]
# Re-scale on commaon available time-span (symmetric output is not available towards sample end)
  if (standardize)
  {
    mplot<-scale(na.exclude(mplot))
  } 
  ax<-rep(NA,nrow(mplot))
  ax[c(1,rep(0,6))+(0:6)*((nrow(mplot)-1)/6)]<-as.integer(c(1,rep(0,6))+(0:6)*((nrow(mplot)-1)/6))
  plot_title<-"Mixed-frequency: tracking a monthly target"
  insamp<-1.e+90
  title_more<-dimnames(mplot)[[2]]
  colo<-rainbow(dim(mplot)[2])
  mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
  return(list(mplot=mplot,perf_mat=perf_mat))  
}






