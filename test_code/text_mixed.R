leni<-100

estim_mixed_freq<-function(L,estimation_length,len,xhh,d,cutoff,Lag,lambda,eta,i1,i2,weight_constraint,
                           lambda_cross,lambda_decay,lambda_smooth,plots,lin_expweight,shift_constraint,
                           grand_mean,b0_H0,c_eta,weights_only,weight_structure, white_noise,synchronicity)
{
  xf_series<-xhh[,-1]#head(xhh)
  for (i in (len-leni):len) #i<-estimation_length #i<-1196
  {
    print(i)
    isData<-xhh[(i-estimation_length+1):i,]  #var(isData[,2],na.rm=T)
    # Data is scaled: to be in accordance with zero-shrinkage!!!
    for (j in 1:ncol(isData))
      isData[,j]<-isData[,j]/sqrt(var(xhh[(i-estimation_length+1):i,j],na.rm=T))
    #    mean(isData[,4],na.rm=T)
    # Re-estimate DFT once per month  
    if (i==max(L,estimation_length)|as.integer(i/20)==i/20)
      weight_func <- spec_comp(estimation_length, isData, d)$weight_func#head(weight_func,100)
    
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
    
    res <- IMDFA_comp(Lag, K, L, lambda, weight_func, Gamma, eta, cutoff,i1, i2, weight_constraint, lambda_cross, lambda_decay, lambda_smooth, isData, plots, lin_expweight, shift_constraint, grand_mean, b0_H0, c_eta, weights_only,weight_structure, white_noise,synchronicity,lag_mat) 
    
    # extract the coefficients
    bn <- res$i_mdfa$b
    # filter the data
    for (j in 1:ncol(bn))#j<-3
    {
      filter_data<-isData[,-1]#ts.plot(isData)
      # Remove NA's from low-frequency data    
      exp_data<-na.exclude(filter_data[,j])
      xf_series[i,j]<-bn[,j] %*% exp_data[length(exp_data):(length(exp_data)-L+1)]#xf_series[i,]
    }
  }
  # Compute aggregate series and re-scale back to original scale of target
  xf<-apply(xf_series,1,sum)*sqrt(var(xhh[(i-estimation_length+1):i,1],na.rm=T))
  xf_h<-apply(xf_series,1,sum)
  
  return(list(xf=xf))
}


ts.plot(diff(xf[(len-leni):len]))
ts.plot(xf_h[(len-leni):len])

ts.plot(xf[(len-20):len])

par(mfrow=c(3,1))
ts.plot(diff(xf[(len-20):len]))
ts.plot(apply(xf_series[(len-20):len,],2,diff),lty=1:ncol(xf_series))
ts.plot(step_isData[(nrow(isData)-20):nrow(isData),-1],lty=1:ncol(isData[,-1]))


step_isData<-na.complete_func(isData)$x_step

