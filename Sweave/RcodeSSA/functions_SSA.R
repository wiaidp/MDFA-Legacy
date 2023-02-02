# Data  function: loads data for BCA in paper

data_load_func<-function(path.data) 
{  
  
  indpro_mat<-read.csv(paste(path.data,"/indpro.csv",sep=""),sep=",",header=T,na.strings="NA",dec=".",row.names=1)
  
  indpro<-indpro_level<-NULL
  for (i in 1:ncol(indpro_mat))
  {
    indpro<-cbind(indpro,diff(log(indpro_mat[,i])))
    indpro_level<-cbind(indpro_level,indpro_mat[,i])
  }
  typeof(indpro)
  colnames(indpro)<-colnames(indpro_level)<-colnames(indpro_mat)
  rownames(indpro)<-rownames(indpro_mat)[2:nrow(indpro_mat)]
  rownames(indpro_level)<-rownames(indpro_mat)
  mean(indpro_mat[,1],na.rm=T)
  
  indpro_mat_eu<-read.csv(paste(path.data,"/indpro_eu_sa.csv",sep=""),sep=",",header=T,na.strings="NA",dec=".",row.names=1)
  
  indpro_eu<-NULL
  for (i in 1:ncol(indpro_mat_eu))
    indpro_eu<-cbind(indpro_eu,diff(log(as.double(indpro_mat_eu[,i]))))
  typeof(indpro_eu)
  tail(indpro_eu)
  colnames(indpro_eu)<-colnames(indpro_mat_eu)
  rownames(indpro_eu)<-rownames(indpro_mat_eu)[2:nrow(indpro_mat_eu)]
  
  indpro<-as.xts(indpro,order.by=as.Date(rownames(indpro),"%d/%m/%Y"))
  indpro_level<-as.xts(indpro_level,order.by=as.Date(rownames(indpro_level),"%d/%m/%Y"))
  return(list(indpro=indpro,indpro_level=indpro_level,indpro_eu=indpro_eu))
} 

#-------------------------------------------------------------------
# HP-functions

# Computes target, MSE, HP-trend and HP-gap original and modified: it is used in BCA-section of paper
HP_target_mse_modified_gap<-function(L,lambda_monthly)
{
  #   MSE relies on white noise assumption while HP-concurrent relies on implicit ARIMA(0,2,2) model
  #   L<-100 is OK i.e. recession datings are nearly identical with L<-200
  setseed<-1
  
  hp_obj<-hp_func(L,lambda_monthly,setseed)
  
  # Concurrent trend
  hp_trend<-hp_obj$concurrent
  ts.plot(hp_trend)
  # Concurrent gap
  hp_gap<-c(1-hp_trend[1],-hp_trend[2:L])
  ts.plot(hp_gap)
  # Modified concurrent gap (as applied to first differences)
  modified_hp_gap<-hp_gap
  for (i in 1:length(hp_gap))
  {
    modified_hp_gap[i]<-sum(hp_gap[1:i])
  }
  ts.plot(modified_hp_gap)
  # Symmetric target
  target<-hp_obj$target
  ts.plot(target)
  # One-sided MSE: must double length in order to retrieve right half of target
  L_target<-2*(L-1)+1
  hp_obj<-hp_func(L_target,lambda_monthly,setseed)
  target_long<-hp_obj$target
  hp_mse<-target_long[(1+(L_target-1)/2):L_target]
  ts.plot(hp_mse)
  return(list(hp_mse=hp_mse,hp_gap=hp_gap,modified_hp_gap=modified_hp_gap,hp_trend=hp_trend,target=target))
}


# Generic HP function relying on R-package. Computes holding-times according to formula in paper
hp_func<-function(L,lambda,setseed)
{
  
  set.seed(setseed)
  eps<-rnorm(L)
  
  hp_filt_obj<-hp_filt_obj <- hpfilter(eps,type="lambda", freq=lambda)
  
  gap_matrix<-hp_filt_obj$fmatrix
  # Extract the coefficients of the symmetric trend:
  #   hpfilter generates coefficients of the HP-gap (see below):
  #   we here transform back to trend filter
  parm_hp<-(diag(rep(1,L))-hp_filt_obj$fmatrix)
  target<-parm_hp[,(L-1)/2+1]
  rho_ht_hp<-compute_holding_time_func(target)
  ht_target<-rho_ht_hp$ht
  concurrent<-parm_hp[,1]
  ht_concurrent<-compute_holding_time_func(concurrent)$ht
  b_mse<-target[((length(target)-1)/2+1):length(target)]
  if (length(b_mse)>L)
    b_mse<-b_mse[1:L]
  if (length(b_mse)<L)
    b_mse<-c(b_mse,rep(0,L-length(b_mse)))
  #  ts.plot(b_mse)
  return(list(target=target,ht_target=ht_target,concurrent=concurrent,ht_concurrent=ht_concurrent,b_mse=b_mse,gap_matrix=gap_matrix))
}

#----------------------------------------------------
# The following three functions compute holding-times and link lag-one acf and holding-times, see section 1 in paper
compute_holding_time_func<-function(b)
{
  rho_ff1<-b[1:(length(b)-1)]%*%b[2:length(b)]/sum(b^2)
  # Mean holding-time
  ht<-1/(2*(0.25-asin(rho_ff1)/(2*pi)))
  # Alternative expression  
  if (F)
    ht<-pi/acos(rho_ff1)
  
  return(list(ht=ht,rho_ff1=rho_ff1))
}

compute_rho_from_ht<-function(ht)
{  
  rho<-sin(-((1/ht)/2-0.25)*2*pi)
  return(list(rho=rho))
}

compute_holding_time_from_rho_func<-function(rho_ff1)
{
  # Mean holding-time
  ht<-1/(2*(0.25-asin(rho_ff1)/(2*pi)))
  return(list(ht=ht))
}

#-------------------------------------------------------
# SSA functions
# This function is used in section 4 of paper
# It computes HP-filters (target, concurrent) and proceeds to SSA-estimation , see SSA_func below
# The assumption for SSA is white noise data: an extension to MA(1)-processes (untested) is provided and this code could be readily extended to ARMA-processes
SSA_compute<-function(ht,L,hp_mse,forecast_horizon_vec,MA1_adjustment,grid_size)
{  
  # SSA hyperparameters
  # Holding-time transformed as lag-one acf
  rho0<-compute_rho_from_ht(ht)$rho
  # Size of grid for determining nu: larger means better but time-consuming
  # Include negative lambda1: default is F (lowpass does not require negative values)
  with_negative_lambda<-F
  # Target gammak_generic is MSE-filter: SSA must be as close to MSE as possible while complying with holding-time   
  gammak_generic<-hp_mse
  
  
  #---------------------------------------------------
  # SSA filter, Assumption: data is white noise
  
  # 3.1 Compute SSA filters
  # Takes a couple seconds
  SSA_obj<-SSA_func(L,forecast_horizon_vec,grid_size,gammak_generic,rho0,with_negative_lambda)
  
  bk_mat=SSA_obj$bk_mat
  colnames(bk_mat)<-paste("SSA(",round(ht,2),",",forecast_horizon_vec,")",sep="")
  
  colo<-c("brown",rainbow(2*ncol(bk_mat)-1)[(ncol(bk_mat)+1):(2*ncol(bk_mat))])
  ts.plot(scale(bk_mat,center=F,scale=rep(T,ncol(bk_mat))),col=colo)
  for (i in 1:ncol(bk_mat))
    mtext(colnames(bk_mat)[i],col=colo[i],line=-i)
  
  # This is not in use in paper: should be checked
  MA1_adjustment<-F
  if (MA1_adjustment)
  {  
    # 4. SSA filter, Assumption: data is MA(1)  
    # 4.1 Fit MA(1) to log-returns
    ma1<-arima(series,order=c(0,0,1))$coef["ma1"]
    
    ma1
    
    # 4.2 Specify new target as applied to white noise
    # Transformed target for white noise (accounts for convolution of original target and MA(1) data-generating filter)
    gammak_generic_ma1<-c(gammak_generic[1],gammak_generic[1:(length(gammak_generic)-1)]*ma1+gammak_generic[2:(length(gammak_generic))])
    
    ts.plot(cbind(gammak_generic,gammak_generic_ma1),col=c("black","blue"))
    
    # 4.3 Fit SSA to new target
    SSA_obj_ma1<-SSA_func(L,forecast_horizon_vec,grid_size,gammak_generic_ma1,rho0,with_negative_lambda)
    
    bk_mat_ma1=SSA_obj_ma1$bk_mat
    colnames(bk_mat_ma1)<-paste("SSA(",ht,"): delta=",forecast_horizon_vec,sep="")
    
    
    # SSA filters: as applied to white noise
    filter_mat_ma1<-cbind(hp_trend,gammak_generic_ma1,bk_mat_ma1)
    colnames(filter_mat_ma1)<-c("HP trend","MSE",paste("SSA(",ht,"): delta=",forecast_horizon_vec,sep=""))
    
    # 4.4 Define filters as applied to data i.e. MA(1)-process (invert convolution of bk and MA(1))
    filter_mat_data<-filter_mat_ma1
    for (i in 1:ncol(filter_mat))
      filter_mat_data[2:nrow(filter_mat_data),i]<-filter_mat_ma1[2:nrow(filter_mat_ma1),i]-ma1*filter_mat_ma1[-1+2:nrow(filter_mat_data),i]
    
    ts.plot(filter_mat_data)
    
    # 4.5 Filter data with transformed filters
    
    filter_obj_ma1<-SSA_filter_func(filter_mat_data,L,series)
    
    y_mat_ma1=filter_obj_ma1$y_mat
    
    colo<-c("brown",rainbow(2*ncol(bk_mat)-1)[(ncol(bk_mat)+1):(2*ncol(bk_mat))])
    selection_vec<-c(1,2,6)
    anf<-1
    anf<-nrow(y_mat_ma1)-400
    enf<-nrow(y_mat_ma1)
    ts.plot(scale(y_mat_ma1[anf:enf,selection_vec],center=F,rep(T,ncol(y_mat_ma1))),col=colo[selection_vec])
    for (i in 1:length(selection_vec))
      mtext(colnames(y_mat_ma1)[selection_vec[i]],col=colo[selection_vec[i]],line=-i)
    abline(h=0)
    
    
    
    # number of crossings
    number_cross<-rep(NA,ncol(filter_mat_ma1))
    names(number_cross)<-colnames(filter_mat_ma1)
    for (i in 1:ncol(y_mat_ma1))
    {
      if (is.xts(y_mat))
      {  
        number_cross[i]<-length(which(sign(y_mat_ma1[,i])!=sign(lag(y_mat_ma1[,i]))))
      } else
      {
        number_cross[i]<-length(which(sign(y_mat_ma1[1:(nrow(y_mat_ma1)-1),i])!=sign(lag(y_mat_ma1[2:nrow(y_mat_ma1),i]))))
      }
    }
    
    number_cross
    # empirical holding time: larger than ht
    nrow(y_mat_ma1)/number_cross
  }  
  
  
  
  return(list(bk_mat=bk_mat))
}



# This function computes exact SSA-filter based on corollary in paper: it relies on find_lambda1_subject_to_holding_time_constraint_func 
#   -grid-search of lambda for given target gammak_generic (target=MSE estimate i.e. not symmetric filter) 
#     rho0 and length L in grid with resolution grid_size
#   -For smoothing one needs positive values of lambda only i.e. with_negative_lambda==T
# It computes optimal estimates for various forecast horizons as specified in forecast_horizon_vec
# forecast_horizon_vec<-delta
SSA_func<-function(L,forecast_horizon_vec,grid_size,gammak_generic,rho0,with_negative_lambda)
{  
  bk_mat<-NULL
  # Loop over all forecast horizons  
  for (i in 1:length(forecast_horizon_vec))#i<-1
  {  
    print(paste(round(100*i/(1+length(forecast_horizon_vec)),0),"%",sep=""))
    forecast_horizon<-forecast_horizon_vec[i]
    opt_obj<-find_lambda1_subject_to_holding_time_constraint_func(grid_size,L,gammak_generic,rho0,forecast_horizon,with_negative_lambda)
    bk_mat<-cbind(bk_mat,opt_obj$bk_best)
    lambda_opt<-opt_obj$lambda_opt
  }
  
# Check: for a lowpass design preserving signs the sum of the filter coefficients should be strictly positive
  apply(bk_mat,2,sum)
  
  
# Check theoretical holding-times of optimum: should match ht-constraint
  apply(bk_mat,2,compute_holding_time_func)
  
  return(list(bk_mat=bk_mat,lambda_opt=lambda_opt))
} 






# This function finds optimal lambda (note that nu=lambda+1/lambda)) through grid search based on function compute_bk_from_ma_expansion_with_boundary_constraint_func above
# It implements solution in corollary 1 of paper
find_lambda1_subject_to_holding_time_constraint_func<-function(grid_size,L,gammak_generic,rho0,forecast_horizon,with_negative_lambda)
{
  Lambda<-(1:grid_size)/(grid_size+1)
  if (with_negative_lambda)
  { 
    # Allow negative lambda1 too    
    Lambda<-((-grid_size):grid_size)/(grid_size+1)
    # remove zero  
    Lambda<-Lambda[-(grid_size+1)]
  }
  rho_mat<-NULL
  crit_rhoyy<-10000
  crit_rhoyz<--2
# M is used for latest solution: lag-one autocovaraince generating function  
  M<-matrix(nrow=L,ncol=L)
  M[L,]<-rep(0,L)
  M[L-1,]<-c(rep(0,L-1),0.5)
  for (i in 1:(L-2))
    M[i,]<-c(rep(0,i),0.5,rep(0,L-1-i))
  M<-M+t(M)
  eigen_obj<-eigen(M)
  maxrho<-max(eigen_obj$values)
# Orthonormal basis of eigenvectors of M  
  V<-eigen_obj$vectors
  nu_vec<-NULL
  if (forecast_horizon==0)
  {  
    gammak_target<-gammak_generic
  } else
  {  
    if (forecast_horizon<length(gammak_generic))
    {  
      gammak_target<-c(gammak_generic[(1+forecast_horizon):length(gammak_generic)],rep(0,forecast_horizon))
    } else
    {
      print("Forecast horizon is larger than length of target: ill-posed problem")
      return()
    }
  } 
  # Add zeroes to target if it is shorter than L  
  if (length(gammak_target)<L)
  {
    gammak_target<-c(gammak_target,rep(0,L-length(gammak_target)))
  }
# Truncate target if it is longer than L  
  if (length(gammak_target)>L)
  {
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print(paste("Length of target ",length(gammak_target)," exceeds L=",L,". Target will be truncated",sep=""))
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    gammak_target<-gammak_target[1:L]  
  }
# Spectral weights of target/MSE  
  w<-solve(V)%*%gammak_target
  
  
  for (i in 1:length(Lambda))#i<-500  #i<-length(Lambda) i<-900
  {
    lambda1<-Lambda[i]
# Multiply with maxrho to have |nu|>2*maxrho instead of |nu|>2 (proof in paper applies for  |nu|>2*maxrho)      
    nu<-(lambda1+1/lambda1)*maxrho
    nu_vec<-c(nu_vec,nu)
    if (F)
    {  
# Old code: numerically inefficient since a matrix inversion is required for each grid-point    
      Nu<-2*M-nu*diag(rep(1,L))
      bk_new<-solve(Nu)%*%gammak_target[1:L]
    } else
    {  
# Frequency-domain convolution for bk: much faster
      bk_new<-(V)%*%(w/(2*eigen_obj$values-nu))
    }
    if (F)
    { 
# Possible checks: discarded when running optimization      
      ts.plot(bk_new)
      k2<-5
      bk_new[k2]-(lambda1+1/lambda1)*bk_new[k2-1]+bk_new[k2-2]-gammak_target[k2-1]
      bk_new/gammak_target
      compute_holding_time_func(bk_new)$rho_ff1
    } 
    rho_yy_best<-compute_holding_time_func(bk_new)$rho_ff1[1,1]
    rho_yz_best<-(t(bk_new)%*%gammak_target[1:length(bk_new)])[1,1]/(sqrt(t(bk_new)%*%bk_new)*sqrt(gammak_target[1:length(bk_new)]%*%gammak_target[1:length(bk_new)]))[1,1]
    rho_mat<-rbind(rho_mat,t(c(rho_yy_best,rho_yz_best)))
# Change sign if correlation with target is negative: the sign of the optimal solution bk below is changed accordingly   
    if (rho_yz_best<0)
    {  
      rho_mat[nrow(rho_mat),2]<--rho_mat[nrow(rho_mat),2]
      bk_new<--bk_new
      rho_yz_best<--rho_yz_best
    }
# Track criterion value for each grid-point and select grid-point which maximizes criterion    
    crit<-abs(rho_yy_best-rho0)
    if (crit<crit_rhoyy)
    {
      bk_best<-bk_new
      crit_rhoyy<-crit
      crit_rhoyz<-rho_yz_best
      lambda_opt<-lambda1
      i_select<-i
    }
  }  
  colnames(rho_mat)<-c("rho_yy_best","rho_yz_best")

  return(list(rho_mat=rho_mat,bk_best=bk_best,crit_rhoyy=crit_rhoyy,
              crit_rhoyz=crit_rhoyz,lambda_opt=lambda_opt,nu_vec=nu_vec))
}










# This function applies filters to data
# The function filt_func can deal with double, matrix or xts objects
SSA_filter_func<-function(filter_mat,L,x)#x<-series
{
  
  y_mat_indicatorh<-NULL
  for (j in 1:ncol(filter_mat))#j<-2
  {
    b<-filter_mat[1:L,j]#/sqrt(var(filter_mat_select[1:L_gdp,j]))
    
    filt_obj<-filt_func(x,b)
    #  if (j==1)
    #    filt_obj$yhat<-filt_obj$yhat+0.01*sqrt(var(filter_mat_select[1:L,j]))
    
    y_mat_indicatorh<-cbind(y_mat_indicatorh,filt_obj$yhat)
    
  }
  
  y_mat<-y_mat_indicatorh[L:nrow(y_mat_indicatorh),]
  colnames(y_mat)<-colnames(filter_mat)
  
  
  return(list(y_mat=y_mat))
} 



# Filter function: applies a filter b to a series x which can be xts or double
#   If x is xts then time ordering of b is reversed
filt_func<-function(x,b)
{
  L<-length(b)
  yhat<-x
  if (is.matrix(x))
  {  
    length_time_series<-nrow(x)
  } else
  {
    if (is.vector(x))
    {
      length_time_series<-length(x)
    } else
    {
      print("Error: x is neither a matrix nor a vector!!!!")
    }
  }
  for (i in L:length_time_series)
  {
    # If x is an xts object then we cannot reorder x in desceding time i.e. x[i:(i-L+1)] is the same as  x[(i-L+1):i]
    #   Therefore, in this case, we have to revert the ordering of the b coefficients.    
    if (is.xts(x))
    {
      yhat[i]<-as.double(b[L:1]%*%x[i:(i-L+1)])#tail(x) x[(i-L+1):i]
    } else
    {
      yhat[i]<-as.double(b%*%x[i:(i-L+1)])#tail(x) x[(i-L+1):i]
    }
  }
  #  names(yhat)<-index(x)#index(yhat)  index(x)
  #  yhat<-as.xts(yhat,tz="GMT")
  return(list(yhat=yhat))
}





#----------------------------------------------------------------------------
# Timeliness function: compute peak correlation and tau-statistic, see appendix in paper
# It relies on new_lead_at_crossing_func for computing the time-shift
compute_timeliness_func<-function(mplot,max_lead,ht,last_crossing_or_closest_crossing,outlier_limit)
{  
  
  # Peak correlation
  
  cor_peak<-NULL
  for (i in 1:max_lead)
  {
    cor_peak<-c(cor_peak,cor(mplot[i:(nrow(mplot)),2],mplot[1:(nrow(mplot)-i+1),1]))
  }
  # Invert time ordering
  cor_peak<-cor_peak[max_lead:1]
  # Compute other tail
  for (i in 1:(max_lead-1))
  {
    cor_peak<-c(cor_peak,cor(mplot[(i+1):(nrow(mplot)),1],mplot[1:(nrow(mplot)-i),2]))
  }
  
  
  plot(cor_peak,col="blue",main="Peak correlations",axes=F,type="l", xlab="Lead/lag",ylab="Correlation")
  abline(v=which(cor_peak==max(cor_peak)),col="blue")
  at_vec<-c(1,11,21,31,41,51,61,71,81)
  axis(1,at=at_vec,labels=at_vec-max_lead)
  axis(2)
  box()
  
  #------------------------------------------------------------
  # Empirical lead/lag at zero-crossings
  # Skip all crossings with lead/lag>outlier_limit
  skip_larger<-outlier_limit
  # Index of series with more crossings: this is measured against the crossings of the reference series
  con_ind<-2
  # Index of reference series: this one has less crossings and shift is measured with reference to thse crossings only
  ref_ind<-1
  # Select closest crossing of same sign (last_crossing_or_closest_crossing<-F) or 
  #   last crossing of same sign in a vicinity of reference crossing (last_crossing_or_closest_crossing<-T)
  # The setting last_crossing_or_closest_crossing<-T is closer to applications though still a bit optimistic #   because one doesn't know that a particular crossing will be the last in the vicinity
  # The setting last_crossing_or_closest_crossing<-F is unrealistic since the contender filter might generate additional noisy crossings after the closest one 
  if (last_crossing_or_closest_crossing)
  {
    # Size of vicinity to look for turning-point: +/- vicinity around a reference crossing: one picks the last
    #  (of correct sign) in this vicinity
    # Select equal to holding-time (beyond that point signs could change, in the mean)  
    vicinity<-ht
  } else
  {
    vicinity<-NULL
  }
  
  
  lead_lag_cross_obj<-new_lead_at_crossing_func(ref_ind,con_ind,mplot,last_crossing_or_closest_crossing,vicinity)
  
  number_cross<-lead_lag_cross_obj$number_crossings_per_sample
  colnames(number_cross)<-colnames(mplot)[c(con_ind,ref_ind)]
  # Summands of Tau statistic in paper  
  tau_vec<-c(lead_lag_cross_obj$cum_ref_con[1],diff(lead_lag_cross_obj$cum_ref_con))
  remove_tp<-which(abs(tau_vec)>skip_larger)
  if (length(remove_tp)>0)
  {  
    tau_vec_adjusted<-tau_vec[-remove_tp]
  } else
  {
    print("no outlier adjustment necessary in tau-statistic")    
    tau_vec_adjusted<-tau_vec
  }
  # Positive drift i.e. lead of SSA filter
  ts.plot(cumsum(tau_vec_adjusted))
  ts.plot(cumsum(tau_vec))
  # Tau-statistic: mean lead (positive) or lag (negative) of reference filter: with outlier removal
  tau_adjusted<-mean(tau_vec_adjusted)
  tau_adjusted
  # Shift without outlier removal
  tau<-lead_lag_cross_obj$mean_lead_ref_con
  # Test for significance of shift
  t_conf_level<-t.test(tau_vec_adjusted,  alternative = "two.sided")$p.value
  # Strongly significant lead
  t_conf_level
  t_test_adjusted<-t.test(tau_vec_adjusted,  alternative = "two.sided")$statistic
  t_test<-t.test(tau_vec,  alternative = "two.sided")$statistic
  
  return(list(cor_peak=cor_peak,tau_vec=tau_vec,tau_vec_adjusted=tau_vec_adjusted,tau=tau,tau_adjusted=tau_adjusted,t_test=t_test,t_test_adjusted=t_test_adjusted,number_cross=number_cross))
  
}


# New lead at crossing function 
#   This function implements the mean-shift statistic in paper
#   It accounts for sign of crossings (up-swing/down-swing)
# There is a new additional feature when compared with earlier version above
#   1. If last_crossing_or_closest_crossing==F then it replicates the function used and described in paper
#   2. If last_crossing_or_closest_crossing==T then for each up- or down-turn of the reference filter 
#     one selects the last up- or down-turn of contender in a vicinity of crossing (otherwise the closest will be selected) 
# The setting last_crossing_or_closest_crossing==T is more realistic but still a bit optimistic because 
#   the user generally does not know that a particular crossing will be the last one.
# Results depend on vicinity which is the size of the neighborhood +/-vicinity about crossing
#   Can be chosen with a link to holding time: if holding time is large, then vicinity could be larger, too
# In general last_crossing_or_closest_crossing==F will lead to smaller lead-times (closest crossing) than
#   last_crossing_or_closest_crossing==T (latest crossing in vicinity)
# In the paper we use last_crossing_or_closest_crossing==F exclusively (closest crossings)
new_lead_at_crossing_func<-function(ref_ind,con_ind,mplot,last_crossing_or_closest_crossing,vicinity)
{
  ref_cross<-which(sign(mplot[1:(nrow(mplot)-1),ref_ind])!=sign(mplot[2:(nrow(mplot)),ref_ind]))
  con_cross<-which(sign(mplot[1:(nrow(mplot)-1),con_ind])!=sign(mplot[2:(nrow(mplot)),con_ind]))
  if (mplot[ref_cross[1],ref_ind]<0)
  {
    ref_cross_sign<-(-1)^(0:(length(ref_cross)-1))*ref_cross  
  } else
  {
    ref_cross_sign<--1*((-1)^(0:(length(ref_cross)-1))*ref_cross)
  }
  if (mplot[con_cross[1],con_ind]<0)
  {
    con_cross_sign<-(-1)^(0:(length(con_cross)-1))*con_cross  
  } else
  {
    con_cross_sign<--1*((-1)^(0:(length(con_cross)-1))*con_cross)
  }
  # Crossings from above and from below  
  con_cross_sign_plus<-con_cross_sign[con_cross_sign>0]
  con_cross_sign_negative<-con_cross_sign[con_cross_sign<0]
  
  disc_vec<-NULL
  # Check if there are crossings    
  if (length(ref_cross)>0&length(con_cross)>0)
  {
    
    
    ref_len<-length(ref_cross)#length(con_cross)
    
    # For the j-th zero-crossing of reference
    for (j in 1:ref_len)#j=12
    {
      if (ref_cross_sign[j]>0)
      {
        # Up-turns      
        if (!last_crossing_or_closest_crossing)
        {    
          # For upturns of reference: select nearest upturn of contender        
          min_i<-min(abs(abs(con_cross_sign_plus)-abs(ref_cross[j])))
          which_min<-which(abs(abs(con_cross_sign_plus)-abs(ref_cross[j]))==min_i)
          disc_vec<-c(disc_vec,max(abs(con_cross_sign_plus[which_min]))-abs(ref_cross[j]))
        } else
        {
          # For upturns of reference: select latest upturn of contender in vicinity of reference up-turn
          # This is closer to real-time application though it is still optimistic because the user doesn't know 
          #   that this will be the last up-turn in practice      
          # 1 All upturns in vicinity of reference up-turn        
          rert<-which(abs(abs(con_cross_sign_plus)-abs(ref_cross[j]))<vicinity)
          # 2. Select latest one or border of vicinity
          if (length(rert)>0)
          {  
            # If there is a TP in vicinity: select last one            
            max_rert<-max(rert)
            last_tp<-abs(con_cross_sign_plus[max_rert])
          } else
          {
            # Otherwise select border of vicinity (one again a bit optimistic)            
            last_tp<-abs(ref_cross[j])+vicinity
          }
          disc_vec<-c(disc_vec,last_tp-abs(ref_cross[j]))
          
        }  
      } else
      {
        # For downturns      
        if (!last_crossing_or_closest_crossing)
        {  
          # For downturns of reference: select nearest downturn of contender        
          min_i<-min(abs(abs(con_cross_sign_negative)-abs(ref_cross[j])))
          which_min<-which(abs(abs(con_cross_sign_negative)-abs(ref_cross[j]))==min_i)
          # Time-difference at down-turn: 
          #   in case of two possible nearest crossings we select the max (because the contrary direction is on at the time point of the crossing of the reference filter)
          disc_vec<-c(disc_vec,max(abs(con_cross_sign_negative[which_min]))-abs(ref_cross[j]))
          
        } else
        {
          # For down-turns of reference: select latest down-turn of contender in vicinity of reference down-turn
          # This is closer to real-time application though it is still optimistic because the user doesn't know 
          #   that this will be the last down-turn in practice      
          # 1 All down-pturns in vicinity of reference down-turn        
          rert<-which(abs(abs(con_cross_sign_negative)-abs(ref_cross[j]))<vicinity)
          # 2. Select latest one or border of vicinity
          if (length(rert)>0)
          {  
            # If there is a TP in vicinity: select last one            
            max_rert<-max(rert)
            last_tp<-abs(con_cross_sign_negative[max_rert])
          } else
          {
            # Otherwise select border of vicinity (one again a bit optimistic)            
            last_tp<-abs(ref_cross[j])+vicinity
          }
          disc_vec<-c(disc_vec,last_tp-abs(ref_cross[j]))
        }  
      }
    }
    # Total number of crossings of contender and of reference
    number_crossings_per_sample<-c(length(con_cross),length(ref_cross))
  } else
  {
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("no zero-crossings observed")
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    number_crossings_per_sample<-rbind(number_crossings_per_sample,c(0,0))
  }
  if (length(names(ref_cross))>0)
    names(disc_vec)<-names(ref_cross)
  number_crossings_per_sample<-t(as.matrix(number_crossings_per_sample,nrow=1))
  cum_ref_con<-cumsum(na.exclude(disc_vec))
  #  ts.plot(cum_ref_con)
  mean_lead_ref_con<-cum_ref_con[length(cum_ref_con)]/length(cum_ref_con)
  mean_lead_ref_con
  tail(cum_ref_con)
  return(list(cum_ref_con=cum_ref_con,mean_lead_ref_con=mean_lead_ref_con,number_crossings_per_sample=number_crossings_per_sample))
}






#------------------------------------------------------------------------
# Plots for BCA-section in paper
plot_paper<-function(y_mat,start_date,end_date,colo_all)
{  
  if (is.null(end_date))
  { 
    # Last data point    
    end_date<-format(Sys.time(), "%Y-%m-%d")
  } 
  if (is.null(start_date))
  { 
    # Last data point    
    start_date<-index(y_mat)[1]
  } 
  select_output<-c(2,4)
  mplot<-scale(y_mat[paste(start_date,"/",end_date,sep=""),select_output],center=F,
               rep(T,length(select_output)))
  coli<-colo_all[select_output]
  q_gap<-plot(mplot,main=paste(colnames(y_mat)[select_output[1]], " vs. ",colnames(y_mat)[select_output[2]],sep=""),col=coli)
  #p<-mtext("HP-gap",col=coli[1],line=-3)
  sel_tp<-c(1,2)
  for (i in 1:length(sel_tp))
  {
    ret<-mplot[,sel_tp[i]]
    tp_last<-index(mplot)[which(sign(ret)!=lag(sign(ret)))]
    events<-xts(rep("",length(tp_last)),tp_last)
    q_gap<-addEventLines(events, srt=90, pos=2,col=coli[sel_tp[i]])
  }
  q_gap
  x_gap<-nber_dates_polygon(start_date,mplot)$x
  y_gap<-nber_dates_polygon(start_date,mplot)$y
  q_gap
  polygon(x_gap, y_gap, xpd = T, col = "grey",density=10)#
  
  
  select_output<-c(1,4)
  mplot<-scale(y_mat[paste(start_date,"/",end_date,sep=""),select_output],center=F,scale=rep(T,length(select_output)))
  coli<-colo_all[select_output]
  q_trend<-plot(mplot,main=paste(colnames(y_mat)[select_output[1]], " vs. ",colnames(y_mat)[select_output[2]],sep=""),col=coli)
  #p<-mtext("HP-gap",col=coli[1],line=-3)
  sel_tp<-c(1,2)
  for (i in 1:length(sel_tp))
  {
    ret<-mplot[,sel_tp[i]]
    tp_last<-index(mplot)[which(sign(ret)!=lag(sign(ret)))]
    events<-xts(rep("",length(tp_last)),tp_last)
    q_trend<-addEventLines(events, srt=90, pos=2,col=coli[sel_tp[i]])
  }
  x_trend<-nber_dates_polygon(start_date,mplot)$x
  y_trend<-nber_dates_polygon(start_date,mplot)$y
  q_trend
  polygon(x_trend, y_trend, xpd = T, col = "grey",density=10)#
  
  select_output<-c(4,3)
  mplot<-scale(y_mat[paste(start_date,"/",end_date,sep=""),select_output],center=F,scale=rep(T,length(select_output)))
  coli<-colo_all[select_output]
  q_SSA<-plot(mplot,main=paste(colnames(y_mat)[select_output[1]], " vs. ",colnames(y_mat)[select_output[2]],sep=""),col=coli)
  #p<-mtext("HP-gap",col=coli[1],line=-3)
  sel_tp<-c(1,2)
  for (i in 1:length(sel_tp))
  {
    ret<-mplot[,sel_tp[i]]
    tp_last<-index(mplot)[which(sign(ret)!=lag(sign(ret)))]
    events<-xts(rep("",length(tp_last)),tp_last)
    q_SSA<-addEventLines(events, srt=90, pos=2,col=coli[sel_tp[i]])
  }
  q_SSA
  x_SSA<-nber_dates_polygon(start_date,mplot)$x
  y_SSA<-nber_dates_polygon(start_date,mplot)$y
  q_SSA
  polygon(x_SSA, y_SSA, xpd = T, col = "grey",density=10)#
  
  
  
  return(list(q_gap=q_gap,q_trend=q_trend,q_SSA=q_SSA,x_trend=x_trend,y_trend=y_trend,x_gap=x_gap,y_gap=y_gap,x_SSA=x_SSA,y_SSA=y_SSA))
}


# Adds shaded areas corresponding to NBER-recession datings
nber_dates_polygon<-function(start_date,mat)
{
  dat<-NULL
  for (i in 1:nrow(nberDates()))
  {
    dat<-c(dat,paste(substr(nberDates()[i,1],1,4),"-",substr(nberDates()[i,1],5,6),"-",substr(nberDates()[i,1],7,8),sep=""))
  }
  
  
  starting_date<-as.POSIXct(strptime(dat, "%Y-%m-%d"),tz="UTC")
  
  starting_date<-starting_date[which(starting_date>start_date)]
  
  dat<-NULL
  for (i in 1:nrow(nberDates()))
  {
    dat<-c(dat,paste(substr(nberDates()[i,2],1,4),"-",substr(nberDates()[i,2],5,6),"-",substr(nberDates()[i,2],7,8),sep=""))
  }
  
  ending_date<-as.POSIXct(strptime(dat, "%Y-%m-%d"),tz="UTC")
  ending_date<-ending_date[which(ending_date>start_date)]
  
  x<-y<-NULL
  for (i in 1:length(starting_date))
  {
    x<-c(x,starting_date[i],starting_date[i],ending_date[i],ending_date[i])
    y<-c(y,min(na.exclude(mat)),max(na.exclude(mat)),max(na.exclude(mat)),min(na.exclude(mat)))
  }
  return(list(x=x,y=y))
}


#-------------------------------------


# Computes amplitude and time shifts (mainly for illustration purposes)
amp_shift_func<-function(K,b,plot_T)
{
  #  if (sum(b)<0)
  #  {
  #    print("Sign of coefficients has been changed")
  #    b<-b*sign(sum(b))
  #  }
  omega_k<-(0:K)*pi/K
  trffkt<-0:K
  for (i in 0:K)
  {
    trffkt[i+1]<-b%*%exp(1.i*omega_k[i+1]*(0:(length(b)-1)))
  }
  amp<-abs(trffkt)
  shift<-Arg(trffkt)/omega_k
  shift[1]<-sum((0:(length(b)-1))*b)/sum(b)
  if (plot_T)
  {
    par(mfrow=c(2,1))
    plot(amp,type="l",axes=F,xlab="Frequency",ylab="Amplitude",main="Amplitude")
    axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
    plot(shift,type="l",axes=F,xlab="Frequency",ylab="Shift",main="Shift",ylim=c(min(min(shift,na.rm=T),0),max(shift,na.rm=T)))
    axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }  
  return(list(trffkt=trffkt,amp=amp,shift=shift))
}


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Older code: numerically much less efficient
# This function finds optimal lambda (note that nu=lambda+1/lambda)) through grid search based on function compute_bk_from_ma_expansion_with_boundary_constraint_func above
# It implements solution in corollary 1 of paper
find_lambda1_subject_to_holding_time_constraint_func_old<-function(grid_size,L,gammak_generic,rho0,forecast_horizon,with_negative_lambda)
{
  Lambda<-(1:grid_size)/(grid_size+1)
  if (with_negative_lambda)
  { 
    # Allow negative lambda1 too    
    Lambda<-((-grid_size):grid_size)/(grid_size+1)
    # remove zero  
    Lambda<-Lambda[-(grid_size+1)]
  }
  rho_mat<-NULL
  crit_rhoyy<-10000
  crit_rhoyz<--2
  # M is used for latest solution: lag-one autocovaraince generating function  
  M<-matrix(nrow=L,ncol=L)
  M[L,]<-rep(0,L)
  M[L-1,]<-c(rep(0,L-1),0.5)
  for (i in 1:(L-2))
    M[i,]<-c(rep(0,i),0.5,rep(0,L-1-i))
  M<-M+t(M)
  maxrho<-max(eigen(M)$values)
  
  nu_vec<-NULL
  if (forecast_horizon==0)
  {  
    gammak_target<-gammak_generic
  } else
  {  
    gammak_target<-c(gammak_generic[(1+forecast_horizon):length(gammak_generic)],rep(0,forecast_horizon))
  } 
  # Add zeroes to target if it is shorter than L  
  if (length(gammak_target)<L)
  {
    gammak_target<-c(gammak_target,rep(0,L-length(gammak_target)))
  }
  if (length(gammak_target)>L)
  {
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print(paste("Lenhth of target ",length(gammak_target)," exceeds L=",L,". Target will be shortened",sep=""))
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    gammak_target<-gammak_target[1:L]  
  }  
  
  for (i in 1:length(Lambda))#i<-1  #i<-length(Lambda) i<-900
  {
    lambda1<-Lambda[i]
    
    # Latest solution based on solve(Nu)
    # Multiply with maxrho to have |nu|>2*maxrho instead of |nu|>2 (proof in paper applies for  |nu|>2*maxrho)      
    nu<-(lambda1+1/lambda1)*maxrho
    nu_vec<-c(nu_vec,nu)
    Nu<-2*M-nu*diag(rep(1,L))
    bk_new<-solve(Nu)%*%gammak_target[1:L]
    if (F)
    { 
      ts.plot(bk_new)
      bk_new[k2]-(lambda1+1/lambda1)*bk_new[k2-1]+bk_new[k2-2]-gammak_target[k2-1]
      bk_new/gammak_target
      compute_holding_time_func(bk_new)$rho_ff1
    } 
    #    ts.plot(bk_new)
    rho_yy_best<-compute_holding_time_func(bk_new)$rho_ff1[1,1]
    rho_yz_best<-(t(bk_new)%*%gammak_target[1:length(bk_new)])[1,1]/(sqrt(t(bk_new)%*%bk_new)*sqrt(gammak_target[1:length(bk_new)]%*%gammak_target[1:length(bk_new)]))[1,1]
    
    rho_mat<-rbind(rho_mat,t(c(rho_yy_best,rho_yz_best)))
    # Change sign if correlation with target is negative: the sign of the optimal solution bk below is changed accordingly   
    if (rho_yz_best<0)
    {  
      rho_mat[nrow(rho_mat),2]<--rho_mat[nrow(rho_mat),2]
      bk_new<--bk_new
      rho_yz_best<--rho_yz_best
    }
    crit<-abs(rho_yy_best-rho0)
    if (crit<crit_rhoyy)
    {
      bk_best<-bk_new
      crit_rhoyy<-crit
      crit_rhoyz<-rho_yz_best
      lambda_opt<-lambda1
      i_select<-i
    }
    
  }  
  colnames(rho_mat)<-c("rho_yy_best","rho_yz_best")
  crit_rhoyy
  crit_rhoyz
  #  ts.plot(bk)
  ts.plot(bk_best)
  
  return(list(rho_mat=rho_mat,bk_best=bk_best,crit_rhoyy=crit_rhoyy,
              crit_rhoyz=crit_rhoyz,lambda_opt=lambda_opt,nu_vec=nu_vec))
}






