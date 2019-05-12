

# Computes regularization matrices: 
# The new function calls Q_reg_func which initializes the bilinear forms.

reg_mat_func<-function(weight_h_exp,L,c_eta,lambda_decay,Lag,grand_mean,lambda_cross,lambda_smooth,lag_mat)
{
  lambda_smooth<-100*tan(min(abs(lambda_smooth),0.999999999)*pi/2)
  lambda_cross<-100*tan(min(abs(lambda_cross),0.9999999999)*pi/2)
  if (c_eta)
  {  
    # Chris replication: imposing non-linear transform to lambda_decay[1] too
    lambda_decay<-c(min(abs(tan(min(abs(lambda_decay[1]),0.999999999)*pi/2)),1),100*tan(min(abs(lambda_decay[2]),0.999999999)*pi/2))
  } else
  {  
    lambda_decay<-c(min(abs(lambda_decay[1]),1),100*tan(min(abs(lambda_decay[2]),0.999999999)*pi/2))
  }
# New regularization function

  Q_obj<-Q_reg_func(L,weight_h_exp,lambda_decay,lambda_smooth,lambda_cross,Lag,lag_mat)

  Q_smooth<-Q_obj$Q_smooth
  Q_decay<-Q_obj$Q_decay
  Q_cross<-Q_obj$Q_cross
# The Q_smooth and Q_decay matrices address regularizations for original unconstrained parameters (Therefore dimension L^2)
# At the end, the matrix des_mat is used to map these regularizations to central-deviance parameters
# accounting for first order constraints!
  
  return(list(Q_smooth=Q_smooth,Q_cross=Q_cross,Q_decay=Q_decay,lambda_decay=lambda_decay,
              lambda_smooth=lambda_smooth,lambda_cross=lambda_cross))
}









Q_reg_func<-function(L,weight_h_exp,lambda_decay,lambda_smooth,lambda_cross,Lag,lag_mat)
{
  Q_smooth<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  Q_decay<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  # Cross-sectional regularization if dimension>1
  if ((length(weight_h_exp[1,])>1))
  {
    # The cross-sectional regularization is conveniently implemented on central-deviance parameters. The regularization is expressed on the
    # unconstrained central-deviance parameters (dimension L), then mapped to the original (unconstrained) parameters (dimension L) with Q_centraldev_original
    # and then maped back to central-deviance with constraint (dim L-1) with des_mat (mathematically unnecessarily complicate but more convenient to implement in code).
    Q_cross<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
    Q_centraldev_original<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  } else
  {
    # 16.08.2012
    Q_cross<-NULL
  }
  for (i in 1:L)
  {
  # For symmetric filters or any historical filter with Lag>0 the decay must be symmetric about b_max(0,Lag) 
  # lambda_decay is a 2-dim vector: the first component controls for the exponential decay and the second accounts for the strength of the regularization
  # with mixed-frequency data the decay should account for the effective lag (of the lower frequency data) as measured on the high-frequency scale
  # The maximal weight is limited to 1e+10
#    Q_decay[i,i]<-(1+lambda_decay[1])^(2*abs(i-1-max(0,Lag)))
    Q_decay[i,i]<-min((1+lambda_decay[1])^(2*abs(lag_mat[i,1]-max(0,Lag))),1e+4)
    
    if (L>4)
    {
      if(i==1)
      {
        Q_smooth[i,i:(i+2)]<-c(1,-2,1)
      } else
      {
        if(i==2)
        {
          Q_smooth[i,(i-1):(i+2)]<-c(-2,5,-4,1)
        } else
        {
          if(i==L)
          {
            Q_smooth[i,(i-2):i]<-c(1,-2,1)
          } else
          {
            if(i==L-1)
            {
              Q_smooth[i,(i-2):(i+1)]<-c(1,-4,5,-2)
            } else
            {
              Q_smooth[i,(i-2):(i+2)]<-c(1,-4,6,-4,1)
            }
          }
        }
      }
    } else
    {
      print("L<=4: no smoothness regularization imposed!!!!!!!!!")
    }
  }
  
  if (length(weight_h_exp[1,])>1)
  {
    
    for (j in 1:max(1,(length(weight_h_exp[1,])-1)))   #j<-1
    {
      if (L>4)
        Q_smooth[j*L+1:L,j*L+1:L]<-Q_smooth[1:L,1:L]
      for (i in 1:L)
        Q_decay[j*L+i,j*L+i]<-min((1+lambda_decay[1])^(2*abs(lag_mat[i,1+j]-max(0,Lag))),1e+4)

    }
    Q_centraldev_original<-diag(rep(1,L*length(weight_h_exp[1,])))
    if (L>1)
    {  
      diag(Q_centraldev_original[1:L,L+1:L])<-rep(-1,L)
      for (i in 2:length(weight_h_exp[1,]))   #i<-2
      {
        diag(Q_centraldev_original[(i-1)*L+1:L,1:L])<-rep(1,L)
        diag(Q_centraldev_original[(i-1)*L+1:L,(i-1)*L+1:L])<-rep(1,L)
        diag(Q_centraldev_original[1:L,(i-1)*L+1:L])<-rep(-1,L)
      }
    }
    Q_centraldev_original<-solve(Q_centraldev_original)
    
    # 06.08.2012: the following if allows for either grand-mean parametrization (I-MDFA version prior to 30.07.2012) or original parameters
    #   If grand_mean==T then the code replicates I-MDFA as released prior 30.07.2012
    #   If grand_mean==F then the new parametrization is used.
    #   Differences between both approaches: see section 7.2 of my elements paper posted on SEFBlog (both codes are identical when no regularization is imposed. Otherwise the later version (grand_mean==F) is logically more consistent becuase it treats all series identically (no asymmetry)).
    if (grand_mean)
    {
      diag(Q_cross[L+1:((length(weight_h_exp[1,])-1)*L),L+1:((length(weight_h_exp[1,])-1)*L)])<-
        rep(1,((length(weight_h_exp[1,])-1)*L))
    } else
    {
      #30.07.2012:new definition (parametrization) of Q_cross (Lambda_{cross} in the elements-paper)
      diag(Q_cross)<-1
      for (i in 1:length(weight_h_exp[1,]))
      {
        for (j in 1:L)
        {
          Q_cross[(i-1)*L+j,j+(0:(length(weight_h_exp[1,])-1))*L]<-Q_cross[(i-1)*L+j,j+(0:(length(weight_h_exp[1,])-1))*L]-1/length(weight_h_exp[1,])
        }
      }
    }
  } else
  {
# define matrix for univariate case
    Q_centraldev_original<-NULL
  }
# Normalizing the troika are new: disentangle the effect by L
  Q_decay<-Q_decay*lambda_decay[2]
  Q_cross<-Q_cross*lambda_cross                   #Qh<-Q_cross
  Q_smooth<-Q_smooth*lambda_smooth
  if (lambda_decay[2]>0)
  {
# The second parameter in lambda_decay accounts for the strength of the regularization
    Q_decay<-lambda_decay[2]*(Q_decay/(sum(diag(Q_decay))))
  }
  if (lambda_cross>0)
  {
    if (c_eta)
    {  
      Q_cross<-lambda_cross^2*(Q_cross/(sum(diag(Q_cross))))
    } else
    {
      Q_cross<-lambda_cross*(Q_cross/(sum(diag(Q_cross))))
    }
  }
  if (lambda_smooth>0&L>4)
  {
    Q_smooth<-lambda_smooth*(Q_smooth/(sum(diag(Q_smooth))))
  }
  return(list(Q_smooth=Q_smooth,Q_decay=Q_decay,Q_cross=Q_cross))
}