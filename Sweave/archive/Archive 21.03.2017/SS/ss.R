
structure_func<-function(weight_structure,weight_func,L,chris_expweight,diff_explanatory,lag_mat,spec_mat)
{
  K<-nrow(weight_func)-1
  if (sum(abs(weight_structure))>0)
  {
    # diff_explanatory<-F implies that DFT is computed as usual (explanatory variables are not differenced: the later is used when targeting trading performances as weighted by diff-filter output)  
    diff_explanatory<-F
    spec_mat_structure<-spec_mat_comp(weight_func,L,-1,chris_expweight,diff_explanatory,lag_mat)$spec_mat     #dim(spec_mat[,1])
    # The target signals are the same for both structural elements: one-step ahead forecasting i.e. the corresponding Gammas are identities
    Gamma_structure<-Gamma_structure_diff<-rep(1,K+1)
    # If the weighting of the investment with the differenced signal is imposed (i.e. if abs(weight_structure)[2]>0) then we must compute the spectral matrix of the differenced series in order to obtain the differenced filter output
    if (abs(weight_structure)[2]>0)
    {
      # diff_explanatory<-T implies that DFT is computed for DIFFERENCED explanatory variables (performance is obtained by weighting investment according to differenced output signal)  
      diff_explanatory<-T  
      spec_mat_structure_diff<-spec_mat_comp(weight_func,L,-1,chris_expweight,diff_explanatory,lag_mat)$spec_mat     #dim(spec_mat[,1])  
    } else
    {
      # If abs(weight_structure)[2]==0 then we may fill the spectrum with any entry since the weight will be zero anyway      
      spec_mat_structure_diff<-spec_mat_structure
    }
  } else
  {
    # We have to initialize spect_mat_structure and spec_mat_structure_diff by an arbitrary matrix of right dimensions    
    spec_mat_structure<-spec_mat_structure_diff<-spec_mat    
    Gamma_structure<-Gamma_structure_diff<-rep(0,K+1)
  }
  return(list(Gamma_structure=Gamma_structure,spec_mat_structure=spec_mat_structure,
              Gamma_structure_diff=Gamma_structure_diff,spec_mat_structure_diff=spec_mat_structure_diff        ))
}







white_noise_synchronicity<-function(weight_func,white_noise,synchronicity)
{
  if (white_noise)
  {
    # In this case we assume that all DFTs are flat but we keep the phase information intact    
    means<-apply(abs(weight_func),2,mean)
    for (j in 1:ncol(weight_func))
    {
      weight_func[,j]<-exp(1.i*Arg(weight_func[,j]))*means[j]  
    }
  }  #ts.plot(Im(weight_func),lty=1:4)
  if (synchronicity)
  {
    # In this case we assume that all series are synchronized at all frequencies but we keep the amplitudes intact    
    weight_func<-abs(weight_func)
  }
  return(list(weight_func=weight_func))
}











reg_mat_func<-function(weight_h_exp,L,chris_expweight,lambda_decay,Lag,grand_mean,lambda_cross,lambda_smooth)
{
  lambda_smooth<-100*tan(min(abs(lambda_smooth),0.999999999)*pi/2)
  lambda_cross<-100*tan(min(abs(lambda_cross),0.9999999999)*pi/2)
  if (chris_expweight)
  {  
    # Chris replication: imposing non-linear transform to lambda_decay[1] too        
    lambda_decay<-c(min(abs(tan(min(abs(lambda_decay[1]),0.999999999)*pi/2)),1),100*tan(min(abs(lambda_decay[2]),0.999999999)*pi/2))
  } else
  {  
    lambda_decay<-c(min(abs(lambda_decay[1]),1),100*tan(min(abs(lambda_decay[2]),0.999999999)*pi/2))
  }
  
  
  # The smoothness and decay regularization are conveniently (rightly) implemented on original parameters
  # The Q_smooth and Q_decay matrices address regularizations for original unconstrained parameters (Therefore dimension L^2)
  # At the end, the matrix des_mat is used to map these regularizations to central-deviance parameters
  # accounting for first order constraints!!!!!!!!!!!!!!!!!!!
  Q_smooth<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  Q_decay<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  # Cross-sectional regularization if dimension>1
  if ((length(weight_h_exp[1,])>1))
  {
    # The cross-sectional regularization is conveniently implemented on central-deviance parameters. The regularization is expressed on the
    # unconstrained central-deviance parameters (dimension L), then mapped to the original (unconstrained) parameters (dimension L) with Q_centraldev_original
    # and then maped back to central-deviance with constraint (dim L-1) with des_mat (mathematically unnecessarily complicate but more convenient to implement in code).
    Q_cross<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  } else
  {
    # 16.08.2012
    Q_cross<-NULL
  }
  
  for (i in 1:L)
  {
    # For symmetric filters or any historical filter with Lag>0 the decay must be symmetric about b_max(0,Lag) zunehmen
    
    if (chris_expweight)
    {  
      # Chris' modification
      Q_decay[i,i]<-(1+lambda_decay[1])^(2*abs(i-1-Lag))
    } else
    {
      # New 07.08.2012: lambda_decay is now a 2-dim vector: the first component controls for the exponential decay and the second accounts for the strength of the regularization
      Q_decay[i,i]<-(1+lambda_decay[1])^(2*abs(i-1-max(0,Lag)))
    }
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
      Q_decay[j*L+1:L,j*L+1:L]<-Q_decay[1:L,1:L]
    }
    
    
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
  } 
  # New 07.08.2012 : the next lines for normalizing the troika are new: disentangle the effect by L
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
    if (chris_expweight)
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
  return(list(Q_smooth=Q_smooth,Q_cross=Q_cross,Q_decay=Q_decay,lambda_decay=lambda_decay,
              lambda_smooth=lambda_smooth,lambda_cross=lambda_cross))
}




# This one has the speed enhancer
mdfa_analytic_new<-function(K,L,lambda,weight_func,Lag,Gamma,expweight,cutoff,i1,i2,weight_constraint,
                            lambda_cross,lambda_decay,lambda_smooth,lin_expweight,shift_constraint,grand_mean,
                            b0_H0,chris_expweight,weights_only=F,weight_structure,white_noise,
                            synchronicity,lag_mat)
{
  
  
  # 07.05.2014: we rotate the DFT (the rotation below will be redundant (but doesn't hurt neither...))
  weight_target<-weight_func[,1]
  weight_func<-weight_func*exp(-1.i*Arg(weight_target))
  
  white_noise_synchronicity_obj<-white_noise_synchronicity(weight_func,white_noise,synchronicity)
  
  weight_func<-white_noise_synchronicity_obj$weight_func
  # 02.08.2013
  if (i2&i1&any(weight_constraint==0))
  {
    print(rep("!",100))
    print("Warning: i2<-T is not meaningful when i1<-T and weight_constraint=0 (bandpass)")
    print(rep("!",100))
  }
  if (!(length(b0_H0)==L*(dim(weight_func)[2]-1))&length(b0_H0)>0)
    print(paste("length of b0_H0 vector is ",length(b0_H0),": it should be ",L*(dim(weight_func)[2]-1)," instead",sep=""))
  # diff_explanatory<-F implies that DFT is computed as usual (explanatory variables are not differenced: the later is used when targeting trading performances as weighted by diff-filter output)  
  diff_explanatory<-F
  # In order to enhance numerical speed this call could be done outside (as long as L and Lag are fixed)
  # The function spect_mat_comp rotates all DFTs such that first column is real!
  spec_mat<-spec_mat_comp(weight_func,L,Lag,chris_expweight,diff_explanatory,lag_mat)$spec_mat     #dim(spec_mat[,1])
  
  # If abs(weight_structure)>0 then additional structure is instilled into MDFA-criterion: one-step ahead MSE (Lag=-1)
  #   In that case the DFTs are calculated, using Lag=-1. Moreover, the spectrum is not affected by customization (expweight)
  #   We consider two terms: 
  #       1.performance is obtained by weighting investment according to filter signal (weight is the first component weight_structure[1] of the vector weight_structure) 
  #       2.performance is obtained by weighting investment according to differenced filter signal (weight is the second component weight_structure[2] of the vector weight_structure) 
  structure_func_obj<-structure_func(weight_structure,weight_func,L,chris_expweight,diff_explanatory,lag_mat)
  
  Gamma_structure<-structure_func_obj$Gamma_structure
  spec_mat_structure<-structure_func_obj$spec_mat_structure
  Gamma_structure_diff<-structure_func_obj$Gamma_structure_diff
  spec_mat_structure_diff<-structure_func_obj$spec_mat_structure_diff
  
  # weighting of amplitude function in stopband
  omega_Gamma<-as.integer(cutoff*K/pi)
  if ((K-omega_Gamma+1)>0)
  {
    #    lin_expweight <- FALSE
    if (lin_expweight)
    {
      expweight_vec<-c(rep(1,omega_Gamma),1+rep(expweight,K-omega_Gamma+1))
    } else
    {
      if (chris_expweight)
      {
        # Chris mod        
        expweight_vec<-c(rep(1,omega_Gamma+1),(1:(K-omega_Gamma))*pi/K +1)^(expweight/10)
      } else
      {
        expweight_vec<-c(rep(1,omega_Gamma),(1:(K-omega_Gamma+1))^(expweight/2))
      }
      #      expweight_vec<-c(rep(1,omega_Gamma),1+(1:(K-omega_Gamma+1))*(expweight/2))
    }
    weight_h<-weight_func*expweight_vec
  } else
  {
    expweight_vec<-rep(1,K+1)
    weight_h<-weight_func* expweight_vec
  }
  #ts.plot(abs(weight_h))       ts.plot(Gamma)
  # Frequency zero receives half weight
  weight_h[1,]<-weight_h[1,]*ifelse(chris_expweight,1,1/sqrt(2))
  # DFT target variable
  weight_target<-weight_h[,1]
  # Rotate all DFT's such that weight_target is real (rotation does not alter mean-square error). Note that target of strutural 
  # additional structural elements and of original MDFA are the same i.e. rotation must be the same and weight_target are identical too.
  weight_h<-weight_h*exp(-1.i*Arg(weight_target))
  weight_target<-Re(weight_target*exp(-1.i*Arg(weight_target)))
  # If abs(weight_structure)>0 then additional structure is instilled into MDFA-criterion: one-step ahead MSE 
  if (sum(abs(weight_structure))>0)
  {
    # The target weight vector must be free of expweight-weighting (no customization imposed in structure)
    weight_target_structure<-weight_target_structure_diff<-weight_target/expweight_vec
  } else
  {
    weight_target_structure<-weight_target_structure_diff<-rep(0,K+1)
  }
  
  
  # DFT's explaining variables: target variable can be an explaining variable too
  weight_h_exp<-as.matrix(weight_h[,2:(dim(weight_h)[2])])
  
  # The spectral matrix is inflated in stopband: effect of expweight 
  spec_mat<-t(t(spec_mat)*expweight_vec) #dim(spec_mat)  as.matrix(Im(spec_mat[150,]))
  
  # Compute design matrix and regularization matrix
  # new 21.06.2012: new vector shift_constraint in function call
  # 06.08.2012: new parameter grand_mean in function call
  # 26.07.2013: added the variable b0_H0 in the function call which corresponds to the null-hypothesis to which the Troika shrinks solutions
  mat_obj<-mat_func(i1,i2,L,weight_h_exp,lambda_decay,lambda_cross,lambda_smooth,Lag,weight_constraint,shift_constraint,grand_mean,b0_H0,chris_expweight)
  
  des_mat<-mat_obj$des_mat
  reg_mat<-mat_obj$reg_mat           #dim(des_mat[1:23,1:24])
  reg_xtxy<-mat_obj$reg_xtxy
  w_eight<-mat_obj$w_eight
  
  # Solve estimation problem 
  mat_x<-des_mat%*%spec_mat          #dim(mat_x)    dim(spec_mat)   length(w_eight)
  X_new<-t(Re(mat_x))+sqrt(1+Gamma*lambda)*1.i*t(Im(mat_x))
  # xtx can be written either in Re(t(Conj(X_new))%*%X_new) or as below:
  xtx<-t(Re(X_new))%*%Re(X_new)+t(Im(X_new))%*%Im(X_new)#dim(xtx)  diag(xtx) det(xtx)
  # If abs(weight_structure)>0 then additional structure is instilled into MDFA-criterion: one-step ahead MSE
  if (sum(abs(weight_structure))>0)
  {
    X_new_structure<-t(abs(weight_structure[1])*des_mat%*%spec_mat_structure)  # dim(X_new_structure) dim(X_new)
    xtx_structure<-t(Re(X_new_structure))%*%Re(X_new_structure)+t(Im(X_new_structure))%*%Im(X_new_structure)
    
    X_new_structure_diff<-t(abs(weight_structure[2])*des_mat%*%spec_mat_structure_diff)  # dim(X_new_structure) dim(X_new)
    xtx_structure_diff<-t(Re(X_new_structure_diff))%*%Re(X_new_structure_diff)+t(Im(X_new_structure_diff))%*%Im(X_new_structure_diff)
    
    Gamma_structure_diff<-abs(weight_structure)[2]*Gamma_structure    
    Gamma_structure<-abs(weight_structure)[1]*Gamma_structure        
    
  } else
  {
    X_new_structure<-x_new_structure_diff<-0*des_mat%*%spec_mat_structure     
    xtx_structure<-xtx_structure_diff<-0*t(Re(X_new_structure))%*%Re(X_new_structure)+t(Im(X_new_structure))%*%Im(X_new_structure)    
  }
  # The filter restrictions (i1<-T and/or i2<-T) appear as constants on the right hand-side of the equation:
  xtxy<-t(Re(t(w_eight)%*%spec_mat)%*%Re(t(spec_mat)%*%t(des_mat))+
            Im(t(w_eight)%*%t(t(spec_mat)*sqrt(1+Gamma*lambda)))%*%Im(t(t(t(spec_mat)*sqrt(1+Gamma*lambda)))%*%t(des_mat)))
  # scaler makes scales of regularization and unconstrained optimization `similar'
  scaler<-mean(diag(xtx))
  
  #  scaler<-min(10000,mean(diag(xtx)))
  
  
  if (sum(abs(weight_structure))>0)
  {  
    X_inv<-solve(xtx+xtx_structure+xtx_structure_diff+scaler*reg_mat)#dim(xtx_structure)  dim(xtx) 
    bh<-as.vector(X_inv%*%(((t(Re(X_new)*weight_target))%*%Gamma)+
                             ((t(Re(X_new_structure)*weight_target_structure))%*%Gamma_structure)+
                             ((t(Re(X_new_structure_diff)*weight_target_structure_diff))%*%Gamma_structure_diff)-
                             xtxy-scaler*reg_xtxy))
  } else
  {
    X_inv<-solve(xtx+scaler*reg_mat)#sum(diag(xtx))   sum(diag(xtx+scaler*reg_mat)) det(xtx+scaler*reg_mat)  xtx%*%solve(reg_mat) det(reg_mat)
    bh<-as.vector(X_inv%*%(((t(Re(X_new)*weight_target))%*%Gamma)-xtxy-scaler*reg_xtxy))    
  }
  # the last two filter weights are functions of the previous ones through the first and second order restrictions
  b<-matrix(nrow=L,ncol=length(weight_h_exp[1,]))
  # Reconstruct original parameters
  bhh<-t(des_mat)%*%bh
  for (k in 1:L) #k<-1
  {      #dim(t(des_mat))
    b[k,]<-bhh[(k)+(0:(length(weight_h_exp[1,])-1))*L]
  }
  
  
  # Modifications 17.04.2012 : the newly defined vector w_eight allows for simple/straightforward adjustment of filter coefficients
  weight_cm<-matrix(w_eight,ncol=(length(weight_h_exp[1,])))
  # Add level and/or time-shift constraints (if i1<-F and i2<-F then this matrix is zero)
  b<-b+weight_cm#*(sum(b))        #apply(b*(0:(L-1)),2,sum)/sum(b)
  
  
  # Transferfunction
  trffkt<-matrix(nrow=K+1,ncol=length(weight_h_exp[1,]))
  trffkth<-trffkt
  trffkt[1,]<-apply(b,2,sum)
  trffkth[1,]<-trffkt[1,]
  #  b<-scale(b,center=T,scale=F)
  
  for (j in 1:length(weight_h_exp[1,]))
  {
    for (k in 0:(K))#k<-1
    {
      trffkt[k+1,j]<-(b[,j]%*%exp(1.i*k*(0:(L-1))*pi/(K)))
    }
  }
  # Speed enhancer: if we need filter coefficients only we may skip all diagnostics computations below
  if (weights_only)
    return(list(b=b,trffkt=trffkt,rever=NULL,degrees_freedom=NULL,aic=NULL,freezed_degrees=NULL,Accuracy=NULL,Smoothness=NULL,Timeliness=NULL,
                MS_error=NULL))
  
  
  
  # The following derivations of the DFA-criterion are all equivalent
  # They are identical with rever (up to normalization by (2*(K+1)^2))  as calculated at the end of the function except if i1<-T and weight_constraint different from zero
  # In the latter case an additional constant interfers with the estimation
  
  # The target Y in the frequency-domain is the real vector weight_target*Gamma (both vectors are real)
  # The Regression estimate (the smoother) of Y is the following expression:
  # Note that we here ignore strucural elements in the target specification i.e. we treat structural elements as regularization elemenets
  trth<-((X_new)%*%(X_inv%*%t(Re(X_new))))%*%(weight_target*Gamma)                    #sum(abs(trth-trt))
  # This expression is identical to trt computed below if lambda=0 (assuming i1<-F or weight_constraint=0); otherwise trth is identical to Re(trt)+1.i*sqrt(1+lambda*Gamma)*Im(trt))
  # The projection matrix (it is a projection for the real part only, see below) is therefore:
  Proj_mat<-((X_new)%*%(X_inv%*%t(Re(X_new))))              #trth-(Proj_mat)%*%trth
  # The residual projection matrix is (it is a projection for the real part, see below)
  res_mat<-diag(rep(1,dim(Proj_mat)[1]))-Proj_mat
  # DFA criterion: first possibility (all three variants are identical)
  sum(abs(res_mat%*%(weight_target*Gamma))^2)
  # Residuals (DFT of filter error): in contrast to OLS this is not iid (weighted regression in frequency-domain)
  resi<-res_mat%*%(weight_target*Gamma)
  # DFA criterion: second possibility
  t(Conj(resi))%*%resi
  t((weight_target*Gamma))%*%(t(Conj(res_mat))%*%(res_mat))%*%(weight_target*Gamma)
  # The interesting `effective degrees of freedom' used here emphasizes an unbiased estimate of the mean-square residual
  #    See  http://en.wikipedia.org/wiki/Degrees_of_freedom_(statistics) (effective degrees of freedom: the expression tr((I-H)(I-H)')
  #    Note that res_mat=I-H see above
  #    Then (Y-Y^)(Y-Y^)/tr((I-H)(I-H)') is an unbiased estimate of the mean-square residual error (in our case Y=weight_target*Gamma. see above)
  #    This correcting `effective degrees of freedom' can then be used to implement a generalized AIC, see below
  degrees_freedom<-2*Re(sum(diag(t(Conj(res_mat))%*%(res_mat))))-1
  # Note that res_mat is a projection matrix with regards to the real part (but not with respect to the imaginary part)
  # Thus we could replace diag(t(Conj(res_mat))%*%(res_mat)) by diag(res_mat) in the degrees_freedom above i.e. the following expression equates to zero:
  Re(t(Conj(res_mat))%*%(res_mat))-Re(res_mat)
  freezed_degrees<-2*K+1-degrees_freedom
  # This is an alternative (identical) expression for the freezed_degreees
  2*Re(sum(diag(Proj_mat)))
  
  # DFA Criterion: third possibility (here an additional normalization by 1/(2*(K+1)^2))
  sum(abs(Gamma*weight_target-trth)^2)
  
  #ts.plot(b)
  trt<-apply(((trffkt)*exp(1.i*(0-Lag)*pi*(0:(K))/K))*weight_h_exp,1,sum)
  # DFA criterion which accounts for customization but not for regularization term
  # MDFA-Legacy : new normalization for all error terms below 
  rever<-sum(abs(Gamma*weight_target-Re(trt)-1.i*sqrt(1+lambda*Gamma)*Im(trt))^2)*2*pi/(K+1)
  # MS-filter error : DFA-criterion without effects by lambda or expweight (one must divide spectrum by expweight_vec)
  MS_error<-sum((abs(Gamma*weight_target-trt)/expweight_vec)^2)*2*pi/(K+1)
  # Definition of Accuracy, time-shift and noise suppression terms
  # Please note that:
  #       1. these terms refer to the original non-linearized criterion: therefore they do not sum up to rever
  #       2. we are interested in decomposing the mean-square error i.e. we ignore expweight and lambda here (we just want to measure the impact of lambda and expweight)
  #               Therefore, we use the DFT without expweight_vec i.e. we have to divide all spectral constents by expweight_vec
  Gamma_cp<-Gamma[1+0:as.integer(K*(cutoff/pi))]
  Gamma_cn<-Gamma[(2+as.integer(K*(cutoff/pi))):(K+1)]
  trt_cp<-(trt/expweight_vec)[1+0:as.integer(K*(cutoff/pi))]
  trt_cn<-(trt/expweight_vec)[(2+as.integer(K*(cutoff/pi))):(K+1)]
  weight_target_cp<-(weight_target/expweight_vec)[1+0:as.integer(K*(cutoff/pi))]
  weight_target_cn<-(weight_target/expweight_vec)[(2+as.integer(K*(cutoff/pi))):(K+1)]
  # define singular observations
  Accuracy<-sum(abs(Gamma_cp*weight_target_cp-abs(trt_cp))^2)*2*pi/(K+1)
  Timeliness<-4*sum(abs(Gamma_cp)*abs(trt_cp)*sin(Arg(trt_cp)/2)^2*weight_target_cp)*2*pi/(K+1)
  Smoothness<-sum(abs(Gamma_cn*weight_target_cn-abs(trt_cn))^2)*2*pi/(K+1)
  Shift_stopband<-4*sum(abs(Gamma_cn)*abs(trt_cn)*sin(Arg(trt_cn)/2)^2*weight_target_cn)*2*pi/(K+1)
  # Check: the following expression should vanish
  Accuracy+Timeliness+Smoothness+Shift_stopband-MS_error
  # Very prototypical: AIC
  #  aic<-ifelse(degrees_freedom<K+1&degrees_freedom>1,log(rever)+2*(K-degrees_freedom)/(K)+2*(K-degrees_freedom)*(K-degrees_freedom+1)/(K*(degrees_freedom-1)),NA)
  aic<-ifelse(degrees_freedom<K+1&degrees_freedom>1,log(rever)+2*(K-degrees_freedom+1)/(degrees_freedom-2),NA)
  
  return(list(b=b,trffkt=trffkt,rever=rever,degrees_freedom=degrees_freedom,aic=aic,freezed_degrees=freezed_degrees,Accuracy=Accuracy,Smoothness=Smoothness,Timeliness=Timeliness,MS_error=MS_error))
  
}