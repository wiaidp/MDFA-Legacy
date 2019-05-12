
#insamp<-min(estimation_length,nrow(isData))
#x<-isData

# Copyright: Marc Wildi
# 15.01.2012
# http://blog.zhaw.ch/idp/sefblog
#
# new 2012-code: computes spectral estimates based on DFT
spec_comp <- function(insamp, x, d) {
  # non-stationary  
  if(d == 1) {
    weight_func <- periodogram_bp(diff(x[1 : insamp, 1]), 1, 
                                  insamp - 1)$fourtrans
    # explaining variables
    if(length(weight_func) > 1) {
      for(j in 2 : ncol(x)) {
        # since the data is integrated one uses the pseudo-periodogram: 
        # diff(data) and d <- 1
        per <- periodogram_bp(diff(x[1 : insamp, j]), 1, insamp - 1)$fourtrans
        weight_func <- cbind(weight_func, per)
      }
    }
  } else 
  {
#target    
    weight_func <- general_spec(x[1 : insamp, 1], 0, insamp)$fourtrans
# explaining variables
    if(length(weight_func) > 1) {
      for(j in 2 : ncol(x)) #j<-3
      {
        per <- general_spec(x[1 : insamp, j], 0, insamp)$fourtrans
        weight_func <- cbind(weight_func, per)
      }
    }
  }
  colnames(weight_func) <- colnames(x)
  # return the spectral estimate
  return(list(weight_func = weight_func))
}


# This function computes the DFT in the case of general mixed_frequency data
#   1. It accounts for NA's
#   2. It folds the DFT of the low-frequency data

#y<-x[,1]
general_spec<-function(y, d, insamp)
{
# Remove NA's
  z<-y[!is.na(y)]
# we could compute the lag until the first observation is available and rotate the DFT accordingly
# But a better/cleaner solution is to postpone these relative lags of the data to the computations 
# of the MDFA spectral matrix
  which(y==y[!is.na(y)][1])[1]-1
# compute the DFT of z  

  dft<-periodogram_bp(z,d, length(z))$fourtrans

# Extending the high-frequency part
#   Assumption: white noise  

  if (length(z)<length(y))
  {
    amp<-mean(abs(dft))
# the shift could be rough towards ferquency zero: therefore we prefer the median over the mean
    shift<-median(Arg(dft)[2:length(dft)]/(pi*(1:(length(dft)-1))/(length(dft)-1)))
# ts.plot(Arg(dft)[2:length(dft)]/(pi*(1:(length(dft)-1))/(length(dft)-1)))
# ts.plot(Arg(dft)[2:length(dft)])
    freq_ext<-(length(dft)):(length(y)/2)
    fourtrans<-c(dft,amp*exp(1.i*shift*pi*freq_ext/(length(y)/2)))
#  ts.plot(Arg(fourtrans))
  } else
  {
    fourtrans<-dft
  }
  return(list(fourtrans=fourtrans))
}






# DFT (old code but still in use for new 2012-version)  n.pg<-insamp dd<-0  x<-x[,3]
periodogram_bp <- function(x, dd, n.pg) {
  # preparations
  n.fit  <- length(x)
  xx     <- as.vector(x[((n.fit - n.pg + 1) : n.fit)])
  npg2   <- n.pg / 2
  perall <- fourtrans <- 0 * 0 : npg2
  
  # case without a seasonal component
  if (dd < 3) {
    for (j in 0 : npg2) {
      fourtrans[j + 1] <- xx %*% exp((1 : (2* npg2)) * 1.i * j * pi / npg2) #length(xx)
      fourtrans[j + 1] <- fourtrans[j + 1] / sqrt(2*pi * (2*npg2-1))
      perall[j + 1] <- abs(fourtrans[j + 1])^2
    }
  }
  
  # case with a seasonal component, special treatment for pi/6
  if (dd >= 3) {
    for (j in (1 : npg2)[(-npg2 / 6) * (1 : 6)]) {
      fourtrans[j + 1] <- xx %*% exp((1 : (2 * npg2)) * 1.i * j * pi / npg2)
      term2 <- abs(1 - exp(j * 1.i * pi / npg2))^2
      term3 <- abs(1 - exp(12 * j * 1.i * pi / npg2))^2
      perall[j + 1] <- abs(fourtrans[j + 1]) / (term2 * term3)
    }
    perall[(npg2 / 6) * (1 : 6) + 1] <- max(perall) * 100000
  }
  
  # output
  return(list(perall = perall, fourtrans = fourtrans))
}
