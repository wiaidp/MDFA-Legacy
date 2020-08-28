mvar.filter <- function(x,filter)
{
  
  #######################################################
  #
  #	mvar.filter by Tucker McElroy
  #
  #	computes multivariate filter output
  #	inputs:
  #   x is T x N multivariate time series 
  #   filter is array 
  #	outputs:
  #		opt.array is array N x N x q of filter coefficients
  #		opt.val is N x N matrix corresponding to minimal MSE
  #
  ##############################################################
  
  N <- dim(spec)[1]
  Grid <- dim(frf)[3]
  m <- floor(Grid/2)
  q <- dim(R)[2]
  M <- dim(R)[4]
  R.mat <- matrix(R,nrow=N*q,ncol=N*M)
  Q.mat <- matrix(Q,ncol=N)
  
  fpsi <- NULL
  fmat <- NULL
  lambda.ft <- exp(-1i*2*pi*Grid^{-1}*(seq(1,Grid) - (m+1)))	## this is e^{-i lambda}
  
  opt.val <- do.call(cbind,lapply(seq(1,Grid),function(i) frf[,,i] %*% spec[,,i] %*% Conj(t(frf[,,i]))))
  opt.val <- Grid^{-1}*opt.val %*% (rep(1,Grid) %x% diag(N))
  for(k in 0:(q-1))
    