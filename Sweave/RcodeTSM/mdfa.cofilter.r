mdfa.cofilter <- function(frf,spec,noisefreqs,beta,R,Q)
{
  
  #######################################################
  #
  #	mdfa.cofilter by Tucker McElroy
  #
  #	computes optimal concurrent moving average filter
  #		as approx of given target filter, using MDFA method,
  #		based upon the moving average filter class of length q
  #	inputs:
  #		frf is array N x N x Grid of complex entries, the target
  #			frequency response function Psi(e^{-i lambda})
  #			for lambda given by Grid number of Fourier frequencies 
  #		spec is array N x N+1 x Grid of complex entries, f(lambda):
  #     the process/data cross spectrum of coint process with data process,
  #			followed by process/data spectral density matrix 
  #		noisefreqs lists frequencies lambda=pi*x, with x in [-1,1], 
  #			with repeats for double roots, 
  #			such that exp(-i*lambda) is a root
  #			of the noise differencing polynomial
  #   beta is N vector of co-integrating relations
  #		R is array N x q x N x M of q-M constraints on the moving 
  #			average filter solution.
  #		Q is array N x q x N of constraints, theta = R phi + Q.
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
  
  spec.xx 
  
  fpsi <- NULL
  fmat <- NULL
  lambda.ft <- exp(-1i*2*pi*Grid^{-1}*(seq(1,Grid) - (m+1)))	## this is e^{-i lambda}
  
  opt.val <- do.call(cbind,lapply(seq(1,Grid),function(i) frf[,,i] %*% spec[,,i] %*% Conj(t(frf[,,i]))))
  opt.val <- Grid^{-1}*opt.val %*% (rep(1,Grid) %x% diag(N))
  for(k in 0:(q-1))
  {
    fpsi.new <- do.call(cbind,lapply(seq(1,Grid),function(i) frf[,,i] %*% spec[,,i]))
    fpsi.new <- Grid^{-1}*fpsi.new %*% (lambda.ft^{-k} %x% diag(N))
    fpsi <- cbind(fpsi,fpsi.new)
    fmat.new <- Grid^{-1}*matrix(spec,nrow=N) %*% (lambda.ft^{-k} %x% diag(N))
    if(k==0) { 
      fmat <- fmat.new 
      fzero <- fmat.new
    } else {
      if(k==1) {
        fmat <- cbind(fmat,fmat.new)
        fmat <- rbind(fmat,cbind(t(fmat.new),fzero))
      } else {
        side.mat <- fmat[1:(dim(fmat)[2]-N),(dim(fmat)[2]+1-N):dim(fmat)[2],drop=FALSE]
        fmat <- cbind(fmat,rbind(fmat.new,side.mat))
        fmat <- rbind(fmat,cbind(t(fmat.new),t(side.mat),fzero))
      }
    }
  }
  fpsi <- Re(fpsi)
  fmat <- Re(fmat)
  
  opt <- solve(t(R.mat) %*% fmat %*% R.mat) %*% t(R.mat) %*% (t(fpsi) - fmat %*% Q.mat)
  opt <- R.mat %*% opt + Q.mat
  opt.val <- Re(opt.val) - t(opt) %*% fmat %*% opt
  opt.array <- array(t(opt),c(N,N,q))
  
  return(list(opt.array,opt.val)) 
}



