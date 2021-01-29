mdfa.whfrf <- function(delta.noise,delta.signal,spec.noise,spec.signal,len)
{
  ########################################
  ## mdfa.whfrf
  #   computes model-based  frequency response function
  #     for optimal concurrent (Wiener-Hopf) signal extraction filter
  #   delta.noise: N x N x (d.noise+1) array of noise differencing
  #   delta.signal: N x N x (d.signal+1) array of signal differencing
  #		spec.noise: array of dimension N x N x grid, consisting of 
  #     spectrum for differenced noise process at Fourier frequencies
  #		spec.signal: array of dimension N x N x grid, consisting of 
  #     spectrum for differenced signal process at Fourier frequencies
  #   len: parameter giving degree of truncation in calculations
  #	Outputs:
  #   frf.wh: array of dimension N x N x grid of frf, where
  #		  grid is the desired number of frequencies 
  #  Note: if no delta.noise or delta.signal portions, pass in identity matrices!
  #  Conventions: normal polynomial conventions for differencing polynomials
  #  Requires: polymulMat.r, mvar.specfact.r, mdfa.frf.r, and mdfa.coeff.r
  #########################################
  
  array_combine <- function(a_array,b_array,inv=FALSE)
  {
    N <- dim(a_array)[1]
    len <- dim(a_array)[3]
    c_array <- array(0,c(N,N,len))
    if(inv) 
    {
      for(i in 1:len) { c_array[,,i] <- a_array[,,i] %*% solve(b_array[,,i]) }
    } else 
    {
      for(i in 1:len) { c_array[,,i] <- a_array[,,i] %*% b_array[,,i] }
    }
    return(c_array)
  }
  
  d.noise <- dim(delta.noise)[3] - 1
  d.signal <- dim(delta.signal)[3] - 1
  N <- dim(spec.noise)[1]
  grid <- dim(spec.noise)[3]
  m <- floor(grid/2)
  lambda <- 2*pi*(seq(1,grid) - (m+1))/grid
  delta.process <- polymulMat(delta.noise,delta.signal)
  d.process <- dim(delta.process)[3] - 1
  
  frf.noise <- t(rep(1,grid)) %x% delta.noise[,,1]  
  if(d.noise > 0) {
    for(i in 1:d.noise)
    {
      frf.noise <- frf.noise + t(exp(-1i*lambda*i)) %x% delta.noise[,,i+1]  
    } }
  frf.noise <- array(frf.noise,c(N,N,grid))
  frf.signal <- t(rep(1,grid)) %x% delta.signal[,,1]  
  if(d.signal > 0) {
    for(i in 1:d.signal)
    {
      frf.signal <- frf.signal + t(exp(-1i*lambda*i)) %x% delta.signal[,,i+1]  
    } }
  frf.signal <- array(frf.signal,c(N,N,grid))
  frf.process <- array_combine(frf.noise,frf.signal,FALSE)
  
  spec.signal.del <- array(t(rep(1,grid)) %x% diag(N),c(N,N,grid))
  spec.noise.del <- array(t(rep(1,grid)) %x% diag(N),c(N,N,grid))
  spec.data.del <- array(t(rep(1,grid)) %x% diag(N),c(N,N,grid))
  frf.wk <- array(t(rep(1,grid)) %x% diag(N),c(N,N,grid))
  for(k in 1:grid)
  {
    spec.noise.del[,,k] <- frf.signal[,,k] %*% spec.noise[,,k] %*% Conj(t(frf.signal[,,k]))
    spec.signal.del[,,k] <- frf.noise[,,k] %*% spec.signal[,,k] %*% Conj(t(frf.noise[,,k]))
    spec.data.del[,,k] <- spec.noise.del[,,k] + spec.signal.del[,,k]
    frf.wk[,,k] <- spec.signal.del[,,k] %*% solve(spec.data.del[,,k])
  }
  
  filter.wk <- mdfa.coeff(frf.wk,len,len)
  acf.data.del <- mdfa.coeff(spec.data.del,0,len)
  data.sf <- mvar.specfact(aperm(acf.data.del,c(1,3,2)))
  data.theta <- data.sf[[1]][,,seq(len+1,1),drop=FALSE]
  data.sigma <- data.sf[[2]]
  frf.theta <- t(rep(1,grid)) %x% data.theta[,,1]
  for(i in 1:len)
  {
    frf.theta <- frf.theta + t(exp(-1i*lambda*i)) %x% data.theta[,,i+1]  
  }
  
  frf.wh <- t(rep(1,grid)) %x% filter.wk[,,len+1]
  phi <- array(0,c(N,N,len+1))
  phi[,,1] <- data.theta[,,1] %*% solve(delta.process[,,1])
  phi_partial <- t(rep(1,grid)) %x% phi[,,1]
  for(i in 1:len)
  {
    phi[,,i+1] <- data.theta[,,i+1] %*% solve(delta.process[,,1])
    for(j in 1:min(i,d.process))
    {
      phi[,,i+1] <- phi[,,i+1] - phi[,,i-j+1] %*% 
        delta.process[,,j+1] %*% solve(delta.process[,,1])
    }
    causal_frf <- array_combine(phi_partial,frf.theta,TRUE)   
    causal_frf <- array_combine(causal_frf,frf.process,FALSE)
    causal_frf <- t(rep(1,grid)) %x% diag(N) - causal_frf
    causal_frf <- array_combine(t(exp(1i*lambda*i)) %x% filter.wk[,,len+1-i],causal_frf)
    frf.wh <- frf.wh + t(exp(-1i*lambda*i)) %x% filter.wk[,,len+1+i]
    frf.wh <- frf.wh + causal_frf
    phi_partial <- phi_partial + t(exp(-1i*lambda*i)) %x% phi[,,i+1]
  }
  frf.wh <- array(frf.wh,c(N,N,grid))
  
  return(frf.wh)
}