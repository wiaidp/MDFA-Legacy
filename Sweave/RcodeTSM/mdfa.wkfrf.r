mdfa.wkfrf <- function(delta.noise,delta.signal,spec.noise,spec.signal)
{
  ########################################
  ## mdfa.wkfrf
  #   computes model-based (WK) frequency response function
  #     for optimal signal extraction filter
  #   delta.noise: N x N x (d.noise+1) array of noise differencing
  #   delta.signal: N x N x (d.signal+1) array of signal differencing
  #		spec.noise: array of dimension N x N x (grid+1), consisting of 
  #     spectrum for differenced noise process 
  #     at frequencies pi*j/grid for 0 <= j <= grid
  #		spec.signal: array of dimension N x N x (grid+1), consisting of 
  #     spectrum for differenced signal process 
  #     at frequencies pi*j/grid for 0 <= j <= grid
  #	Outputs:
  #   frf.wk: array of dimension N x N x (grid+1) of frf, where
  #		  grid is the desired number of frequencies 
  #  Note: if no delta.noise or delta.signal portions, pass in identity matrices!
  #  Conventions: normal polynomial conventions for differencing polynomials
  #########################################
  
  d.noise <- dim(delta.noise)[3] - 1
  d.signal <- dim(delta.signal)[3] - 1
  N <- dim(spec.noise)[1]
  grid <- dim(spec.noise)[3]-1
  lambda <- pi*seq(0,grid)/grid
  
  frf.noise <- t(rep(1,(grid+1))) %x% delta.noise[,,1]  
  if(d.noise > 0) {
    for(i in 1:d.noise)
    {
      frf.noise <- frf.noise + t(exp(-1i*lambda*i)) %x% delta.noise[,,i+1]  
    } }
  frf.noise <- array(frf.noise,c(N,N,(grid+1)))
  
  frf.signal <- t(rep(1,(grid+1))) %x% delta.signal[,,1]  
  if(d.signal > 0) {
    for(i in 1:d.signal)
    {
      frf.signal <- frf.signal + t(exp(-1i*lambda*i)) %x% delta.signal[,,i+1]  
    } }
  frf.signal <- array(frf.signal,c(N,N,(grid+1)))
  
  spec.signal.del <- array(t(rep(1,(grid+1)) %x% diag(N)),c(N,N,(grid+1)))
  spec.noise.del <- array(t(rep(1,(grid+1)) %x% diag(N)),c(N,N,(grid+1)))
  spec.data.del <- array(t(rep(1,(grid+1)) %x% diag(N)),c(N,N,(grid+1)))
  frf.wk <- array(t(rep(1,(grid+1)) %x% diag(N)),c(N,N,(grid+1)))
  for(k in 1:(grid+1))
  {
    spec.noise.del[,,k] <- frf.signal[,,k] %*% spec.noise[,,k] %*% Conj(t(frf.signal[,,k]))
    spec.signal.del[,,k] <- frf.noise[,,k] %*% spec.signal[,,k] %*% Conj(t(frf.noise[,,k]))
    spec.data.del[,,k] <- spec.noise.del[,,k] + spec.signal.del[,,k]
    frf.wk[,,k] <- spec.signal.del[,,k] %*% solve(spec.data.del[,,k])
  }
  
  return(frf.wk)
}