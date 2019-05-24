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
  #  Requires: mvar.specfact.r, mdfa.frf.r, and mdfa.coeff.r
  #########################################
  
  d.noise <- dim(delta.noise)[3] - 1
  d.signal <- dim(delta.signal)[3] - 1
  N <- dim(spec.noise)[1]
  grid <- dim(spec.noise)[3]
  m <- floor(grid/2)
  lambda <- 2*pi*(seq(1,grid) - (m+1))/grid
  
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
  
  spec.signal.del <- array(t(rep(1,grid) %x% diag(N)),c(N,N,grid))
  spec.noise.del <- array(t(rep(1,grid) %x% diag(N)),c(N,N,grid))
  spec.data.del <- array(t(rep(1,grid) %x% diag(N)),c(N,N,grid))
  frf.wh <- array(t(rep(1,grid) %x% diag(N)),c(N,N,grid))
  for(k in 1:grid)
  {
    spec.noise.del[,,k] <- frf.signal[,,k] %*% spec.noise[,,k] %*% Conj(t(frf.signal[,,k]))
    spec.signal.del[,,k] <- frf.noise[,,k] %*% spec.signal[,,k] %*% Conj(t(frf.noise[,,k]))
    spec.data.del[,,k] <- spec.noise.del[,,k] + spec.signal.del[,,k]
  }
  
  acf.data.del <- mdfa.coeff(spec.data.del,0,len)
  data.sf <- mvar.specfact(aperm(acf.data.del,c(1,3,2)))
  data.theta <- data.sf[[1]][,,seq(len+1,1),drop=FALSE]
  data.sigma <- data.sf[[2]]
  
  f.ma <- t(rep(1,grid)) %x% data.theta[,,1]  
  if(len > 0) {
    for(i in 1:len)
    {
      f.ma <- f.ma + t(exp(-1i*lambda*i)) %x% data.theta[,,i+1]  
    } }
  f.ma <- array(f.ma,c(N,N,grid))
  
  g.spec <- array(t(rep(1,grid) %x% diag(N)),c(N,N,grid))
  for(k in 1:grid)
  {
    g.spec[,,k] <- spec.signal.del[,,k] %*% solve(Conj(t(f.ma[,,k])))
  }
  causalg.acf <- mdfa.coeff(g.spec,0,len)
  causalg.frf <- mdfa.frf(causalg.acf,0,grid)
  for(k in 1:grid)
  {
    frf.wh[,,k] <- causalg.frf[,,k] %*% solve(data.sigm) %*% solve(f.ma[,,k])
   }
 
  return(frf.wh)
}