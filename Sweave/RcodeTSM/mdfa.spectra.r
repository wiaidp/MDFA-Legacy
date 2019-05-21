mdfa.spectra <- function(phi,theta,sigma,grid)
{
  ########################################
  ## mdfa.spectra
  #   computes VARMA spectral density matrix
  #   phi: N x N x (p+1) array of VAR(p)
  #   theta: N x N x (q+1) array of VMA(q)
  #   sigma: N x N covariance matrix
  #		grid: desired number of frequencies for output spectrum
  #	Outputs:
  #		f.spec: array of dimension N x N x (grid+1), consisting of spectrum
  #			at frequencies pi*j/grid for 0 <= j <= grid
  #  Note: if no VAR or VMA portions, pass in identity matrices!
  #  Conventions: normal polynomial conventions for phi and theta!
  #########################################
  
  p <- dim(phi)[3] - 1
  q <- dim(theta)[3] - 1
  N <- dim(sigma)[1]
  lambda <- pi*seq(0,grid)/grid
  
  f.ma <- t(rep(1,(grid+1))) %x% theta[,,1]  
  if(q > 0) {
    for(i in 1:q)
    {
      f.ma <- f.ma + t(exp(-1i*lambda*i)) %x% theta[,,i+1]  
    } }
  f.ma <- array(f.ma,c(N,N,(grid+1)))
  
  f.ar <- t(rep(1,(grid+1))) %x% phi[,,1]  
  if(p > 0) {
    for(i in 1:p)
    {
      f.ar <- f.ar + t(exp(-1i*lambda*i)) %x% phi[,,i+1] 
    } 
  }
  f.ar <- array(f.ar,c(N,N,(grid+1)))
  
  f.wold <- array(t(rep(1,(grid+1)) %x% diag(N)),c(N,N,(grid+1)))
  f.spec <- f.wold
  for(k in 1:(grid+1))
  {
    f.wold[,,k] <- solve(f.ar[,,k]) %*% f.ma[,,k]
    f.spec[,,k] <- f.wold[,,k] %*% sigma %*% Conj(t(f.wold[,,k]))
  }
  
  return(f.spec)
}  
  
