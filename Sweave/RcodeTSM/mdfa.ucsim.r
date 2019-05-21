mdfa.ucsim <- function(delta,innovations)
{
  
  ###########################
  #	mdfa.ucsim
  #		by Tucker McElroy
  #
  #	generates (T-d) x N simulations from input innovations,
  #   which are T x N dimensional white noise.
  #  delta is the degree d differencing polynomial
  #		
  ##############################
  
  N <- dim(innovations)[2]
  T <- dim(innovations)[1]
  d <- length(delta) - 1
  delta.recurse <- -delta[-1]/delta[1]
  sims <- as.matrix(filter(innovations[(d+1):T,,drop=FALSE],
                           delta.recurse,method="recursive",
                           innovations[1:d,,drop=FALSE])[-seq(1,d),])
  
  return(sims)
}


