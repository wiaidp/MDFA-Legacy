mdfa.hi2low <- function(filter.hi,hi.freq,low.freq,shift.hi)
{
  
  #######################################################
  #
  #	mdfa.hi2low by Tucker McElroy
  #
  #  embeds a given high frequency filter as a low frequency filter
  #     for a vector embedding
  #	inputs:
  #  filter.hi is given scalar high frequency filter of length L
  #  hi.freq is sampling frequency of high frequency time series
  #  low.freq is sampling frequency of low frequency time series
  #  shift.hi gives the integer offset for the high frequency filter:
  #     filter coefficients have indices -shift.hi,...,0,...,L-1-shift.hi
  #     set shift.hi = 0 for a causal filter
  #	outputs:
  #   filter.low is array s x s x M,
  #     where s = hi.freq/low.freq (must be an integer) 
  #     and M is length, which depends on s, L, and shift.hi
  #   shift.low is integer offset for the low frequency filter
  ##############################################################
  
  s.embed <- hi.freq/low.freq
  len.hi <- length(filter.hi)
  left.pad <- NULL
  if(shift.hi > 0) { left.pad <- rep(0,s.embed - (shift.hi %% s.embed)) }
  right.pad <- NULL
  if((len.hi-1-shift.hi) > 0) { right.pad <- 
    rep(0,2*s.embed - ((len.hi-1-shift.hi) %% s.embed)-1) }
  filter.embed <- NULL
  next.column <- c(left.pad,filter.hi,right.pad)
  len.low <- length(next.column)/s.embed
  for(k in 1:s.embed)
  {
    next.column <- c(0,rev(rev(next.column)[-1]))
    filter.embed <- cbind(filter.embed,next.column)
  }
  filter.low <- array(filter.embed,c(s.embed,len.low,s.embed))
  filter.low <- aperm(filter.low,c(1,3,2))
  shift.low <- (length(left.pad) + shift.hi)/s.embed - 1
  
  return(list(filter.low,shift.low))
}

  
  