mdfa.getconstraints <- function(frf,sigfreqs,noisefreqs,q)
{

	#######################################################
	#
	#	mdfa.getconstraints by Tucker McElroy
	#
	#	computes filter constraints needed for nonstationary series,
	#		based on sigfreqs and noisefreqs.
	#	inputs:
	#		sigfreqs lists frequencies lambda=pi*x, with x in [-1,1], 
	#			with repeats for double roots, 
	#			such that exp(-i*lambda) is a root
	#			of the signal differencing polynomial
	#		noisefreqs lists frequencies lambda=ppi*x, with x in [-1,1], 
	#			with repeats for double roots, 
	#			such that exp(-i*lambda) is a root
	#			of the noise differencing polynomial
	#		frf is array N x N x Grid of complex entries, the target
	#			frequency response function Psi(e^{-i lambda})
	#			for lambda given by Grid number of Fourier frequencies 
  	#           q: desired order of MA filter
	#	outputs:
	#		R is array N x q x N x M of q-M constraints on the moving 
	#			average filter solution.
	#		Q is array N x q x N of constraints, theta = R phi + Q.
	#
	##############################################################

	N <- dim(frf)[1]
	Grid <- dim(frf)[3]
	m <- floor(Grid/2)

	sig.lambdas <- unique(sigfreqs)
	sig.mults <- NULL
	for(j in 1:length(sig.lambdas))
	{
		sig.mults <- c(sig.mults,sum(sigfreqs == sig.lambdas[j]))
	}

	noise.lambdas <- unique(noisefreqs)
	noise.mults <- NULL
	for(j in 1:length(noise.lambdas))
	{
		noise.mults <- c(noise.mults,sum(noisefreqs == noise.lambdas[j]))
	}

	sig.mat <- NULL
	sig.vec <- NULL
	if(length(sig.lambdas) > 0) {
	for(k in 1:length(sig.lambdas))
	{
		j.star <- floor(sig.lambdas[k]*Grid/2) + m+1
		if(j.star==(Grid+1)) j.star <- 1
		zeta <- exp(-1i*sig.lambdas[k]*pi)
		sig.mat.new <- zeta^(seq(0,q-1))
		sig.vec.new <- frf[,,j.star]
		sig.mat <- rbind(sig.mat,Re(sig.mat.new))
		sig.vec <- rbind(sig.vec,Re(sig.vec.new))
		if( (sig.lambdas[k] != 0) && (sig.lambdas[k] != 1) )
#		if(Im(zeta) != 0) 
		{ 
			sig.mat <- rbind(sig.mat,Im(sig.mat.new))
			sig.vec <- rbind(sig.vec,Im(sig.vec.new))
		}
		if(sig.mults[k]==2)
		{
			sig.mat.new <- seq(0,q-1)*zeta^(seq(-1,q-2))
			if(j.star==Grid) 
			{ 
				sig.vec.new <- 1i*exp(1i*sig.lambdas[k]*pi)*
					(frf[,,1]-frf[,,j.star])*Grid/(2*pi) 
			} else 
			{ 
				sig.vec.new <- 1i*exp(1i*sig.lambdas[k]*pi)*
					(frf[,,j.star+1]-frf[,,j.star])*Grid/(2*pi) 
			}
			sig.mat <- rbind(sig.mat,Re(sig.mat.new))
			sig.vec <- rbind(sig.vec,Re(sig.vec.new))
			if( (sig.lambdas[k] != 0) && (sig.lambdas[k] != 1) )
#			if(Im(zeta) != 0) 
			{ 
				sig.mat <- rbind(sig.mat,Im(sig.mat.new))
				sig.vec <- rbind(sig.vec,Im(sig.vec.new))
			}
		}
 	} }
	
	noise.mat <- NULL
	noise.vec <- NULL
	if(length(noise.lambdas) > 0) {
	for(k in 1:length(noise.lambdas))  
	{
		j.star <- floor(noise.lambdas[k]*Grid/2) + m+1
		if(j.star==(Grid+1)) j.star <- 1
		zeta <- exp(-1i*noise.lambdas[k]*pi)
		noise.mat.new <- zeta^(seq(0,q-1))
		noise.vec.new <- frf[,,j.star]
		noise.mat <- rbind(noise.mat,Re(noise.mat.new))
		noise.vec <- rbind(noise.vec,Re(noise.vec.new))
		if( (noise.lambdas[k] != 0) && (noise.lambdas[k] != 1) )
#		if(Im(zeta) != 0) 
		{ 
			noise.mat <- rbind(noise.mat,Im(noise.mat.new))
			noise.vec <- rbind(noise.vec,Im(noise.vec.new))
		}
		if(noise.mults[k]==2)
		{
			noise.mat.new <- seq(0,q-1)*zeta^(seq(-1,q-2))
			if(j.star==Grid) 
			{ 
				noise.vec.new <- 1i*exp(1i*noise.lambdas[k]*pi)*
					(frf[,,1]-frf[,,j.star])*Grid/(2*pi) 
			} else 
			{ 
				noise.vec.new <- 1i*exp(1i*noise.lambdas[k]*pi)*
					(frf[,,j.star+1]-frf[,,j.star])*Grid/(2*pi) 
			}	
			noise.mat <- rbind(noise.mat,Re(noise.mat.new))
			noise.vec <- rbind(noise.vec,Re(noise.vec.new))
			if( (noise.lambdas[k] != 0) && (noise.lambdas[k] != 1) )
#			if(Im(zeta) != 0) 
			{ 
				noise.mat <- rbind(noise.mat,Im(noise.mat.new))
				noise.vec <- rbind(noise.vec,Im(noise.vec.new))
			}
		}
 	} }
 
	constraint.mat <- rbind(sig.mat,noise.mat)
	constraint.vec <- rbind(sig.vec,noise.vec)

	## compute decomposition
	constraint.qr <- qr(constraint.mat)
	constraint.q <- qr.Q(constraint.qr)
	constraint.r <- qr.R(constraint.qr)
	constraint.pivot <- constraint.qr$pivot
	constraint.ipivot <- sort.list(constraint.pivot)

	M <- q - dim(constraint.r)[2] + dim(constraint.q)[2]
	R.mat <- rbind(-solve(constraint.r[,1:(dim(constraint.q)[2]),drop=FALSE],
		constraint.r[,(dim(constraint.q)[2]+1):q,drop=FALSE]),diag(q-M))
	R.mat <- R.mat[constraint.ipivot,] %x% diag(N)
	Q.mat <- rbind(solve(constraint.r[,1:(dim(constraint.q)[2]),drop=FALSE]) %*% 
		solve(constraint.q),matrix(0,q-M,M)) 
	Q.mat <- (Q.mat[constraint.ipivot,] %x% diag(N)) %*% constraint.vec
	R <- array(R.mat,c(N,q,N,q-M))
	Q <- array(Q.mat,c(N,q,N))
 
	return(list(R,Q))
}



