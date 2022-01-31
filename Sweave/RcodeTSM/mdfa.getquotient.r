mdfa.getquotient <- function(psivals,rootfreqs)
{
  
  #######################################################
  #
  #	mdfa.getquotient by Tucker McElroy
  #
  #	computes quotient polynomial after having divided by Delta(z).
  #
  # Background: consider Psi(z) = Delta(z)Psi^{sharp}(z) + Psi^{flat}(z),
  #   where Delta(z) has degree d and Psi^{flat}(z) has degree d-1.
  #   Suppose there are r roots zeta_j with multiplicity m_j.
  #   Then Psi^(k) (zeta_j) = {Psi^{flat}}^k (zeta_j) for 0 <= k < m_j.
  #   We compute Psi^{sharp}(z) at the roots of Delta(z)
  #	inputs:
  #   psivals: vector of d complex numbers, corresponding to Psi^(k) (zeta_j),
  #     ordered by first root and all multiplicities, and then the next root...
  #		rootfreqs: lists frequencies x in [-1,1] such that lambda=pi*x, 
  #			with repeats for double roots, 
  #			such that exp(-i*lambda) is a root of Delta(z).
  #     Presumes bot x and -x are included root is non-real.
  #	outputs:
  #   remvals: vector of d real numbers, corresponding to the coefficients
  #     of Psi^{flat} (z), starting with constant coefficient.
  #
  #######################################################
  
  root.lambdas <- unique(rootfreqs)
  root.mults <- NULL
  if(length(rootfreqs)>0) 
  {
    for(j in 1:length(root.lambdas))
    {
      root.mults <- c(root.mults,sum(rootfreqs == root.lambdas[j]))
    } 
  }
  
  d <- length(psivals)
  remvals <- mdfa.getquotient(psivals,rootfreqs)
  for(j in 1:length(root.lambdas))
  {
    zeta <- exp(1i*pi*root.lambdas[j])
    k <- root.mults[j]
    weights <- c(rep(0,k),factorial(k)*choose(seq(k,d-1),k))
    print(weights*zeta^(seq(0,d-1)-k))
    
  }
  
  
  
}  