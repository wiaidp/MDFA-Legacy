
set.seed(1)
len<-120
eps1<-rnorm(len)
eps2<-rnorm(len)
# Define the data-matrix:
# The first column must be the target series. 
data_matrix<-cbind(eps1,eps1,eps2)
insample<-nrow(data_matrix)
# Compute the DFT: d=0 for stationary data (default settings)
weight_func<-spec_comp(insample, data_matrix, d)$weight_func
# Target
Gamma<-(1:nrow(weight_func))<=(nrow(weight_func)-1)/6+1
L<-12
# Source the default (MSE-) parameter settings
source(file=paste(path_MDFA.pgm,"control_default.r",sep=""))
# Estimate filter coefficients: MSE
mdfa_obj<-mdfa_analytic(K,L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,b0_H0,c_eta,weights_only=F,weight_structure,white_noise,synchronicity,lag_mat)





Proj_mat<-((X_new)%*%(X_inv%*%t(Re(X_new))))
res_mat<-diag(rep(1,dim(Proj_mat)[1]))-Proj_mat


H<-Proj_mat
res_mat<-diag(rep(1,dim(Proj_mat)[1]))-H
Y<-(weight_target*Gamma)
resi<-res_mat%*%Y
Yhat<-Proj_mat%*%Y


# projection for real part only
Re(H)-Re(t(Conj(H))%*%H)
Im(H)-Im(t(Conj(H))%*%H)

# Eigenvalues
Re(eigen(H)$values)
eigen(H%*%H)$values
eigen(t(Conj(H))%*%H)$values

eigen(res_mat)$values
# This is wrong because conjugation is omitted, but transposition is done
eigen((res_mat)%*%res_mat)$values
# This is correct
eigen(t(Conj(res_mat))%*%res_mat)$values

Im(H)+Im(Conj(t(H)))


sum(eigen(Re(res_mat))$values^2)
sum(eigen(Re(res_mat))$values)
sum(diag((res_mat)))
sum(abs(eigen(Im(res_mat))$values)^2)
sum(eigen(Conj(t(res_mat))%*%(res_mat))$values)
sum(diag(Conj(t(res_mat))%*%(res_mat)))


Re(Conj(t(H))%*%H)-Re(H)%*%Re(H)-Im(t(H))%*%Im(H)
Re(H)-t(Re(H))  
Re(H)%*%Re(H)-Re(H)  




abs(eigen(H)$values)


sum((eigen(t(Conj(H))%*%H)$values))

sum(eigen(H)$values)

# res_mat is not symmetric but Re(res_mat) is
res_mat-t(res_mat)
Re(res_mat)-t(Re(res_mat))



# Orthogonality applies to real part only
t(resi)%*%Conj(Yhat)


Im(H)+Im(t(Conj(H)))-Im((t(Conj(H)))%*%H)


2*K+1-2*Re(sum(diag(t(Conj(res_mat))%*%(res_mat))))
eigen(Conj(res_mat)%*%(res_mat))$values

sum(diag((t(Conj(H))%*%H)))

sum(diag(H))


sum(diag(Proj_mat))

2*sum(diag(t(Proj_mat)%*%Conj(Proj_mat)))-1

abs(diag((t(Conj(H))%*%H)))


#------------------------------------------------------------------------------------

set.seed(2)
len<-120
X<-cbind(rnorm(len),rnorm(len),rnorm(len))
Y<-cbind(rnorm(len),rnorm(len),rnorm(len))

Z<-X
Z<-X+1.i*Y

# 1. Unconstrained complex mean-square (not MDFA): inverse and end account for imaginary parts
M<-Z%*%solve(t(Conj(Z))%*%Z)%*%t(Conj(Z))
# 2. MDFA minus Imaginary part at start
M<-Re(Z)%*%solve(Re(t(Conj(Z))%*%Z))%*%t(Re(Z))
# Augmented hat matrix (MDFA plus imaginary part at the end): sums the correct number of degrees of freedom!
M<-Z%*%solve(Re(t(Conj(Z))%*%Z))%*%t(Conj(Z))
#M<-Re(Z%*%solve(Re(t(Conj(Z))%*%Z))%*%t(Conj(Z)))
# MDFA estimate
M<-Z%*%solve(Re(t(Conj(Z))%*%Z))%*%t(Re(Z))

# Interesting result when M is augmented hat-matrix: all eigenvalues vanish except first and last one.   L
eigen(Conj(t(M))-Conj(t(M))%*%M)$values
sum(eigen(Conj(t(M))-Conj(t(M))%*%M)$values)



diag(M)

#t(Conj(M))%*%M-M
#eigen(M)$values
sum(eigen(Re(M))$values)
sum(eigen(Im(M))$values)
sum(eigen(M)$values)

res<-diag(rep(1,dim(M)[1]))-M
diag(res)
eigen(res)$values
sum(eigen(Re(res))$values)
sum(eigen(Im(res))$values)
sum(eigen(res)$values)

sum(eigen(res%*%Conj(t(res)))$values)
sum(eigen(Conj(t(res))%*%res)$values)

res%*%Z



# Split real and imaginary parts

M<-Z%*%solve(Re(t(Conj(Z))%*%Z))%*%t(Re(Z))
N<-Z%*%solve(Re(t(Conj(Z))%*%Z))%*%t(Im(Z))

sum(eigen(Re(M))$values)
sum(eigen(Im(M))$values)

sum(eigen(Re(N))$values)
sum(eigen(Im(N))$values)

sum(eigen(Re(M))$values)+sum(eigen(Im(N))$values)
sum(eigen(Re(M)+Im(N))$values)
eigen(Re(M))$values+eigen(Im(N))$values


Re(solve(t(Conj(Z))%*%Z))-solve(Re(t(Conj(Z))%*%Z))



