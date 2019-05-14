# MDFA
q <- 20	# length of moving average filter
r <- T/2	# number of autocovariances used in MDFA
psidfa.array <- array(0,c(2,2,(2*r+q-2)))
psidfa.array[,,(r+q+1)] <- diag(2)
fore.mdfa <- mdfa.maunconstrained(psidfa.array,x.acf[,,1:r],q)

#plot(ts(fore.mdfa[[1]][1,1,]),xlab="Index",ylab="")
#plot(ts(fore.mdfa[[1]][1,2,]),xlab="Index",ylab="")



# MDFA
q <- 20	# length of moving average filter
r <- T/2	# number of autocovariances used in MDFA
psidfa.array <- array(0,c(2,2,(2*r+q-2)))
psidfa.array[,,(r+q+1)] <- diag(2)
fore.mdfa <- mdfa.maunconstrained(psidfa.array,x.acf[,,1:r],q)
 
## MDFA
q <- 12	# length of moving average filter
r <- T/2
dfa.lp <- lp.filter[(r+q-2+len+1):(1-r+len+1)]
psidfa.array <- array(t(dfa.lp %x% diag(2)),c(2,2,length(dfa.lp)))
fore.mdfa <- mdfa.maunconstrained(psidfa.array,x.acf[,,1:r],q)

# get MDFA concurrent filter
q <- 12	# length of moving average filter
r <- T/2
dfa.lp <- lp.filter[(r+q-2+len+1):(1-r+len+1)]
psidfa.array <- array(t(dfa.lp %x% diag(2)),c(2,2,length(dfa.lp)))
fore.mdfa <- mdfa.maunconstrained(psidfa.array,x.acf[,,1:r],q)
