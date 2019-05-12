

\subsubsection{Low-Pass}

The (mean) cycle-length in our model is \Sexpr{round(2*pi/Arg(polyroot(c(1,-fit2_07$par[1:2])))[1],2)} quarters. Isolating the trend from the cycle is almost impossible, in a real-time setting, because the spectral peaks are confounded. Instead of relying on a more or less unstable model-based trend-signal\footnote{Whose signature strongly depends on whether the great recession is in-sample or not.}, obtained \emph{ex cathedra} by optimizing short-term forecast performances, we suggest that the DFA is more amenable to reconcile research priorities and optimization principle. For that purpose we interface with the optimization criterion by specifying the interesting signal: an ideal trend with cutoff $\pi/32$ (16 years). Note that the peak of the AR(2)-cycle lies in the stopband of the specified trend target. 
\begin{enumerate}
\item Specify the filter design: cutoff, restriction ($i1==T$), MSE-design ($\lambda=\eta=0$)
<<echo=True>>=
  cutoff_len<-32
cutoff<-pi/cutoff_len
L<-100
# Spectrum: MDFA requires DFT i.e. square-root (imaginary part is irrelevant since design is univariate)
weight_func<-cbind(sqrt(trffkt_GDP),sqrt(trffkt_GDP))
# remove singularity in frequency zero (we impose a first order constraint)
weight_func[1,]<-0
# Target
Gamma<-(0:K)<=as.integer(cutoff*K/pi)+1
# Restrictions: i1 constraint
i1<-T
i2<-F
weight_constraint<-Gamma[1]
# MSE-design
lambda<-eta<-0
@
<<echo=False>>=
  # Additional configuration settings for MDFA
  d<-0
weight_constraint<-Gamma[1]
lambda_cross<-lambda_smooth<-0
lambda_decay<-c(0,0)
lin_expweight<-F
shift_constraint<-rep(0,ncol(weight_func)-1)
grand_mean<-F
b0_H0<-NULL
c_eta<-F
weights_only<-F
weight_structure<-c(0,0)
white_noise<-F
synchronicity<-F
lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)-1),nrow=L)  
@
\item Proceed to estimation: we select $L=100$\footnote{Mainly for convenience: the filter is not subject to overfitting since we know the true DGP and since the frequency-grid is fine enough ($K=$\Sexpr{K}).}.
<<echo=True>>=
  # Estimate MDFA MSE filter coefficients  
  mdfa_obj<-mdfa_analytic(K,L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,  lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,b0_H0,c_eta,weights_only=F,weight_structure,white_noise,synchronicity,lag_mat)
@      
\item Plot amplitude and time-shift functions, see fig.\ref{z_us_real_log_gdp_detrended_amp_shift}.
<<echo=False>>=
  file = paste("z_us_real_log_gdp_detrended_amp_shift.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)

mplot<-cbind(Gamma,abs(mdfa_obj$trffkt),Arg(mdfa_obj$trffkt)/omega_k)#head(mplot)
# Complete by shift in frequency zero
mplot[1,3]<-as.double(t(mdfa_obj$b)%*%(0:(L-1)))/sum(mdfa_obj$b)
dimnames(mplot)[[2]]<-c("Target","Estimate","Shift")
ax<-rep(NA,nrow(mplot))
ax[1+(0:6)*((nrow(mplot)-1)/6)]<-c(0,"pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi")
par(mfrow=c(1,2))
plot_title<-"Amplitude"
insamp<-1.e+90
title_more<-dimnames(mplot)[[2]][1:2]
colo<-c("black","blue")
mplot_func(as.matrix(mplot[,1:2]), ax, plot_title, title_more, insamp, colo)
plot_title<-"Shift"
insamp<-1.e+90
title_more<-NULL
colo<-"blue"
mplot_func(as.matrix(mplot[,3]), ax, plot_title, title_more, insamp, colo)

dev.off()
@
<<label=z_us_real_log_gdp_detrended_amp_shift.pdf,echo=FALSE,results=tex>>=
  file = paste("z_us_real_log_gdp_detrended_amp_shift", sep = "")
cat("\\begin{figure}[H]")
cat("\\begin{center}")
cat("\\includegraphics[height=3in, width=6in]{", file, "}\n",sep = "")
cat("\\caption{Amplitude (left) and shift (right) of real-time lowpass MSE-design, i1=T", sep = "")
cat("\\label{z_us_real_log_gdp_detrended_amp_shift}}", sep = "")
cat("\\end{center}")
cat("\\end{figure}")
@
\item Filter the data, see fig.\ref{z_us_real_log_gdp_detrended_filt}: note that we centered the series.
<<echo=True>>=
  xf<-rep(NA,length(detrended))
for (i in L:length(detrended))  
  xf[i]<-t(mdfa_obj$b)%*%(detrended[i:(i-L+1)]-mean(detrended))  
@
<<echo=False>>=
  file = paste("z_us_real_log_gdp_detrended_filt.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
mplot<-scale(cbind(detrended,xf),center=T,scale=F)
dimnames(mplot)[[2]]<-c("Detrended US GDP","Low-pass filtered")
plot(mplot,ylim=c(min(mplot,na.rm=T),max(mplot,na.rm=T)),xlim=c(start_year,end_year),plot.type="s",col=c("black","blue"),main="Detrended Log Real US GDP (black) and lowpass MSE (blue)")
nberShade()
lines(mplot[,2],col="blue")
lines(mplot[,1])
mtext("detrended GDP", side = 3, line = -1,at=mean(c(start_year,end_year)),col="black")
mtext("Lowpass MSE", side = 3, line = -2,at=mean(c(start_year,end_year)),col="blue")
abline(h=0)
dev.off()
@
<<label=z_us_real_log_gdp_detrended_filt.pdf,echo=FALSE,results=tex>>=
  file = paste("z_us_real_log_gdp_detrended_filt", sep = "")
cat("\\begin{figure}[H]")
cat("\\begin{center}")
cat("\\includegraphics[height=3in, width=6in]{", file, "}\n",sep = "")
cat("\\caption{Detrende (log-real) US-GDP and low-pass MSE output   ", sep = "")
cat("\\label{z_us_real_log_gdp_detrended_filt}}", sep = "")
cat("\\end{center}")
cat("\\end{figure}")
@
\end{enumerate}
