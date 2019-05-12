\section{Customization: a Turning-Point Perspective}

In the previous sections we emphasized the connection between customization, Timeliness, Smoothness, and the scale-invariant Peak Correlation and Curvature statistics. We here widen the scope by linking customization to MSE-performances at the important turning-points of a time series. This section relies on \href{https://www.dropbox.com/s/69ihcnsyy6hqypd/DFA_beyond_max_lik.pdf?dl=0}{Wildi (2008)}
, chapter 5. The interested reader is referred to the cited literature for more comprehensive treatment.

\subsection{Definition}

For the ideal trend, turning-points can be uniquely identified as local extrema of the signal\footnote{If the trend signal is noisy i.e. if the target filter is leaking (Beveridge-Nelson trend for example), then turning-points cannot be defined in such a simple way.}, see fig.\ref{z_turning_point_sym}, top graph (we rely on the first AR(1)-process, $a_1=0.9$): vertical lines correspond to the local extrema of the trend with cutoff $\pi/12$. Alternatively, turning-points can be defined according to zero-crossings of the \emph{differenced} signal (sign-changes of first differences), see fig.\ref{z_turning_point_sym}, bottom graph. Both definitions are equivalent.

<<echo=False>>=
  # Specify the processes: ar(1) with coefficients -0.9,0 and 0.9
  a_vec<-0.9
# Ordinary ATS-components
scaled_ATS<-F
# Specify the lambdas
lambda_vec<-c(0,30)
# Specify the etas
eta_vec<-c(0,0)
# Specify filter length
L<-24
# Use periodogram
mba<-F
estim_MBA<-T
M<-len/2
# Length of symmetric target filter (for computing MSEs)
L_sym<-2*200
# Length of long data
len1<-2000
# Length of estimation sample
len<-120
# cutoff
cutoff<-pi/12
# Real-time design
Lag<-0
# no constraints
i1<-i2<-F
# difference data
dif<-F
# Single simulation run
anzsim<-1
@
  <<echo=True>>=
  # Proceed to simulation run
  for_sim_obj<-for_sim_out(a_vec,len1,len,cutoff,L,mba,estim_MBA,L_sym,Lag,i1,i2,scaled_ATS,lambda_vec,eta_vec,anzsim,M,dif)
@
  
<<echo=False>>=
  amp_shift_mat_sim<-for_sim_obj$amp_shift_mat_sim
amp_sim_per<-for_sim_obj$amp_sim_per
shift_sim_per<-for_sim_obj$shift_sim_per
xff_sim<-for_sim_obj$xff_sim
xff_sim_sym<-for_sim_obj$xff_sim_sym
ats_sym<-for_sim_obj$ats_sym
dim_names<-for_sim_obj$dim_names
x1<-for_sim_obj$x1
@
  <<echo=FALSE>>=
  file = paste("z_turning_point_sym.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
# extract all series from first realization (there is only one realization here)
anf<-300
enf<-1700
mplot_level<-cbind(x1,xff_sim_sym[,1,1,anzsim])[anf:enf,]#head(xf_per)
anf<-900
enf<-1100
par(mfrow=c(2,1))
mplot<-cbind(x1,xff_sim_sym[,1,1,anzsim])[anf:enf,]#head(xf_per)
plot(as.ts(mplot[,1]),type="l",axes=F,col="blue",ylim=c(min(na.exclude(mplot)),
                                                        max(na.exclude(mplot))),ylab="",xlab="",
     main=paste("Ideal trend and turning-points",sep=""),lwd=1)
lines(mplot[,2])
axis(1,at=c(1,rep(0,6))+as.integer((0:6)*(enf-anf)/6),
     labels=as.integer(c(1,rep(0,6))+(0:6)*(enf-anf)/6))
abline(v=1+which(diff(mplot[1:(nrow(mplot)-1),2])*diff(mplot[2:nrow(mplot),2]-1)<0))
axis(2)
box()
mplot_diff<-as.matrix(diff(xff_sim_sym[,1,1,1])[-1+anf:enf])#head(xf_per)
plot(as.ts(mplot_diff[,1]),type="l",axes=F,col="black",ylim=c(min(na.exclude(mplot_diff)),
                                                              max(na.exclude(mplot_diff))),ylab="",xlab="",
     main=paste("Zero-crossings of the differenced signal ",sep=""),lwd=1)
axis(1,at=c(1,rep(0,6))+as.integer((0:6)*(enf-anf)/6),
     labels=as.integer(c(1,rep(0,6))+(0:6)*(enf-anf)/6))
abline(v=1+which(diff(mplot[1:(nrow(mplot)-1),2])*diff(mplot[2:nrow(mplot),2]-1)<0))
abline(h=0)
axis(2)
box()
dev.off()
@
  <<label=z_turning_point_sym.pdf,echo=FALSE,results=tex>>=
  file = paste("z_turning_point_sym", sep = "")
cat("\\begin{figure}[H]")
cat("\\begin{center}")
cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
cat("\\caption{Ideal trend and turning points (a1=0.9)", sep = "")
cat("\\label{z_turning_point_sym}}", sep = "")
cat("\\end{center}")
cat("\\end{figure}")
@
  
  
  
  
  
  \subsection{Customization: Efficiency Gains in Turning-Points}


For illustration we here rely on zero-crossings of the differenced signal (bottom graph in previous figure). Accordingly, we compute real-time filters for the \emph{differenced} data and we compare the ordinary MSE-filter $\lambda=\eta=0$ with the customized design $\lambda=30,\eta=0$. To be clear, the (quadratic) criterion \ref{idfatp} used here is
\begin{eqnarray}
&&\frac{2\pi}{T} \sum_{k=-[T/2]}^{[T/2]}
\left|\Gamma(\omega_k)-\hat{\Gamma}(\omega_k)\right|^2 I_{T\Delta X}(\omega_k)\nonumber\\
&+&\lambda\frac{2\pi}{T} \sum_{k=-[T/2]}^{[T/2]}
A(\omega_k)\hat{A}(\omega_k)^2\sin\left(\hat{\Phi}(\omega_k)-\Phi(\omega_k)\right)^2I_{T\Delta X}(\omega_k)\to\min_{\mathbf{b}}\label{idfatp_diff}
\end{eqnarray}
where $I_{T\Delta X}(\omega_k)$ is the periodogram of the \emph{differenced} data, $\lambda=(0,30)$, and $\eta=0$, for both filters; the latter implies $W(\omega_k,\eta)\simeq 1$. We now proceed to estimation
<<echo=False>>=
  # Difference the data
dif<-T
anzsim<-10
for_sim_obj<-for_sim_out(a_vec,len1,len,cutoff,L,mba,estim_MBA,L_sym,Lag,i1,i2,scaled_ATS,lambda_vec,eta_vec,anzsim,M,dif)
@
  and plot the filter-errors $y_t-\hat{y}_t^{MSE}$ and $y_t-\hat{y}_t^{Cust}$, see fig.\ref{z_turning_point_sym_diff}. 

<<echo=False>>=
amp_shift_mat_sim<-for_sim_obj$amp_shift_mat_sim
amp_sim_per<-for_sim_obj$amp_sim_per
shift_sim_per<-for_sim_obj$shift_sim_per
xff_sim<-for_sim_obj$xff_sim
xff_sim_sym<-for_sim_obj$xff_sim_sym
ats_sym<-for_sim_obj$ats_sym
dim_names<-for_sim_obj$dim_names
x1<-for_sim_obj$x1
@
  
  <<echo=FALSE>>=
  file = paste("z_turning_point_sym_diff.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
# extract all series from first realization (there is only one realization here)
anf<-300
enf<-1700
index_var<-!is.na(xff_sim_sym[anf:enf,1,1,anzsim])
scale_sym<-sqrt(var((xff_sim_sym[anf:enf,1,1,anzsim])[index_var]))
scale_xf_2<-sqrt(var((xff_sim[anf:enf,2,1,anzsim])[index_var]))
scale_xf_3<-sqrt(var((xff_sim[anf:enf,3,1,anzsim])[index_var]))

#var((xff_sim[anf:enf,2,1,anzsim])[index_var]/scale_xf_2)


mplot_err_long<-cbind(xff_sim_sym[,1,1,anzsim],xff_sim_sym[,1,1,anzsim]-xff_sim[,2:3,1,anzsim])[anf:enf,]#head(xf_per)
# we should scale because of zero-shrinkage of customized design
mplot_err_long<-cbind(xff_sim_sym[,1,1,anzsim],
                      xff_sim_sym[,1,1,anzsim]-xff_sim[,2,1,anzsim],
                      xff_sim_sym[,1,1,anzsim]-xff_sim[,3,1,anzsim]*scale_xf_2/scale_xf_3)[anf:enf,]#head(xf_per)


tp_loc_long<-which(mplot_err_long[1:(nrow(mplot_err_long)-1),1]*mplot_err_long[2:nrow(mplot_err_long),1]<0)
anf<-900
enf<-1100
mplot_err<-cbind(xff_sim_sym[,1,1,anzsim],xff_sim_sym[,1,1,anzsim]-xff_sim[,2:3,1,anzsim])[anf:enf,]#head(xf_per)


index_var<-!is.na(xff_sim_sym[anf:enf,1,1,anzsim])
scale_sym<-sqrt(var((xff_sim_sym[anf:enf,1,1,anzsim])[index_var]))
scale_xf_2<-sqrt(var((xff_sim[anf:enf,2,1,anzsim])[index_var]))
scale_xf_3<-sqrt(var((xff_sim[anf:enf,3,1,anzsim])[index_var]))


mplot_err<-cbind(xff_sim_sym[,1,1,anzsim]/scale_sym,
                 xff_sim_sym[,1,1,anzsim]/scale_xf_2-xff_sim[,2,1,anzsim]/scale_xf_2,
                 xff_sim_sym[,1,1,anzsim]/scale_xf_3-xff_sim[,3,1,anzsim]/scale_xf_3)[anf:enf,]#head(xf_per)
tp_loc<-which(mplot_err[1:(nrow(mplot_err)-1),1]*mplot_err[2:nrow(mplot_err),1]<0)
plot(as.ts(mplot_err[,2]),type="l",axes=F,col="red",ylim=c(min(na.exclude(mplot_err)),
                                                           max(na.exclude(mplot_err))),ylab="",xlab="",
     main=paste("Filter error",sep=""),lwd=1)
mtext(paste("MSE-diff (",lambda_vec[1],",",eta_vec[1],")",sep=""), side = 3, line = -1,at=(enf-anf)/2,col="red")
lines(mplot_err[,3],col="green")
mtext(paste("Customized-diff (",lambda_vec[2],",",eta_vec[2],")",sep=""), side = 3, line = -2,at=(enf-anf)/2,col="green")
axis(1,at=c(1,rep(0,6))+as.integer((0:6)*(enf-anf)/6),
     labels=as.integer(c(1,rep(0,6))+(0:6)*(enf-anf)/6))
#tp_loc<-which(diff(mplot_level[1:(nrow(mplot_level)-1),2])*diff(mplot_level[2:nrow(mplot_level),2]-1)<0)
abline(v=tp_loc)
abline(h=0)
axis(2)
box()
dev.off()
@
  <<label=z_turning_point_sym_diff.pdf,echo=FALSE,results=tex>>=
  file = paste("z_turning_point_sym_diff", sep = "")
cat("\\begin{figure}[H]")
cat("\\begin{center}")
cat("\\includegraphics[height=3in, width=6in]{", file, "}\n",sep = "")
cat("\\caption{Filter error MSE-diff (red) vs. customized-diff (green), a1=0.9", sep = "")
cat("\\label{z_turning_point_sym_diff}}", sep = "")
cat("\\end{center}")
cat("\\end{figure}")
@
  We observe that the customized filter (green) seems to outperform the MSE-design in the turning-points (vertical lines): the filter-error of the former is closer to zero, in the mean. In contrast, the MSE-design seems to outperform the customized design between the turning-points. In order to verify our conjecture we compute the corresponding MSEs, see table \ref{MSE_tp_nontp_short}.  
<<echo=True>>=
  mse_mat<-rbind(apply(mplot_err[,2:3]^2,2,mean),apply(mplot_err[tp_loc,2:3]^2,2,mean))
dimnames(mse_mat)[[2]]<-c("MSE-diff","Cust-diff")
dimnames(mse_mat)[[1]]<-c("MSE all time points","MSE at turning-points")
@
  <<label=print_tab_cor,echo=FALSE,results=tex>>=
  library(Hmisc)
require(xtable)
#latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
#center = "centering", file = "", floating = FALSE)
xtable(mse_mat, dec = 1,digits=6, caption = paste("MSE-performances of MSE-diff (left) and customized-diff (right) in all time points (top) and in the turning-points only (bottom)",sep=""),label=paste("MSE_tp_nontp_short",sep=""),
       center = "centering", file = "", floating = FALSE)
@
  Since the turning-points are scarce, we use a much longer series (length 2000: approximately 100 TPs) and we recompute MSEs, see table \ref{MSE_tp_nontp_long}: note that these figures represent out-of-sample performances (the estimation span is of length $T=120$, only). 
<<echo=True>>=
  mse_mat_long<-rbind(apply(mplot_err_long[,2:3]^2,2,mean),apply(mplot_err_long[tp_loc_long,2:3]^2,2,mean))
dimnames(mse_mat_long)[[2]]<-c("MSE-diff","Cust-diff")
dimnames(mse_mat_long)[[1]]<-c("MSE all time points","MSE at turning-points")
@
  <<label=print_tab_cor,echo=FALSE,results=tex>>=
  library(Hmisc)
require(xtable)
#latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
#center = "centering", file = "", floating = FALSE)
xtable(mse_mat_long, dec = 1,digits=6, caption = paste("MSE-performances of MSE-diff (left) and customized-diff (right) in all time points (top) and in the turning-points only (bottom): very long time span",sep=""),label=paste("MSE_tp_nontp_long",sep=""),
       center = "centering", file = "", floating = FALSE)
@
  The interesting effect is clearly confirmed by a comparison of left and right columns in the last table: the customized design (right) outperforms the MSE-design at the turning-points (TPs), at costs of degraded performances between TPs. A top-down comparison is evocative too: the performance of the customized design improves at the TPs. We can rely on simple geometric arguments for decoding these visual hints: 
  \begin{enumerate}
\item The slope of the differenced signal or, equivalently, the curvature of the original signal, are strongest in the TPs.
\item To a first-order approximation, the error imputable to the time-shift is proportional to the slope of the signal: a shift of $\delta$ time-units, by the real-time filter, generates a filter-error of $s\delta$, where $s$ is the instantaneous slope of the target signal.
\end{enumerate}
By combining these two arguments we infer that the customized filter  improves MSE-performances specifically in the TPs, where the  Timeliness effect on the MSE, being proportional to the slope, is the most impactful\footnote{The effect of the inflated Accuracy component is negligible, in the TPs, because the signal is vanishing: amplitude mismatches in the passband are negligible.}. Interestingly, this localized selective efficiency-gain is obtained without any a priori knowledge about the true location of the TPs: it applies `mechanically' or, more precisely, `geometrically' (whenever the curvature of the signal in the TPs is large\footnote{The latter property is typical for many economic time series, because TPs are often disruptive, unexpected events, whose `curvature dynamics' can be magnified by asymmetric risk-aversion (fear psychology) and/or by adverse (pro-cyclical) regulatory market-interventions.}). We infer that customization emphasizes MSE-performances in those important time points which are most relevant for economic analysis, namely the turning points\footnote{The interested reader is referred to \href{https://www.dropbox.com/s/69ihcnsyy6hqypd/DFA_beyond_max_lik.pdf?dl=0}{Wildi (2008)}, section 5.2, for a comparative analysis of logit-models and customized designs.}.




\subsection{A Link to Timeliness and Smoothness}

We can re-write  criterion \ref{idfatp_diff} by relying on the discrete convolution result \ref{conv_per}:
  \begin{eqnarray}
&&\frac{2\pi}{T} \sum_{k=-[T/2]}^{[T/2]}
\left|\Gamma(\omega_k)-\hat{\Gamma}(\omega_k)\right|^2 I_{T\Delta X}(\omega_k)\nonumber\\
&+&\lambda\frac{2\pi}{T} \sum_{k=-[T/2]}^{[T/2]}
A(\omega_k)\hat{A}(\omega_k)^2\sin\left(\hat{\Phi}(\omega_k)-\Phi(\omega_k)\right)^2I_{T\Delta X}(\omega_k)\nonumber\\
&\approx&\frac{2\pi}{T} \sum_{k=-[T/2]}^{[T/2]}
\left|\Gamma(\omega_k)-\hat{\Gamma}(\omega_k)\right|^2 |1-\exp(-i\omega_k)|^2I_{TX}(\omega_k)\nonumber\\
&+&\lambda\frac{2\pi}{T} \sum_{k=-[T/2]}^{[T/2]}
A(\omega_k)\hat{A}(\omega_k)^2\sin\left(\hat{\Phi}(\omega_k)-\Phi(\omega_k)\right)^2|1-\exp(-i\omega_k)|^2I_{T X}(\omega_k)\label{idfatp_diff_2}  
\end{eqnarray}
where we used the fact that 
\[I_{T\Delta X}(\omega_k)\approx |1-\exp(-i\omega_k)|^2I_{T X}(\omega_k)\]
see \ref{conv_per}. We infer that the resulting customized design can be interpreted in two different ways, according to which input-series is fed to the filter: 
  \begin{itemize}
\item If the differenced data is used, as we did above, then the filter-design emphasizes Timeliness only ($\lambda=30$, $\eta=0$). MSE-performances improve specifically in the important turning-points of the data (zero-crossings of the differenced signal).
\item If the original (undifferenced) data is used, then \ref{idfatp_diff_2}  suggests that customization emphasizes Timeliness ($\lambda=30$) \emph{and} Smoothness. Indeed, the difference operator $|1-\exp(-i\omega_k)|^2$, in \ref{idfatp_diff_2} can be assimilited to a weighting-function  $W(\omega_k,\eta)$, in \ref{idfatp}\footnote{There are minor (negligible) differences in the passband, where $W(\omega_k,\eta)$ is supposed to be `neutral'.}. 
\end{itemize}
If we interpret the filter-design in terms of original (unfiltered) data, then the proposed customized design emphasizes Timeliness \emph{and} Smoothness much like the balanced design in sections \ref{double_score_ats} and \ref{ucdvbmseli}.   
