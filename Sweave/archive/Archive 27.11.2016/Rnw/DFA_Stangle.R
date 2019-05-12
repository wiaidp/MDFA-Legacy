### R code from vignette source ''
### Encoding: ISO8859-1

###################################################
### code chunk number 1: init
###################################################
library(tseries)
library(fGarch)

# Set paths
path.main<-"C:\\wia_desktop\\Projekte\\2013\\Unterricht\\Eco2\\"
path.pgm<-paste(path.main,"R\\",sep="")
path.out<-paste(path.main,"Latex\\",sep="")
path.dat<-paste(path.main,"Data\\",sep="")

#-----------------------------------------------------------------
# Source estimation routine
#source(paste(path.pgm,"examples.r",sep=""))


###################################################
### code chunk number 2: eco2.rnw:154-171
###################################################
cutoff<-pi/6
K<-1200
Gamma_h<-((0:K)*pi/K)<cutoff
Gamma<-rep(c(Gamma_h[length(Gamma_h):2],Gamma_h),3)
file = paste("z_Gamma_t.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
plot(Gamma,col="blue",main="Periodic Gamma",xlab="",ylab="",type="l",axes="F")
axis(1,at=c(c(1,0,0,0)+(length(Gamma)/3)*(0:(length(Gamma)/(length(Gamma)/3))),
length(Gamma)/2),labels=c("-3pi","-pi","pi","3pi","0"))
axis(2)
box()
plot(Gamma_h,col="blue",main="Gamma in [0,pi]",xlab="",ylab="",type="l",axes="F")
axis(1,at=c(1,(K/6)*1:6),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 3: z_Gamma_t.pdf
###################################################
  file = paste("z_Gamma_t.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Ideal lowpass with cutoff pi/6", sep = "")
  cat("\\label{z_Gamma_t}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 4: eco2.rnw:235-259
###################################################
cutoff<-pi/6
# Order of approximation : 10, 100, 1000, 10000
ord<-K
# Compute coefficients gamma
gamma<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))
sum(gamma)+sum(gamma[2:ord])
# Plot
file = paste("z_gamma.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,1))
plot(c(gamma[min(1000,ord):2],gamma[1:min(1000,ord)]),col="black",
main="Fourier coefficients of Gamma: -1000:1000",xlab="",ylab="",
type="l",axes="F")
axis(1,at=c(0:4)*ord/2,labels=c(-ord,-ord/2,0,ord/2,ord))
axis(2)
box()
plot(c(gamma[min(1000,ord/10):2],gamma[1:min(1000,ord/10)]),col="black",
main="Fourier coefficients of Gamma: -100:100",xlab="",ylab="",
type="l",axes="F")
axis(1,at=c(0:4)*ord/20,labels=c(-ord/10,-ord/20,0,ord/20,ord/10))
axis(2)
box()
dev.off()


###################################################
### code chunk number 5: z_gamma.pdf
###################################################
  file = paste("z_gamma.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Fourier coefficients of Gamma: from -1000 to 1000 (top) and -100 to 100 (bottom)", sep = "")
  cat("\\label{z_gamma}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 6: eco2.rnw:283-308
###################################################
cutoff<-pi/6
# Order of approximation : 10, 100, 1000, 10000
ord<-10000
# Compute coefficients gamma
gamma<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))
# Compute finite sum
# len1: Number of frequency ordinates (resolution of discrete grid in [-pi,pi])
len1<-K
Gamma_hat<-0:len1
for (k in 0:len1)#k<-0
{
	omegak<-k*pi/len1
	Gamma_hat[k+1]<-gamma%*%(cos(omegak*0:ord)*c(1,rep(2,ord)))
}
file = paste("z_Gamma_a.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
plot(Gamma_hat,type="l",axes=F,
main="Gamma (blue) and finite approximation (black)",xlab="",ylab="")
lines(Gamma_h,col="blue")
axis(1, at=c(0,1+1:6*len1/6),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 7: z_Gamma_a.pdf
###################################################
  file = paste("z_Gamma_a.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Finite approximation m=1000 (black) vs. Gamma (blue)", sep = "")
  cat("\\label{z_Gamma_a}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 8: eco2.rnw:325-326
###################################################
Gamma_hat[1+round(len1/6)]


###################################################
### code chunk number 9: eco2.rnw:334-354
###################################################
len<-300
set.seed(10)
eps<-rnorm(len)
y<-rep(NA,len)
for (k in 100:200)#k<-0
{
	y[k]<-gamma[1:100]%*%eps[k+(0:99)]+gamma[2:100]%*%eps[k-(1:99)]
}
file = paste("z_Gamma_n.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
ts.plot(as.ts(eps),type="l",axes=F,
main="Noise (blue) and lowpassed noise(red)",xlab="",
ylab="",col="blue")
lines(y,col="red",lwd=2)
axis(1, at=c(0,1+1:6*len1/6),labels=c("0","pi/6","2pi/6",
"3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 10: z_Gamma_n.pdf
###################################################
  file = paste("z_Gamma_n.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Noise (blue) and lowpassed noise (red)", sep = "")
  cat("\\label{z_Gamma_n}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 11: eco2.rnw:388-401
###################################################
cutoff_1<-pi/6
cutoff_2<-2*pi/6
K<-1200
Gamma_pb<-((0:K)*pi/K)>cutoff_1&((0:K)*pi/K)<cutoff_2
file = paste("z_Gamma_pb.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
plot(Gamma_pb,col="blue",main="Gamma passband in [0,pi]",xlab="",
ylab="",type="l",axes="F")
axis(1,at=c(1,(K/6)*1:6),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 12: z_Gamma_pb.pdf
###################################################
  file = paste("z_Gamma_pb.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Ideal passband with cutoffs pi/60 and pi/6", sep = "")
  cat("\\label{z_Gamma_pb}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 13: eco2.rnw:424-461
###################################################
len1<-1200
f<-cos(2*pi*(1:len1)/len1)
file = paste("z_cos_pb.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,1))
resolution<-10
plot(f,col="black",
main=paste("Approximation of Cosine by bp-functions: resolution=",resolution,sep=""),
xlab="",ylab="",type="l",axes="F")
for (i in 0:resolution)
  lines(cos(2*pi*(i+0.5)/resolution)*((1:len1)>(i)*len1/resolution&(1:len1)<
  (i+1)*len1/resolution),col="blue")
resolution<-33
axis(1,at=c(1,(K/6)*1:6),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
resolution<-33
plot(f,col="black",
main=paste("Approximation of Cosine by bp-functions: resolution=",resolution,sep=""),
xlab="",ylab="",type="l",axes="F")
for (i in 0:resolution)
  lines(cos(2*pi*(i+0.5)/resolution)*((1:len1)>(i)*len1/resolution&(1:len1)<
  (i+1)*len1/resolution),col="red")
axis(1,at=c(1,(K/6)*1:6),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
resolution<-100
plot(f,col="black",
main=paste("Approximation of Cosine by bp-functions: resolution=",resolution,sep=""),
xlab="",ylab="",type="l",axes="F")
for (i in 0:resolution)
  lines(cos(2*pi*(i+0.5)/resolution)*((1:len1)>(i)*len1/resolution&(1:len1)<
  (i+1)*len1/resolution),col="green")
axis(1,at=c(1,(K/6)*1:6),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 14: z_cos_pb.pdf
###################################################
  file = paste("z_cos_pb.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Approximation of cosine by bandpass functions", sep = "")
  cat("\\label{z_cos_pb}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 15: eco2.rnw:492-515
###################################################
# Compute coefficients gamma of the two lowpass filters
gamma_1<-c(cutoff_1/pi,(1/pi)*sin(cutoff_1*1:ord)/(1:ord))
gamma_2<-c(cutoff_2/pi,(1/pi)*sin(cutoff_2*1:ord)/(1:ord))
# Compute finite sum (approximation)
Gamma_hat_pb<-0:len1
for (k in 0:len1)#k<-0
{
	omegak<-k*pi/len1
	Gamma_hat_pb[k+1]<-(gamma_2-gamma_1)%*%
  (cos(omegak*0:ord)*c(1,rep(2,ord)))
}
file = paste("z_Gammapb_a.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
plot(Gamma_hat_pb,type="l",axes=F,
main="Passband Gamma (blue) and finite approximation (black)",xlab="",
ylab="")
lines(Gamma_pb,col="blue")
axis(1, at=c(0,1+1:6*len1/6),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 16: z_Gammapb_a.pdf
###################################################
  file = paste("z_Gammapb_a.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Finite approximation m=1000 (black) vs. passband Gamma (blue)", sep = "")
  cat("\\label{z_Gammapb_a}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")



#--------------------------------
# Section on DFT/IDFT


###################################################
### code chunk number 17: eco2.rnw:591-623
###################################################

# Exercises 2.2.1 Section 2.2

len<-301
set.seed(10)
# Seasonal AR(1)
model<-list(order=c(12,0,0),ar=c(rep(0,11),0.9))
y<-arima.sim(model,n=len)
# DFT
DFT<-0:(len/2)
DFT_c<-0:(len/2)
for (k in 0:(len/2))
{
	cexp <- complex(arg=-(1:len)*2*pi*k/len)
# Complex conjugate: this will be used for computing the IDFT
	cexpt <- complex(arg=(1:len)*2*pi*k/len)
	four<-sum(cexp*y*sqrt(1/(2*pi*len)))
# Complex conjugate: this will be used for computing the IDFT
	fourc<-sum(cexpt*y*sqrt(1/(2*pi*len)))
	DFT[k+1]<-four
# Complex conjugate: this will be used for computing the IDFT
	DFT_c[k+1]<-fourc
}
file = paste("z_sar1.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(3,1))
ts.plot(y,main="Data: SAR(1)")
acf(y)
plot(abs(DFT),type="l",axes=F,col="blue",main="DFT")
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 18: z_sar1.pdf
###################################################
  file = paste("z_sar1.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Seasonal process: data (top), acf (middle) and DFT (bottom)", sep = "")
  cat("\\label{z_sar1}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 19: eco2.rnw:640-660
###################################################
z<-1:len
zh<-z
for (k in 1:len)#k<-1
{
	cexp <- complex(arg=(0:(len/2))*2*pi*k/len)
# Complex canjugate
	cexpt <- complex(arg=-(0:(len/2))*2*pi*k/len)
# Account for weights w_k
	if (abs(len/2-as.integer(len/2))<0.1)
	{
		cexp<-cexp*c(rep(1,length(cexp)-1),(1/2))
		cexpt<-cexpt*c(rep(1,length(cexp)-1),(1/2))
	}
# IDFT
	z[k]<-sum(cexp[1:length(cexp)]*DFT[1:length(DFT)]*sqrt((2*pi)/len))+
  sum(cexpt[2:length(cexpt)]*DFT_c[2:length(DFT_c)]*sqrt((2*pi)/len))
  zh[k]<-sum(c(cexp[1]*DFT[1]*sqrt((2*pi)/len),2*Re(cexp[2:length(cexp)]*
  DFT[2:length(DFT)]*sqrt((2*pi)/len))))
}
cbind(y,z,zh)[(len-30):len,]


###################################################
### code chunk number 20: eco2.rnw:671-711
###################################################
# MA(1)
model_1<-list(order=c(0,0,1),ma=-0.9)
model_3<-list(order=c(0,0,1),ma=0.9)
set.seed(10)
y_1<-arima.sim(model_1,n=len)
set.seed(10)
y_2<-rnorm(len)
set.seed(10)
y_3<-arima.sim(model_3,n=len)
# DFT-----------
DFT_1<-0:(len/2)
DFT_2<-0:(len/2)
DFT_3<-0:(len/2)
for (k in 0:(len/2))
{
	cexp <- complex(arg=-(1:len)*2*pi*k/len)
	DFT_1[k+1]<-sum(cexp*y_1*sqrt(1/(2*pi*len)))
	DFT_2[k+1]<-sum(cexp*y_2*sqrt(1/(2*pi*len)))
	DFT_3[k+1]<-sum(cexp*y_3*sqrt(1/(2*pi*len)))
}
file = paste("z_ma1_dft.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(3,1))
plot(abs(DFT_1),type="l",axes=F,col="blue",main="|DFT|: b_1=-0.9",xlab="")
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
plot(abs(DFT_2),type="l",axes=F,col="blue",main="|DFT|: b_1=0",xlab="")
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
plot(abs(DFT_3),type="l",axes=F,col="blue",main="|DFT|: b_1=0.9",xlab="")
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 21: z_ma1_dft.pdf
###################################################
  file = paste("z_ma1_dft.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{MA(1): positive (top), zero (middle) and negative (bottom) MA-coefficients", sep = "")
  cat("\\label{z_ma1_dft}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 22: eco2.rnw:729-789
###################################################
# MA(1)
model_1<-list(order=c(0,0,1),ma=-0.9)
model_3<-list(order=c(0,0,1),ma=0.9)
# AR(2)
pers=0.99;
phi=pi/2;
a1=2*pers*cos(phi);
a2= -pers^2;
model<-list(order=c(2,0,0),ar=c(a1,a2))
simanz<-100
DFT<-matrix(nrow=len/2+1,ncol=3)
abs_dft<-matrix(rep(0,3*(1+(len-1)/2)),nrow=len/2+1,ncol=3)
# Simulation runs
for (i in 1:simanz)
{
  set.seed(i)
  y_1<-arima.sim(model_1,n=len)
  set.seed(i)
  y_2<-rnorm(len)
  set.seed(i)
  y_3<-arima.sim(model_3,n=len)
  # DFT-----------
  for (k in 0:(len/2))
  {
  	cexp <- complex(arg=-(1:len)*2*pi*k/len)
  	DFT[k+1,1]<-sum(cexp*y_1*sqrt(1/(2*pi*len)))
  	DFT[k+1,2]<-sum(cexp*y_2*sqrt(1/(2*pi*len)))
  	DFT[k+1,3]<-sum(cexp*y_3*sqrt(1/(2*pi*len)))
  }
  abs_dft<-abs_dft+abs(DFT)
}
abs_dft<-abs_dft/simanz
  
file = paste("z_ma1s_dft.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(3,1))
ymin<-0
ymax<-max(abs_dft[,1])
plot(abs_dft[,1],type="l",axes=F,col="blue",main="|DFT|: b_1=-0.9",
xlab="",ylim=c(ymin,ymax))
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6",
"5pi/6","pi"))
axis(2)
box()
ymax<-max(abs_dft[,2])
plot(abs_dft[,2],type="l",axes=F,col="blue",main="|DFT|: b_1=0",xlab="",
ylim=c(ymin,ymax))
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6",
"5pi/6","pi"))
axis(2)
box()
ymax<-max(abs_dft[,3])
plot(abs_dft[,3],type="l",axes=F,col="blue",main="|DFT|: b_1=0.9",
xlab="",ylim=c(ymin,ymax))
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6",
"5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 23: z_ma1s_dft.pdf
###################################################
  file = paste("z_ma1s_dft.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{MA(1): positive (top), zero (middle) and negative (bottom) MA-coefficients", sep = "")
  cat("\\label{z_ma1s_dft}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 24: eco2.rnw:808-811
###################################################
abs(polyroot(c(0.9801,-0.99,1)))
Arg(polyroot(c(0.9801,-0.99,1)))/pi
2*pi/Arg(polyroot(c(0.9801,-0.99,1)))


###################################################
### code chunk number 25: eco2.rnw:815-842
###################################################
# AR(2)
pers=0.99;
phi=pi/3;
a1=2*pers*cos(phi)
a2= -pers^2
model<-list(order=c(2,0,0),ar=c(a1,a2))
set.seed(10)
y_ar2<-arima.sim(model,n=len)
# DFT-----------
DFT_ar2<-0:(len/2)
for (k in 0:(len/2))
{
	cexp <- complex(arg=-(1:len)*2*pi*k/len)
	DFT_ar2[k+1]<-sum(cexp*y_ar2*sqrt(1/(2*pi*len)))
}
file = paste("z_ar2_dft.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,1))
ts.plot(y,xlab="",main="Data: AR(2)")
plot(abs(DFT_ar2),type="l",axes=F,col="blue",
main="|DFT| cyclical AR(2)",xlab="")
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6",
"5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 26: z_ar2_dft.pdf
###################################################
  file = paste("z_ar2_dft.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{AR(2): data (top) and abs(DFT) (bottom)", sep = "")
  cat("\\label{z_ar2_dft}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 27: eco2.rnw:861-889
###################################################
len<-1000
z_1<-1:len
z_2<-z
DFT_1<-rep(0,1+len/2)
DFT_1[10]<-1
DFT_2<-DFT_1
DFT_2[20]<-2
for (k in 1:len)#k<-1
{
	cexp <- complex(arg=(0:(len/2))*2*pi*k/len)
# Weights w_k
	if (abs(len/2-as.integer(len/2))<0.1)
	{
		cexp<-cexp*c(rep(1,length(cexp)-1),(1/2))
	}
# IDFT
  z_1[k]<-sum(c(cexp[1]*DFT_1[1]*sqrt((2*pi)/len),
  2*Re(cexp[2:length(cexp)]*DFT_1[2:length(DFT_1)]*sqrt((2*pi)/len))))
  z_2[k]<-sum(c(cexp[1]*DFT_2[1]*sqrt((2*pi)/len),
  2*Re(cexp[2:length(cexp)]*DFT_2[2:length(DFT_2)]*sqrt((2*pi)/len))))
}

file = paste("z_ls.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
ts.plot(z_1,xlab="",main="IDFT of single line spectrum")
ts.plot(z_2,xlab="",main="IDFT of double line spectrum")

dev.off()


###################################################
### code chunk number 28: z_ls.pdf
###################################################
  file = paste("z_ls.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{IDFT applied to single (top) and double (bottom) line spectra", sep = "")
  cat("\\label{z_ls}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 29: eco2.rnw:909-946
###################################################
set.seed(10)
len<-300
y<-cumsum(rnorm(len))
# DFT-----------
DFT<-0:(len/2)
for (k in 0:(len/2))
{
	cexp <- complex(arg=-(1:len)*2*pi*k/len)
	DFT[k+1]<-sum(cexp*y*sqrt(1/(2*pi*len)))
}
# IDFT
z<-y
for (k in 1:len)#k<-1
{
	cexp <- complex(arg=(0:(len/2))*2*pi*k/len)
# Weights w_k
	if (abs(len/2-as.integer(len/2))<0.1)
	{
		cexp<-cexp*c(rep(1,length(cexp)-1),(1/2))
	}
# IDFT
  z[k]<-sum(c(cexp[1]*DFT[1]*sqrt((2*pi)/len),2*Re(cexp[2:length(cexp)]*
  DFT[2:length(DFT)]*sqrt((2*pi)/len))))
}
cbind(y,z)[(len-20):len,]
file = paste("z_ns.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,1))
ts.plot(y,xlab="",main="Non-stationary data")
plot(abs(DFT),type="l",axes=F,col="blue",
main="|DFT| non-stationary data",xlab="")
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 30: z_ns.pdf
###################################################
  file = paste("z_ns.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Non-stationary data (top) and absolute DFT (bottom)", sep = "")
  cat("\\label{z_ns}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 31: eco2.rnw:962-984
###################################################
# IDFT
z<-rep(y,3)
for (k in len+1:(2*len))#k<-1
{
	cexp <- complex(arg=(0:(len/2))*2*pi*k/len)
# Weights w_k
	if (abs(len/2-as.integer(len/2))<0.1)
	{
		cexp<-cexp*c(rep(1,length(cexp)-1),(1/2))
	}
# IDFT
  z[k]<-sum(c(cexp[1]*DFT[1]*sqrt((2*pi)/len),
  2*Re(cexp[2:length(cexp)]*DFT[2:length(DFT)]*sqrt((2*pi)/len))))
}
z<-Re(z)
file = paste("z_nsf.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
ts.plot(c(z[1:len],rep(NA,2*len)),xlab="",
main="Non-stationary data (black) and forecasts based on IDFT (blue)")
lines(c(rep(NA,len),z[len+1:(2*len)]),col="blue",lwd=1)

dev.off()


###################################################
### code chunk number 32: z_nsf.pdf
###################################################
  file = paste("z_nsf.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Non-stationary data (black) and IDFT forecasts", sep = "")
  cat("\\label{z_nsf}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")



#------------------------------------------------
# Section on Periodogram

###################################################
### code chunk number 33: eco2.rnw:1079-1094
###################################################
len<-length(y_ar2)
Rk<-1:len
# Compute sample autocovariances
for (j in 1:len)
{
  Rk[j]<-y_ar2[j:len]%*%y_ar2[1:(len-j+1)]/len
}
# Compute periodogram based on sample autocovariances
per<-rep(0,1+len/2)
for (k in 0:as.integer(len/2))  #k<-1
{
  omega_k<-k*pi*2/len
  per[k+1]<-(1/(2*pi))*(Rk[1]+2*Re(Rk[2:len]%*%exp(-1.i*omega_k*1:(len-1))))
}
cbind(abs(DFT_ar2)^2,per)[1:30,]


###################################################
### code chunk number 34: eco2.rnw:1100-1130
###################################################
per<-function(x,plot_T)
{
  len<-length(x)
  per<-0:(len/2)
  DFT<-per

  for (k in 0:(len/2))
  {
    cexp <- complex(arg=-(1:len)*2*pi*k/len)
    DFT[k+1]<-sum(cexp*x*sqrt(1/(2*pi*len)))
  }
  per<-abs(DFT)^2
  if (plot_T)
  {
    par(mfrow=c(2,1))
    plot(per,type="l",axes=F,xlab="Frequency",ylab="Periodogram",
    main="Periodogram")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
    "4pi/6","5pi/6","pi"))
    axis(2)
    box()
    plot(log(per),type="l",axes=F,xlab="Frequency",ylab="Log-periodogram",
    main="Log-periodogram")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
    "4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }
  return(list(DFT=DFT,per=per))
}


###################################################
### code chunk number 35: eco2.rnw:1134-1136
###################################################
plot_T<-F
cbind(per(y_ar2,plot_T)$per,abs(DFT_ar2)^2)[1:30,]


###################################################
### code chunk number 36: eco2.rnw:1142-1163
###################################################
len<-300
set.seed(10)
eps<-rnorm(len)
y<-rep(NA,len)
for (k in 100:200)#k<-0
{
	y[k]<-gamma[1:100]%*%eps[k+(0:99)]+gamma[2:100]%*%eps[k-(1:99)]
}
file = paste("z_per_n_lpn.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
plot_T<-F
plot(per(na.exclude(y),plot_T)$per,type="l",axes=F,
main="Periodograms of noise (blue) and lowpassed noise(red)",xlab="",
ylab="",col="red",lwd=2)
lines(per(eps[!is.na(y)],plot_T)$per,col="blue",lwd=2)
K<-length(per(na.exclude(y),plot_T)$per)-1
axis(1, at=c(0,1+1:6*K/6),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 37: z_per_n_lpn.pdf
###################################################
  file = paste("z_per_n_lpn.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Noise (blue) and lowpassed noise (red)", sep = "")
  cat("\\label{z_per_n_lpn}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 38: eco2.rnw:1192-1205
###################################################
len<-100
set.seed(10)
model<-list(order=c(2,0,0),ar=c(1.4,-0.7))
y<-arima.sim(model,n=len)

file = paste("z_acf_ar2.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(3,1))
ts.plot(y,main="Data AR(2)")
acf(y)
acf(y,type="partial")

dev.off()


###################################################
### code chunk number 39: z_acf_ar2.pdf
###################################################
  file = paste("z_acf_ar2.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Data AR(2) and sample acf/pacf", sep = "")
  cat("\\label{z_acf_ar2}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 40: eco2.rnw:1219-1260
###################################################
# false models
model_ar1<-list(order=c(1,0,0))
model_ma2<-list(order=c(0,0,2))
model_ma5<-list(order=c(0,0,5))

y.true<-arima(y,order=model$order)
res_true<-y.true$residuals
y.ar1<-arima(y,order=model_ar1$order)
res_ar1<-y.ar1$residuals
y.ma2<-arima(y,order=model_ma2$order)
res_ma2<-y.ma2$residuals
y.ma5<-arima(y,order=model_ma5$order)
res_ma5<-y.ma5$residuals

file = paste("z_per_ar2_diag.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,1))
plot_T<-F
plot(per(res_true,plot_T)$per,type="l",axes=F,
main="Periodograms of residuals : true AR(2)",xlab="",
ylab="",col="black",lwd=2)
K<-length(per(res_true,plot_T)$per)-1
axis(1, at=c(0,1+1:6*K/6),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
plot(per(res_ar1,plot_T)$per,type="l",axes=F,
main="Periodograms of residuals: false AR(1), MA(2), MA(5)",xlab="",
ylab="",col="red",lwd=2)
lines(per(res_ma2,plot_T)$per,col="green",lwd=2)
lines(per(res_ma5,plot_T)$per,col="blue",lwd=2)
mtext("AR(1)", side = 3, line = -1,at=len/4,col="red")
mtext("MA(2)", side = 3, line = -2,at=len/4,col="green")
mtext("MA(5)", side = 3, line = -3,at=len/4,col="blue")
K<-length(per(res_true,plot_T)$per)-1
axis(1, at=c(0,1+1:6*K/6),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 41: z_per_ar2_diag.pdf
###################################################
  file = paste("z_per_ar2_diag.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Periodograms of true AR(2) (top) and false AR(1), MA(2), MA(5)", sep = "")
  cat("\\label{z_per_ar2_diag}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 42: eco2.rnw:1276-1277
###################################################
polyroot(c(0.7,-1.4,1))


###################################################
### code chunk number 43: eco2.rnw:1280-1281
###################################################
Arg(polyroot(c(0.7,-1.4,1)))



#--------------------------------------------------------

# Section on filters Filters


###################################################
### code chunk number 44: eco2.rnw:1573-1618
###################################################
len<-100
b0<-1
file = paste("z_cos_in_out.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,2))
b1<--1
omega1<-pi/20
y<-rep(0,len)
x<-cos((1:len)*omega1)
y[2:len]<-b0*x[2:len]+b1*x[1:(len-1)]
ts.plot(x,type="l",main=paste("MA(1): b1 = ",b1,
", omega=",round(omega1,3),sep=""),xlab="",
ylab="",col="blue",lwd=1)
lines(y,col="red",lwd=1)
mtext("Input", side = 3, line = -1,at=len/2,col="blue")
mtext("Output", side = 3, line = -2,at=len/2,col="red")
b1<-1
y[2:len]<-b0*x[2:len]+b1*x[1:(len-1)]
ts.plot(y,type="l",main=paste("MA(1):  b1 = ",b1,
", omega=",round(omega1,3),sep=""),xlab="",
ylab="",col="red",lwd=1)
lines(x,col="blue",lwd=1)
mtext("Input", side = 3, line = -1,at=len/2,col="blue")
mtext("Output", side = 3, line = -2,at=len/2,col="red")
b1<--1
omega2<-pi/5
y<-rep(0,len)
x<-cos((1:len)*omega2)
y[2:len]<-b0*x[2:len]+b1*x[1:(len-1)]
ts.plot(x,type="l",main=paste("MA(1):  b1 = ",b1,
", omega=",round(omega2,3),sep=""),xlab="",
ylab="",col="blue",lwd=1)
lines(y,col="red",lwd=1)
mtext("Input", side = 3, line = -1,at=len/2,col="blue")
mtext("Output", side = 3, line = -2,at=len/2,col="red")
b1<-1
y[2:len]<-b0*x[2:len]+b1*x[1:(len-1)]
ts.plot(y,type="l",main=paste("MA(1): b1 = ",b1,
", omega=",round(omega2,3),sep=""),xlab="",
ylab="",col="red",lwd=1)
lines(x,col="blue",lwd=1)
mtext("Input", side = 3, line = -1,at=len/2,col="blue")
mtext("Output", side = 3, line = -2,at=len/2,col="red")
dev.off()


###################################################
### code chunk number 45: z_cos_in_out.pdf
###################################################
  file = paste("z_cos_in_out.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Filter effect of MA(1) on trigonometric input", sep = "")
  cat("\\label{z_cos_in_out}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 46: eco2.rnw:1632-1710
###################################################
file = paste("z_amp_pha_ma1.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,2))
# Resolution of discrete frequency grid
len1<-1001
omega_k<-(0:len1)*pi/len1
# Compute transfer function, amplitude, phase and time-shift
b1<-1
trffkt1<-1+b1*complex(arg=-omega_k)
amplitude1<-abs(trffkt1)
phase1<-Arg(trffkt1)
shift1<--phase1/omega_k
plot(amplitude1,type="l",main=paste("Amplitude MA(1): b1 = ",b1,sep=""),
axes=F,xlab="Frequency",ylab="Amplitude",col="black")
abline(v=(len1-1)/(pi/omega1),lty=2,col="blue")
abline(v=(len1-1)/(pi/omega2),lty=2,col="orange")
mtext(paste("pi/",pi/omega1,sep=""), side = 3, line = -4,
at=len1/20,col="blue")
mtext(paste("pi/",pi/omega2,sep=""), side = 3, line = -6,
at=len1/5,col="orange")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6",
"3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(shift1,type="l",main=paste("Time-shift MA(1): b1 = ",b1,sep=""),
axes=F,xlab="Frequency",ylab="Shift",col="green",
ylim=c(min(na.exclude(shift1))-0.5,max(na.exclude(shift1))+0.5))
abline(v=(len1-1)/(pi/omega1),lty=2,col="blue")
abline(v=(len1-1)/(pi/omega2),lty=2,col="orange")
mtext(paste("pi/",pi/omega1,sep=""), side = 3, line = -4,
at=len1/20,col="blue")
mtext(paste("pi/",pi/omega2,sep=""), side = 3, line = -6,
at=len1/5,col="orange")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6",
"3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
b1<--1
trffkt2<-1+b1*complex(arg=-omega_k)
amplitude2<-abs(trffkt2)
phase2<-Arg(trffkt2)
shift2<--phase2/omega_k
plot(amplitude2,type="l",main=paste("Amplitude MA(1): b1 = ",b1,sep=""),
axes=F,xlab="Frequency",ylab="Amplitude",col="black")
abline(v=(len1-1)/(pi/omega1),lty=2,col="blue")
abline(v=(len1-1)/(pi/omega2),lty=2,col="orange")
mtext(paste("pi/",pi/omega1,sep=""), side = 3, line = -4,
at=len1/20,col="blue")
mtext(paste("pi/",pi/omega2,sep=""), side = 3, line = -6,
at=len1/5,col="orange")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6",
"3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(shift2,type="l",main=paste("Time-shift MA(1): b1 = ",b1,sep=""),
axes=F,xlab="Frequency",ylab="Shift",col="green")
abline(v=(len1-1)/(pi/omega1),lty=2,col="blue")
abline(v=(len1-1)/(pi/omega2),lty=2,col="orange")
mtext(paste("pi/",pi/omega1,sep=""), side = 3, line = -4,
at=len1/20,col="blue")
mtext(paste("pi/",pi/omega2,sep=""), side = 3, line = -6,
at=len1/5,col="orange")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6",
"3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()
# Compute table with shifts and amplitudes in omega=pi/20 and omega=pi/5
table_amp_shift<-
rbind(c(amplitude1[(len1-1)/(pi/omega1)+1],shift1[(len1-1)/(pi/omega1)+1]),
c(amplitude2[(len1-1)/(pi/omega1)+1],shift2[(len1-1)/(pi/omega1)+1]),
c(amplitude1[(len1-1)/(pi/omega2)+1],shift1[(len1-1)/(pi/omega2)+1]),
c(amplitude2[(len1-1)/(pi/omega2)+1],shift2[(len1-1)/(pi/omega2)+1]))
dimnames(table_amp_shift)[[2]]<-c("Amplitude","Shift")
dimnames(table_amp_shift)[[1]]<-c(paste("pi/",pi/omega1,",b1=1",sep=""),
paste("pi/",pi/omega1,",b1=-1",sep=""),
paste("pi/",pi/omega2,",b1=1",sep=""),paste("pi/",pi/omega2,",b1=-1",sep=""))


###################################################
### code chunk number 47: z_amp_pha_ma1.pdf
###################################################
  file = paste("z_amp_pha_ma1.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Amplitude and time-shifts of MA(1)-filters", sep = "")
  cat("\\label{z_amp_pha_ma1}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 48: table_amp_shift
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(table_amp_shift, dec = 1,digits=rep(3,dim(table_amp_shift)[2]+1),
  paste("Amplitudes and time-shifts of MA(1)-filters in pi/20 and pi/5",sep=""),
  label=paste("table_amp_shift",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 49: eco2.rnw:1758-1779
###################################################
# We verify the claim for omega=pi/20 and b1=1
amp_scaling<-table_amp_shift[1,1]
shift<-table_amp_shift[1,2]
len<-100
b0<-1
b1<-1
omega1<-pi/20
y<-rep(0,len)
x<-cos((1:len)*omega1)
y[2:len]<-b0*x[2:len]+b1*x[1:(len-1)]
# Here We scale and shift the input signal
z<-amp_scaling*cos((1:len-shift)*omega1)
file = paste("z_out_ssout.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
ts.plot(y,type="l",main=paste("Output vs. scaled and shifted input: b0 = b1 = ",
b1,", omega=",round(omega1,3),sep=""),xlab="",
ylab="",col="red",lwd=1)
lines(z,col="green",lwd=1)
mtext("Scaled and shifted input", side = 3, line = -1,at=len/2,col="green")
mtext("Output", side = 3, line = -2,at=len/2,col="red")
dev.off()


###################################################
### code chunk number 50: z_out_ssout.pdf
###################################################
  file = paste("z_out_ssout.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Output of MA(1) vs. scaled and shifted input", sep = "")
  cat("\\label{z_out_ssout}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 51: eco2.rnw:1800-1850
###################################################
len<-1000
a1<-0.9
file = paste("z_cos_in_out_ar1.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
len<-1000
par(mfrow=c(2,2))
b0<-1
omega1<-pi/20
y<-rep(0,len)
x<-cos((1:len)*omega1)
for (i in 2:len)
  y[i]<-b0*x[i]+a1*y[i-1]
ts.plot(y[(len-99):len],type="l",main=paste("AR(1): b0 = ",
b0, ", a1 = ",a1,", omega=",round(omega1,3),sep=""),xlab="",
ylab="",col="red",lwd=1)
lines(x[(len-99):len],col="blue",lwd=1)
mtext("Input", side = 3, line = -1,at=50,col="blue")
mtext("Output", side = 3, line = -2,at=50,col="red")
b0<-0.1
y<-rep(0,len)
for (i in 2:len)
  y[i]<-b0*x[i]+a1*y[i-1]
ts.plot(x[(len-99):len],type="l",main=paste("AR(1): b0 = ",
b0, ", a1 = ",a1,", omega=",round(omega1,3),sep=""),xlab="",
ylab="",col="blue",lwd=1)
lines(y[(len-99):len],col="red",lwd=1)
mtext("Input", side = 3, line = -1,at=50,col="blue")
mtext("Output", side = 3, line = -2,at=50,col="red")
b0<-1
omega2<-pi/5
y<-rep(0,len)
x<-cos((1:len)*omega2)
for (i in 2:len)
  y[i]<-b0*x[i]+a1*y[i-1]
ts.plot(y[(len-99):len],type="l",main=paste("AR(1): b0 = ",
b0, ", a1 = ",a1,", omega=",round(omega2,3),sep=""),xlab="",
ylab="",col="red",lwd=1)
lines(x[(len-99):len],col="blue",lwd=1)
mtext("Input", side = 3, line = -1,at=50,col="blue")
mtext("Output", side = 3, line = -2,at=50,col="red")
b0<-0.1
for (i in 2:len)
  y[i]<-b0*x[i]+a1*y[i-1]
ts.plot(x[(len-99):len],type="l",main=paste("AR(1): b0 = ",
b0, ", a1 = ",a1,", omega=",round(omega2,3),sep=""),xlab="",
ylab="",col="blue",lwd=1)
lines(y[(len-99):len],col="red",lwd=1)
mtext("Input", side = 3, line = -1,at=50,col="blue")
mtext("Output", side = 3, line = -2,at=50,col="red")
dev.off()


###################################################
### code chunk number 52: z_cos_in_out_ar1.pdf
###################################################
  file = paste("z_cos_in_out_ar1.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Filter effect of AR(1) on trigonometric input", sep = "")
  cat("\\label{z_cos_in_out_ar1}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 53: eco2.rnw:1866-1948
###################################################
file = paste("z_amp_pha_ar1.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,2))
# Resolution of discrete frequency grid
len1<-1001
omega_k<-(0:len1)*pi/len1
# Compute transfer function, amplitude, phase and time-shift
a1<-0.9
b0<-1
trffkt1<-b0/(1-a1*complex(arg=-omega_k))
amplitude1<-abs(trffkt1)
phase1<-Arg(trffkt1)
shift1<--phase1/omega_k
plot(amplitude1,type="l",
main=paste("Amplitude AR(1): b0 = ",b0," ,a1 = ",a1,sep=""),
axes=F,xlab="Frequency",ylab="Amplitude",col="black")
abline(v=(len1-1)/(pi/omega1),lty=2,col="blue")
abline(v=(len1-1)/(pi/omega2),lty=2,col="orange")
mtext(paste("pi/",pi/omega1,sep=""), side = 3,
line = -4,at=len1/20,col="blue")
mtext(paste("pi/",pi/omega2,sep=""), side = 3,
line = -6,at=len1/5,col="orange")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6",
"3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(shift1,type="l",
main=paste("Time-shift AR(1): b0 = ",b0," ,a1 = ",a1,sep=""),
axes=F,xlab="Frequency",ylab="Shift",col="green")
abline(v=(len1-1)/(pi/omega1),lty=2,col="blue")
abline(v=(len1-1)/(pi/omega2),lty=2,col="orange")
mtext(paste("pi/",pi/omega1,sep=""), side = 3,
line = -4,at=len1/20,col="blue")
mtext(paste("pi/",pi/omega2,sep=""), side = 3,
line = -6,at=len1/5,col="orange")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6",
"3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
b0<-0.1
trffkt2<-b0/(1-a1*complex(arg=-omega_k))
amplitude2<-abs(trffkt2)
phase2<-Arg(trffkt2)
shift2<--phase2/omega_k
plot(amplitude2,type="l",
main=paste("Amplitude AR(1): b0 = ",b0," ,a1 = ",a1,sep=""),
axes=F,xlab="Frequency",ylab="Amplitude",col="black")
abline(v=(len1-1)/(pi/omega1),lty=2,col="blue")
abline(v=(len1-1)/(pi/omega2),lty=2,col="orange")
mtext(paste("pi/",pi/omega1,sep=""), side = 3,
line = -4,at=len1/20,col="blue")
mtext(paste("pi/",pi/omega2,sep=""), side = 3,
line = -6,at=len1/5,col="orange")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6",
"3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(shift2,type="l",
main=paste("Time-shift AR(1): b0 = ",b0," ,a1 = ",a1,sep=""),
axes=F,xlab="Frequency",ylab="Shift",col="green")
abline(v=(len1-1)/(pi/omega1),lty=2,col="blue")
abline(v=(len1-1)/(pi/omega2),lty=2,col="orange")
mtext(paste("pi/",pi/omega1,sep=""), side = 3,
line = -4,at=len1/20,col="blue")
mtext(paste("pi/",pi/omega2,sep=""), side = 3,
line = -6,at=len1/5,col="orange")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6",
"3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()
# Compute table with shifts and amplitudes in omega=pi/20 and omega=pi/5
table_amp_shift_ar1<-rbind(c(amplitude1[(len1-1)/(pi/omega1)+1],
shift1[(len1-1)/(pi/omega1)+1]),
c(amplitude2[(len1-1)/(pi/omega1)+1],shift2[(len1-1)/(pi/omega1)+1]),
c(amplitude1[(len1-1)/(pi/omega2)+1],shift1[(len1-1)/(pi/omega2)+1]),
c(amplitude2[(len1-1)/(pi/omega2)+1],shift2[(len1-1)/(pi/omega2)+1]))
dimnames(table_amp_shift_ar1)[[2]]<-c("Amplitude","Shift")
dimnames(table_amp_shift_ar1)[[1]]<-c(paste("pi/",pi/omega1,",b0=1",sep=""),
paste("pi/",pi/omega1,",b0=0.1",sep=""),
paste("pi/",pi/omega2,",b0=1",sep=""),paste("pi/",pi/omega2,",b0=0.1",sep=""))


###################################################
### code chunk number 54: z_amp_pha_ar1.pdf
###################################################
  file = paste("z_amp_pha_ar1.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Amplitude and time-shifts of AR(1)-filters", sep = "")
  cat("\\label{z_amp_pha_ar1}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 55: table_amp_shift_ar1
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(table_amp_shift_ar1, dec = 1,digits=rep(3,dim(table_amp_shift)[2]+1),
  paste("Amplitudes and time-shifts of AR(1)-filters in pi/20 and pi/5",sep=""),
  label=paste("table_amp_shift_ar1",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 56: eco2.rnw:1994-2015
###################################################
# We verify the claim for omega=pi/20 and b1=1
amp_scaling<-table_amp_shift_ar1[1,1]
shift<-table_amp_shift_ar1[1,2]
len<-200
b0<-1
omega1<-pi/20
y<-rep(0,len)
x<-cos((1:len)*omega1)
for (i in 2:len)
  y[i]<-b0*x[i]+a1*y[i-1]
# Here We scale and shift the input signal
z<-amp_scaling*cos((1:len-shift)*omega1)
file = paste("z_out_ssout_ar1.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
ts.plot(y,type="l",main=paste("Output vs. scaled and shifted input: b0 = ",b0,
", a1 = ", a1,", omega=",round(omega1,3),sep=""),xlab="",
ylab="",col="red",lwd=1)
lines(z,col="green",lwd=1)
mtext("Scaled and shifted input", side = 3, line = -1,at=len/2,col="green")
mtext("Output", side = 3, line = -2,at=len/2,col="red")
dev.off()


###################################################
### code chunk number 57: z_out_ssout_ar1.pdf
###################################################
  file = paste("z_out_ssout_ar1.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Output of AR(1) vs. scaled and shifted input", sep = "")
  cat("\\label{z_out_ssout_ar1}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 58: eco2.rnw:2039-2086
###################################################
len<-120
b0<-1
b12<--1
file = paste("z_cos_in_out_sd.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,2))
omega1<-pi/20
y<-rep(0,len)
x<-cos((1:len)*omega1)
y[13:len]<-b0*x[13:len]+b12*x[1:(len-12)]
ts.plot(y,type="l",main=paste("Seasonal difference MA(12): b0 =1, b12 = ",
b12,", omega=",round(omega1,3),sep=""),xlab="",
ylab="",col="red",lwd=1)
lines(x,col="blue",lwd=1)
mtext("Input", side = 3, line = -1,at=len/2,col="blue")
mtext("Output", side = 3, line = -2,at=len/2,col="red")
omega2<-pi/5
y<-rep(0,len)
x<-cos((1:len)*omega2)
y[13:len]<-b0*x[13:len]+b12*x[1:(len-12)]
ts.plot(y,type="l",main=paste("Seasonal difference MA(12): b0 =1, b12 = ",
b12,", omega=",round(omega1,3),sep=""),xlab="",
ylab="",col="red",lwd=1)
lines(x,col="blue",lwd=1)
mtext("Input", side = 3, line = -1,at=len/2,col="blue")
mtext("Output", side = 3, line = -2,at=len/2,col="red")
omega3<-pi/6
y<-rep(0,len)
x<-cos((1:len)*omega3)
y[13:len]<-b0*x[13:len]+b12*x[1:(len-12)]
ts.plot(x,type="l",main=paste("Seasonal difference MA(12): b0 =1, b12 = ",
b12,", omega=",round(omega1,3),sep=""),xlab="",
ylab="",col="blue",lwd=1)
lines(y,col="red",lwd=1)
mtext("Input", side = 3, line = -1,at=len/2,col="blue")
mtext("Output", side = 3, line = -2,at=len/2,col="red")
omega4<-4*pi/6
y<-rep(0,len)
x<-cos((1:len)*omega4)
y[13:len]<-b0*x[13:len]+b12*x[1:(len-12)]
ts.plot(x,type="l",main=paste("Seas. diff. MA(12), omega=",
round(omega1,3),sep=""),xlab="",
ylab="",col="blue",lwd=1)
lines(y,col="red",lwd=1)
mtext("Input", side = 3, line = -1,at=len/2,col="blue")
mtext("Output", side = 3, line = -2,at=len/2,col="red")
dev.off()


###################################################
### code chunk number 59: z_cos_in_out_sd.pdf
###################################################
  file = paste("z_cos_in_out_sd.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Filter effect of seasonal difference MA(12) on trigonometric inputs", sep = "")
  cat("\\label{z_cos_in_out_sd}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 60: eco2.rnw:2103-2141
###################################################
file = paste("z_amp_pha_ma12.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
# Resolution of discrete frequency grid
len1<-1200
omega_k<-(0:len1)*pi/len1
# Compute transfer function, amplitude, phase and time-shift
b1<-1
b12<--1
trffkt<-1+b12*complex(arg=-12*omega_k)
amplitude<-abs(trffkt)
phase<-Arg(trffkt)
shift<--phase/omega_k
par(mfrow=c(2,1))
plot(amplitude,type="l",main="Amplitude seasonal difference MA(12)",
axes=F,xlab="Frequency",ylab="Amplitude",col="black")
abline(v=(len1-1)/(pi/omega1),lty=2,col="blue")
abline(v=(len1-1)/(pi/omega2),lty=2,col="orange")
mtext(paste("pi/",pi/omega1,sep=""), side = 3, line = -2,
at=len1/20,col="blue")
mtext(paste("pi/",pi/omega2,sep=""), side = 3, line = -4,
at=len1/5,col="orange")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
plot(shift,type="l",main="Time shift seasonal difference MA(12)",
axes=F,xlab="Frequency",ylab="Amplitude",col="black")
abline(v=(len1-1)/(pi/omega1),lty=2,col="blue")
abline(v=(len1-1)/(pi/omega2),lty=2,col="orange")
mtext(paste("pi/",pi/omega1,sep=""), side = 3, line = -2,
at=len1/20,col="blue")
mtext(paste("pi/",pi/omega2,sep=""), side = 3, line = -4,
at=len1/5,col="orange")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 61: z_amp_pha_ma12.pdf
###################################################
  file = paste("z_amp_pha_ma12.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Amplitude and time-shifts of seasonal difference MA(12)-filter", sep = "")
  cat("\\label{z_amp_pha_ma12}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 62: eco2.rnw:2272-2334
###################################################
set.seed(10)
len<-100
x<-rnorm(len)
y<-c(0,diff(x))
# Resolution of discrete frequency grid
len1<-len/2
omega_k<-(0:len1)*pi/len1
# Compute transfer function, amplitude, phase and time-shift
b1<--1
trffkt<-1+b1*complex(arg=-omega_k)
amplitude<-abs(trffkt)
phase<-Arg(trffkt)
shift<-phase/omega_k
plot_T<-F

file = paste("z_convo_diff.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,2))
plot(amplitude,type="l",main="Convolution: ordinary difference filter",
axes=F,xlab="Frequency",ylab="periodograms",col="black",ylim=c(0,2))
lines(amplitude^2*per(x,plot_T)$per,col="green")
lines(per(y,plot_T)$per,col="red",lty=2)
lines(per(x,plot_T)$per,col="blue")
mtext("Amplitude", side = 3, line = -1,at=len1/2,col="black")
mtext("Convolution", side = 3, line = -2,at=len1/2,col="green")
mtext("Periodogram output", side = 3, line = -3,at=len1/2,col="red")
mtext("Periodogram input", side = 3, line = -4,at=len1/2,col="blue")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
ts.plot(x,col="blue",main="Input and output series",ylab="")
lines(y,col="red")
mtext("Input series", side = 3, line = -1,at=len/2,col="blue")
mtext("Output series", side = 3, line = -2,at=len/2,col="red")
plot(Re(trffkt),type="l",main="Convolution: real-parts of DFT",
axes=F,xlab="Frequency",ylab="Real parts of DFT",col="black",ylim=c(-1,2))
lines(Re(trffkt*per(x,plot_T)$DFT),col="green")
lines(Re(per(y,plot_T)$DFT),col="red",lty=2)
lines(Re(per(x,plot_T)$DFT),col="blue")
mtext("Re(transfer function)", side = 3, line = -1,at=len1/2,col="black")
mtext("Re(convolution of DFT)", side = 3, line = -2,at=len1/2,col="green")
mtext("Re(DFT output)", side = 3, line = -3,at=len1/2,col="red")
mtext("Re(DFT input)", side = 3, line = -4,at=len1/2,col="blue")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
plot(Im(trffkt),type="l",main="Convolution: imaginary-parts of DFT",
axes=F,xlab="Frequency",ylab="Imaginary parts of DFT",col="black",ylim=c(-1,1))
lines(Im(trffkt*per(x,plot_T)$DFT),col="green")
lines(Im(per(y,plot_T)$DFT),col="red",lty=2)
lines(Im(per(x,plot_T)$DFT),col="blue")
mtext("Im(transfer function)", side = 3, line = -1,at=len1/2,col="black")
mtext("Im(convolution of DFT)", side = 3, line = -2,at=len1/2,col="green")
mtext("Im(DFT output)", side = 3, line = -3,at=len1/2,col="red")
mtext("Im(DFT input)", side = 3, line = -4,at=len1/2,col="blue")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 63: z_convo_diff.pdf
###################################################
  file = paste("z_convo_diff.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Convolution of difference filter (green), periodogram of input (blue) and output (red)
  signals and amplitude function of ordinary differences (black)", sep = "")
  cat("\\label{z_convo_diff}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 64: eco2.rnw:2361-2433
###################################################
set.seed(10)
len2<-1000
a1<-0.95
b0<-1-a1
omega1<-pi/20
y1<-rep(0,len2)
x1<-rnorm(len2)
for (i in 2:len2)
  y1[i]<-b0*x1[i]+a1*y1[i-1]
# The first 900 observations are discarded (initialization)
len<-100
y<-y1[(len2-len+1):len2]
x<-x1[(len2-len+1):len2]
# Resolution of discrete frequency grid
len1<-len/2
omega_k<-(0:len1)*pi/len1
# Compute transfer function, amplitude, phase and time-shift
trffkt<-b0/(1-a1*complex(arg=-omega_k))
amplitude<-abs(trffkt)
phase<-Arg(trffkt)
shift<-phase/omega_k
plot_T<-F

file = paste("z_convo_ar1.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,2))
plot(amplitude,type="l",main=paste("Convolution AR(1)-filter: b0 = ",
b0," a1 = ",a1,sep=""),
axes=F,xlab="Frequency",ylab="periodograms",col="black",ylim=c(0,1))
lines(amplitude^2*per(x,plot_T)$per,col="green")
lines(per(y,plot_T)$per,col="red",lty=2)
lines(per(x,plot_T)$per,col="blue")
mtext("Amplitude", side = 3, line = -1,at=len1/2,col="black")
mtext("Convolution", side = 3, line = -2,at=len1/2,col="green")
mtext("Periodogram output", side = 3, line = -3,at=len1/2,col="red")
mtext("Periodogram input", side = 3, line = -4,at=len1/2,col="blue")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
ts.plot(x,col="blue",main="Input and output series",ylab="")
lines(y,col="red")
mtext("Input series", side = 3, line = -1,at=len/2,col="blue")
mtext("Output series", side = 3, line = -2,at=len/2,col="red")
plot(Re(trffkt),type="l",main="Convolution: real-parts of DFT",
axes=F,xlab="Frequency",ylab="Real parts of DFT",col="black",ylim=c(-1,1))
lines(Re(trffkt*per(x,plot_T)$DFT),col="green")
lines(Re(per(y,plot_T)$DFT),col="red",lty=2)
lines(Re(per(x,plot_T)$DFT),col="blue")
mtext("Re(transfer function)", side = 3, line = -1,at=len1/2,col="black")
mtext("Re(convolution of DFT)", side = 3, line = -2,at=len1/2,col="green")
mtext("Re(DFT output)", side = 3, line = -3,at=len1/2,col="red")
mtext("Re(DFT input)", side = 3, line = -4,at=len1/2,col="blue")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
plot(Im(trffkt),type="l",main="Convolution: imaginary-parts of DFT",
axes=F,xlab="Frequency",ylab="Imaginary parts of DFT",col="black",ylim=c(-1,1))
lines(Im(trffkt*per(x,plot_T)$DFT),col="green")
lines(Im(per(y,plot_T)$DFT),col="red",lty=2)
lines(Im(per(x,plot_T)$DFT),col="blue")
mtext("Im(transfer function)", side = 3, line = -1,at=len1/2,col="black")
mtext("Im(convolution of DFT)", side = 3, line = -2,at=len1/2,col="green")
mtext("Im(DFT output)", side = 3, line = -3,at=len1/2,col="red")
mtext("Im(DFT input)", side = 3, line = -4,at=len1/2,col="blue")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 65: z_convo_ar1.pdf
###################################################
  file = paste("z_convo_ar1.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Convolution of AR(1)-filter (green), periodogram of input (blue) and output (red)
  signals and amplitude function of AR(1) filter (black)", sep = "")
  cat("\\label{z_convo_ar1}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 66: eco2.rnw:2460-2461
###################################################
polyroot(c(1,rep(0,11),-1))


###################################################
### code chunk number 67: eco2.rnw:2464-2465
###################################################
abs(polyroot(c(1,rep(0,11),-1)))


###################################################
### code chunk number 68: eco2.rnw:2472-2546
###################################################
set.seed(10)
len2<-1000
a12<-1
b0<-1
omega1<-pi/20
y1<-rep(0,len2)
x1<-rnorm(len2)
for (i in 13:len2)
  y1[i]<-b0*x1[i]+a12*y1[i-12]
# The first 900 observations are discarded (initialization)
len<-120
y<-y1[(len2-len+1):len2]
x<-x1[(len2-len+1):len2]
# Resolution of discrete frequency grid
len1<-len/2
omega_k<-(0:len1)*pi/len1
# Compute transfer function, amplitude, phase and time-shift
trffkt<-b0/(1-a12*complex(arg=-12*omega_k))
# The transfer function is infinite at the seasonal frequencies
trffkt[c(1,1+1:6*(len1/6))]<-Inf
amplitude<-abs(trffkt)
phase<-Arg(trffkt)
shift<-phase/omega_k
plot_T<-F

file = paste("z_convo_ar12.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,2))
plot(amplitude,type="l",main="Convolution seasonal AR(12)",
axes=F,xlab="Frequency",ylab="periodograms",col="black",ylim=c(0,5))
lines(amplitude^2*per(x,plot_T)$per,col="green")
lines(per(y,plot_T)$per,col="red",lty=2)
lines(per(x,plot_T)$per,col="blue")
mtext("Amplitude", side = 3, line = -1,at=len1/2,col="black")
mtext("Convolution", side = 3, line = -2,at=len1/2,col="green")
mtext("Periodogram output", side = 3, line = -3,at=len1/2,col="red")
mtext("Periodogram input", side = 3, line = -4,at=len1/2,col="blue")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
ts.plot(x,col="blue",main="Input and output series",ylab="",
ylim=c(min(y),max(y)))
lines(y,col="red")
mtext("Input series", side = 3, line = -1,at=len/2,col="blue")
mtext("Output series", side = 3, line = -2,at=len/2,col="red")
plot(Re(trffkt),type="l",main="Convolution: real-parts of DFT",
axes=F,xlab="Frequency",ylab="Real parts of DFT",col="black",ylim=c(-5,5))
lines(Re(trffkt*per(x,plot_T)$DFT),col="green")
lines(Re(per(y,plot_T)$DFT),col="red",lty=2)
lines(Re(per(x,plot_T)$DFT),col="blue")
mtext("Re(transfer function)", side = 3, line = -1,at=len1/2,col="black")
mtext("Re(convolution of DFT)", side = 3, line = -2,at=len1/2,col="green")
mtext("Re(DFT output)", side = 3, line = -3,at=len1/2,col="red")
mtext("Re(DFT input)", side = 3, line = -4,at=len1/2,col="blue")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
plot(Im(trffkt),type="l",main="Convolution: imaginary-parts of DFT",
axes=F,xlab="Frequency",ylab="Imaginary parts of DFT",col="black",ylim=c(-5,5))
lines(Im(trffkt*per(x,plot_T)$DFT),col="green")
lines(Im(per(y,plot_T)$DFT),col="red",lty=2)
lines(Im(per(x,plot_T)$DFT),col="blue")
mtext("Im(transfer function)", side = 3, line = -1,at=len1/2,col="black")
mtext("Im(convolution of DFT)", side = 3, line = -2,at=len1/2,col="green")
mtext("Im(DFT output)", side = 3, line = -3,at=len1/2,col="red")
mtext("Im(DFT input)", side = 3, line = -4,at=len1/2,col="blue")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 69: z_convo_ar12.pdf
###################################################
  file = paste("z_convo_ar12.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Convolution of unstable seasonal AR(12)-filter (green), periodogram of input (blue) and output (red)
  signals and amplitude function of seasonal AR(12) (black)", sep = "")
  cat("\\label{z_convo_ar12}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 70: eco2.rnw:2661-2713
###################################################
len<-1000
set.seed(10)
u1<-arima.sim(list(ar=-0.9),n=len)
set.seed(10)
u2<-arima.sim(list(ar=0.),n=len)
set.seed(10)
u3<-arima.sim(list(ar=0.9),n=len)
# Compute transferfunction of difference filter
len1<-length(per(u1,plot_T)$per)-1
trffkt_diff<-(1-exp(1.i*(0:len1)*pi/len1))
# Amplitude function of difference filter
amp_diff<-abs(trffkt_diff)

file = paste("z_convo_bias.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
plot_T<-F
par(mfrow=c(3,1))
per_pseudo<-per(u1,plot_T)$per
# compute non-stationary x
x1<-cumsum(u1)
per_direct<-amp_diff^2*per(x1,plot_T)$per
plot(per_pseudo,type="l",main="Bias effect: a1=-0.9",
axes=F,xlab="Frequency",ylab="periodograms",col="black",ylim=c(0,max(per_direct)))
lines(per_direct,col="blue")
mtext("Direct (biased)", side = 3, line = -1,at=len1/2,col="blue")
mtext("Pseudo (unbiased)", side = 3, line = -2,at=len1/2,col="black")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
per_pseudo<-per(u2,plot_T)$per
x2<-cumsum(u2)
per_direct<-amp_diff^2*per(x2,plot_T)$per
plot(per_pseudo,type="l",main="Bias effect: a1=0",
axes=F,xlab="Frequency",ylab="periodograms",col="black",ylim=c(0,max(per_direct)))
lines(per_direct,col="blue")
mtext("Direct (biased)", side = 3, line = -1,at=len1/2,col="blue")
mtext("Pseudo (unbiased)", side = 3, line = -2,at=len1/2,col="black")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
per_pseudo<-per(u3,plot_T)$per
x3<-cumsum(u3)
per_direct<-amp_diff^2*per(x3,plot_T)$per
plot(per_pseudo,type="l",main="Bias effect: a1=0.9",
axes=F,xlab="Frequency",ylab="periodograms",col="black",ylim=c(0,max(per_direct)))
lines(per_direct,col="blue")
mtext("Direct (biased)", side = 3, line = -1,at=len1/2,col="blue")
mtext("Pseudo (unbiased)", side = 3, line = -2,at=len1/2,col="black")
axis(1,at=c(0,1:6*len1/6+1),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 71: z_convo_bias.pdf
###################################################
  file = paste("z_convo_bias.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Bias effect depending on DGP: a1=-0.9 (top), a1=0 (middle) and a1=0.9 (bottom)", sep = "")
  cat("\\label{z_convo_bias}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


#-------------------------
# DFA


###################################################
### code chunk number 72: eco2.rnw:2846-2879
###################################################
# This function computes mean-square DFA-solutions
# L is the length of the MA filter,
# periodogram is the frequency weighting function in the DFA
# Gamma is the transferfunction of the symmetric filter (target) and
# Lag is the lag-parameter: Lag=0 implies real-time filtering, Lag=L/2
#     implies symmetric filter
# The function returns optimal coefficients as well as the transfer function of the
#     optimized real-time filter
dfa_ms<-function(L,periodogram,Lag,Gamma)
{

  K<-length(periodogram)-1
  X<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)*sqrt(periodogram)
  X_y<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)
  for (l in 2:L)          #l<-L<-21
  {
    X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))+
    1.i*sin((l-1-Lag)*pi*(0:(K))/(K)))*sqrt(periodogram))
    X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))+
    1.i*sin((l-1-Lag)*pi*(0:(K))/(K))))
  }
  xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
# MA-Filtercoefficients
  b<-as.vector(solve(xtx)%*%(t(Re(X_y))%*%(Gamma*periodogram)))
# Transferfunction
  trffkt<-1:(K+1)
  trffkt[1]<-sum(b)
  for (k in 1:(K))#k<-1
  {
    trffkt[k+1]<-(b%*%exp(1.i*k*(0:(length(b)-1))*pi/(K)))
  }
  return(list(b=b,trffkt=trffkt))
}


###################################################
### code chunk number 73: eco2.rnw:2911-2952
###################################################
# Generate series
set.seed(10)
len<-120
a_vec<-c(0.9,0,-0.9)
x<-matrix(nrow=len,ncol=3)
plot_T<-F
yhat<-x
periodogram<-matrix(ncol=3,nrow=len/2+1)
trffkt<-periodogram
# Generate series
for (i in 1:3)
{
  set.seed(10)
  x[,i]<-arima.sim(list(ar=a_vec[i]),n=len)
}
# Specify filter settings
L<-12
Lag<-0
Gamma<-c(1,(1:(len/2))<len/12)
b<-matrix(nrow=L,ncol=3)
# Compute real-time filters
for (i in 1:3)
{
  periodogram[,i]<-per(x[,i],plot_T)$per
# Optimize filters
  filt<-dfa_ms(L,periodogram[,i],Lag,Gamma)
  trffkt[,i]<-filt$trffkt
  b[,i]<-filt$b
# Compute outputs
  for (j in L:len)
    yhat[j,i]<-filt$b%*%x[j:(j-L+1),i]
}

file = paste("z_dfa_ar1_output.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,1))
for (i in 1:3)
{
  ts.plot(x[,i],main=paste("a1 = ",a_vec[i],sep=""),col="blue")
  lines(yhat[,i],col="red")
}


dev.off()


###################################################
### code chunk number 74: z_dfa_ar1_output.pdf
###################################################
  file = paste("z_dfa_ar1_output", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Inputs (blue) and real-time outputs (red) for a1=0.9 (top), a1=0 (middle) and a1=-0.9 (bottom)", sep = "")
  cat("\\label{z_dfa_ar1_output}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 75: eco2.rnw:2967-3013
###################################################
omega_k<-pi*0:(len/2)/(len/2)
file = paste("z_dfa_ar1_amp_shift.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,2))
amp<-abs(trffkt)
shift<-Arg(trffkt)/omega_k
plot(amp[,1],type="l",main="Amplitude functions",
axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
lines(amp[,2],col="orange")
lines(amp[,3],col="green")
lines(Gamma,col="violet")
mtext("Amplitude a1=0.9", side = 3, line = -1,at=len/4,col="black")
mtext("Amplitude a1=0", side = 3, line = -2,at=len/4,col="orange")
mtext("Amplitude a1=-0.9", side = 3, line = -3,at=len/4,col="green")
mtext("Target", side = 3, line = -4,at=len/4,col="violet")
axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
plot(shift[,1],type="l",main="Time-shifts",
axes=F,xlab="Frequency",ylab="Shift",col="black",
ylim=c(0,max(na.exclude(shift[,3]))))
lines(shift[,2],col="orange")
lines(shift[,3],col="green")
lines(rep(0,len/2+1),col="violet")
mtext("Shift a1=0.9", side = 3, line = -1,at=len/4,col="black")
mtext("Shift a1=0", side = 3, line = -2,at=len/4,col="orange")
mtext("Shift a1=-0.9", side = 3, line = -3,at=len/4,col="green")
mtext("Target", side = 3, line = -4,at=len/4,col="violet")
axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
plot(periodogram[,1],type="l",main="Periodograms",
axes=F,xlab="Frequency",ylab="Periodogram",col="black",
ylim=c(0,max(periodogram[,3])/6))
lines(periodogram[,2],col="orange")
lines(periodogram[,3],col="green")
mtext("Periodogram a1=0.9", side = 3, line = -1,at=len/4,col="black")
mtext("Periodogram a1=0", side = 3, line = -2,at=len/4,col="orange")
mtext("Periodogram a1=-0.9", side = 3, line = -3,at=len/4,col="green")
axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 76: z_dfa_ar1_amp_shift.pdf
###################################################
  file = paste("z_dfa_ar1_amp_shift.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Amplitude (top left), time-shifts (top-right) and periodograms (bottom left) for
  a1=0.9 (black), a1=0 (orange) and a1=-0.9 (green)", sep = "")
  cat("\\label{z_dfa_ar1_amp_shift}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 77: eco2.rnw:3045-3059
###################################################
# Compute criterion values for all (3*3) combinations of filters and periodograms
perf_mat<-matrix(nrow=3,ncol=3)
for (i in 1:3)
{
  for (j in 1:3)
  {
# Criterion value
    perf_mat[j,i]<-(2*pi/length(Gamma))*
    abs(Gamma-trffkt[,i])^2%*%periodogram[,j]
  }
}
dimnames(perf_mat)[[2]]<-c("filter: a1=0.9","filter: a1=0","filter: a1=-0.9")
dimnames(perf_mat)[[1]]<-c("periodogram: a1=0.9",
" periodogram: a1=0","periodogram: a1=-0.9")


###################################################
### code chunk number 78: table_perf_mat
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(perf_mat, dec = 1,digits=rep(3,dim(perf_mat)[2]+1),
  paste("Criterion values for all 3*3 combinations of periodograms and filters",sep=""),
  label=paste("perf_mat",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 79: eco2.rnw:3081-3152
###################################################
# Define all relevant variables
set.seed(10)
lenh<-2000
len<-120
a_vec<-c(0.9,0,-0.9)
xh<-matrix(nrow=lenh,ncol=length(a_vec))
x<-matrix(nrow=len,ncol=length(a_vec))
plot_T<-F
yhat<-x
y<-x
periodogram<-matrix(ncol=3,nrow=len/2+1)
trffkt<-periodogram
# Generate series for each process
for (i in 1:length(a_vec))
{
  set.seed(10)
  xh[,i]<-arima.sim(list(ar=a_vec[i]),n=lenh)
}
# We extract 120 observations in the midddle of xh: these will be used
# for real-time filtering
x<-xh[lenh/2+(-len/2):((len/2)-1),]
# Compute the coefficients of the symmetric target filter
cutoff<-pi/6
# Order of approximation
ord<-1000
gamma<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))
# Compute the outputs yt of the symmetric target filter
for (i in 1:length(a_vec))
{
  for (j in 1:120)
  {
    y[j,i]<-gamma[1:900]%*%xh[lenh/2+(-len/2)-1+(j:(j-899)),i]+
    gamma[2:900]%*%xh[lenh/2+(-len/2)+(j:(j+898)),i]
  }
}
# Specify real-time filter settings
L<-12
Lag<-0
Gamma<-c(1,(1:(len/2))<len/12)
b<-matrix(nrow=L,ncol=3)
# Compute real-time filters
for (i in 1:3)
{
  periodogram[,i]<-per(x[,i],plot_T)$per
# Optimize filters
  filt<-dfa_ms(L,periodogram[,i],Lag,Gamma)
  trffkt[,i]<-filt$trffkt
  b[,i]<-filt$b
# Compute outputs: We can use the longer series in order to obtain
# outputs for j=1,...,len
  for (j in 1:len)
    yhat[j,i]<-filt$b%*%xh[lenh/2+(-len/2)-1+j:(j-L+1),i]
  perf_mat[i,i]<-(2*pi/length(Gamma))*abs(Gamma-trffkt[,i])^2%*%periodogram[,i]
}
# Compute time-domain MSE
mse<-apply(na.exclude((yhat-y))^2,2,mean)
file = paste("z_dfa_ar1_sym_output.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,1))
for (i in 1:3)   #i<-1
{
  ymin<-min(min(y[,i]),min(na.exclude(yhat)[,i]))
  ymax<-max(max(y[,i]),max(na.exclude(yhat)[,i]))
  ts.plot(yhat[,i],main=paste("Time-domain MSE = ",
  round(mse[i],3)," , Frequency-domain MSE = ",
  round(perf_mat[i,i],3),", a1 = ",a_vec[i],sep=""),col="blue",ylim=c(ymin,ymax))
  lines(y[,i],col="red")
  mtext("Real-time", side = 3, line = -1,at=len/2,col="blue")
  mtext("target", side = 3, line = -2,at=len/2,col="red")
}

dev.off()


###################################################
### code chunk number 80: z_dfa_ar1_sym_output.pdf
###################################################
  file = paste("z_dfa_ar1_sym_output", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Real-time filter output (blue) vs. targets (red) for a1=0.9 (top), a1=0 (middle) and a1=-0.9 (bottom)", sep = "")
  cat("\\label{z_dfa_ar1_sym_output}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 81: eco2.rnw:3176-3255
###################################################
L<-13
yhat_Lag<-array(dim=c(len,3,L/2+2))
trffkt<-array(dim=c(len/2+1,3,L/2+2))
b<-array(dim=c(L,3,L/2+2))
# Compute real-time filters for Lag=,...,L/2 and for the above three AR-processes
for (i in 1:3)
{
  periodogram[,i]<-per(x[,i],plot_T)$per
  for (Lag in 0:((L/2)+1))
  {
# Optimize filters
    filt<-dfa_ms(L,periodogram[,i],Lag,Gamma)
    trffkt[,i,Lag+1]<-filt$trffkt
    b[,i,Lag+1]<-filt$b
# Compute outputs
    for (j in L:len)
      yhat_Lag[j,i,Lag+1]<-filt$b%*%x[j:(j-L+1),i]
  }
}
omega_k<-pi*0:(len/2)/(len/2)
colo<-rainbow(L/2+2)
file = paste("z_dfa_ar1_amp_shift_Lag.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,4))
amp<-abs(trffkt)
shift<-Arg(trffkt)/omega_k
for (i in 1:3)
{
  ymin<-min(amp[,i,],na.rm=T)
  ymax<-max(amp[,i,],na.rm=T)
  plot(amp[,i,1],type="l",main=paste("Amplitude functions, a1 = ",a_vec[i],sep=""),
  axes=F,xlab="Frequency",ylab="Amplitude",col=colo[1],ylim=c(ymin,ymax))
  mtext("Lag=0", side = 3, line = -1,at=len/4,col=colo[1])
  for (j in 2:(L/2+2))
  {
    lines(amp[,i,j],col=colo[j])
    mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=len/4,col=colo[j])
  }
  axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
  axis(2)
  box()
  ymin<-min(shift[,i,],na.rm=T)
  ymax<-max(shift[,i,],na.rm=T)
  plot(shift[,i,1],type="l",main=paste("Time-Shifts, a1 = ",a_vec[i],sep=""),
  axes=F,xlab="Frequency",ylab="Shift",col=colo[1],ylim=c(ymin,ymax))
  mtext("Lag=0", side = 3, line = -1,at=len/4,col=colo[1])
  for (j in 2:(L/2+2))
  {
    lines(shift[,i,j],col=colo[j])
    mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=len/4,col=colo[j])
  }
  axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
  axis(2)
  box()
  ymin<-min(b[,i,],na.rm=T)
  ymax<-max(b[,i,],na.rm=T)
  plot(b[,i,1],col=colo[1],ylim=c(ymin,ymax),main=paste("Filter coefficients"),
  ylab="Output",xlab="lag",axes=F,typ="l")
  for (j in 2:(L/2+2))
  {
    lines(b[,i,j],col=colo[j],type="l")
    mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=L/2,col=colo[j])
  }
  axis(1,at=1:L,labels=-1+1:L)
  axis(2)
  box()
  ymin<-min(yhat_Lag[,i,],na.rm=T)
  ymax<-max(yhat_Lag[,i,],na.rm=T)
  ts.plot(yhat_Lag[,i,1],col=colo[1],ylim=c(ymin,ymax),main=paste("Output series"),
  ylab="Output")
  mtext("Lag=0", side = 3, line = -1,at=len/4,col=colo[1])
  for (j in 2:(L/2+2))
  {
    lines(yhat_Lag[,i,j],col=colo[j])
    mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=len/4,col=colo[j])
  }

}
dev.off()


###################################################
### code chunk number 82: z_dfa_ar1_amp_shift_Lag.pdf
###################################################
  file = paste("z_dfa_ar1_amp_shift_Lag.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Amplitude (left) and time-shift (right) functions as a function of Lag (rainbow) for
  a1=0.9 (top), a1=0 (middle) and a1=-0.9 (bottom)", sep = "")
  cat("\\label{z_dfa_ar1_amp_shift_Lag}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 83: eco2.rnw:3269-3329
###################################################
file = paste("z_dfa_ar1_amp_shift_Lag_0.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,2))
amp<-abs(trffkt)
shift<-Arg(trffkt)/omega_k
for (i in 2:2)
{
  ymin<-min(amp[,i,],na.rm=T)
  ymax<-max(amp[,i,],na.rm=T)
  plot(amp[,i,1],type="l",main=paste("Amplitude functions, a1 = ",a_vec[i],sep=""),
  axes=F,xlab="Frequency",ylab="Amplitude",col=colo[1],ylim=c(ymin,ymax))
  mtext("Lag=0", side = 3, line = -1,at=len/4,col=colo[1])
  for (j in 2:(L/2+2))
  {
    lines(amp[,i,j],col=colo[j])
    mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=len/4,col=colo[j])
  }
  axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
  "4pi/6","5pi/6","pi"))
  axis(2)
  box()
  ymin<-min(shift[,i,],na.rm=T)
  ymax<-max(shift[,i,],na.rm=T)
  plot(shift[,i,1],type="l",main=paste("Time-Shifts, a1 = ",a_vec[i],sep=""),
  axes=F,xlab="Frequency",ylab="Shift",col=colo[1],ylim=c(ymin,ymax))
  mtext("Lag=0", side = 3, line = -1,at=len/4,col=colo[1])
  for (j in 2:(L/2+2))
  {
    lines(shift[,i,j],col=colo[j])
    mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=len/4,col=colo[j])
  }
  axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
  "4pi/6","5pi/6","pi"))
  axis(2)
  box()
  ymin<-min(b[,i,],na.rm=T)
  ymax<-max(b[,i,],na.rm=T)
  plot(b[,i,1],col=colo[1],ylim=c(ymin,ymax),main=paste("Filter coefficients"),
  ylab="Output",xlab="lag",axes=F,typ="l")
  for (j in 2:(L/2+2))
  {
    lines(b[,i,j],col=colo[j],type="l")
    mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=L/2,col=colo[j])
  }
  axis(1,at=1:L,labels=-1+1:L)
  axis(2)
  box()

  ymin<-min(yhat_Lag[,i,],na.rm=T)
  ymax<-max(yhat_Lag[,i,],na.rm=T)
  ts.plot(yhat_Lag[,i,1],col=colo[1],ylim=c(ymin,ymax),
  main=paste("Output series"),ylab="Output")
  for (j in 2:(L/2+2))
  {
    lines(yhat_Lag[,i,j],col=colo[j])
    mtext(paste("Lag=",j-1,sep=""), side = 3, line = -j,at=len/4,col=colo[j])
  }

}
dev.off()


###################################################
### code chunk number 84: z_dfa_ar1_amp_shift_Lag_0.pdf
###################################################
  file = paste("z_dfa_ar1_amp_shift_Lag_0.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Amplitude (left) and time-shift (right) functions as a function of Lag (rainbow colors) for
  the white noise process (a1=0)", sep = "")
  cat("\\label{z_dfa_ar1_amp_shift_Lag_0}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 85: eco2.rnw:3372-3442
###################################################
# Compute real-time filters of length 3,6,12,24 and 48
Lag<-0
L_vec<-c(3,6,12,24,48)
# Compute the coefficients of the symmetric target filter
cutoff<-pi/6
# Order of approximation
ord<-999
gamma<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))
# Determine full-length and in-sample length
len<-120
lenh<-2240
set.seed(10)
anzsim<-100
perf_outsample<-matrix(nrow=anzsim,ncol=length(L_vec))
perf_insample<-perf_outsample
b_sim<-as.list(1:anzsim)
b<-as.list(1:length(L_vec))
trffkt<-array(dim=c(len/2+1,length(L_vec),anzsim))
ymin<-1:length(L_vec)
ymax<-1:length(L_vec)
# Target in frequency-domain
Gamma<-c(1,(1:(len/2))<len/12)
#
# Start simulation
for (k in 1:anzsim)
{
# Full length data
  xh<-rnorm(lenh)
# in-sample data-set: used for estimation of periodogram
  x<-xh[(lenh/2-len):(lenh/2-1)]
  y<-x
# compute output of symmetric filter
  for (j in 1:(2*len))
  {
    y[j]<-gamma%*%xh[lenh/2-len-1+(j:(j-(ord)))]+
    gamma[2:ord]%*%xh[lenh/2-len+(j:(j+ord-2))]
  }
  yhat_insample<-array(dim=c(2*len,length(L_vec)))
# Compute Periodogram based on in-sample data
  per_wn<-per(x,plot_T)$per
  for (i in 1:length(L_vec))
  {
# Optimize filters in sample
    filt<-dfa_ms(L_vec[i],per_wn,Lag,Gamma)
    trffkt[,i,k]<-filt$trffkt
    b[[i]]<-filt$b
    if (k==1)
    {
      ymin[i]<-min(b[[i]])
      ymax[i]<-max(b[[i]])
    } else
    {
      ymin[i]<-min(ymin[i],min(b[[i]]))
      ymax[i]<-max(ymax[i],max(b[[i]]))
    }
# Compute the outputs in and out-of-sample
    for (j in 1:(2*len))#j<-1
      yhat_insample[j,i]<-b[[i]]%*%xh[lenh/2-len-1+(j:(j-(L_vec[i])+1))]
  }
# Store filter coefficients
  b_sim[[k]]<-b
# MSE performances: in and out-of-sample
  perf_insample[k,]<-apply(na.exclude(y[1:len]-yhat_insample[1:len,])^2,2,mean)
  perf_outsample[k,]<-
  apply(na.exclude(y[len+1:len]-yhat_insample[len+1:len,])^2,2,mean)
}

perf_mat<-rbind(apply(perf_insample,2,mean),apply(perf_outsample,2,mean))
dimnames(perf_mat)[[2]]<-paste("L=",L_vec,sep="")
dimnames(perf_mat)[[1]]<-c("In-sample","Out-of-sample")


###################################################
### code chunk number 86: perf_mat_dfa_i_o
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(perf_mat, dec = 1,digits=rep(3,dim(perf_mat)[2]+1),
  paste("In- and out-of-sample MSE-performances of real-time trend extraction filters of length 3-48",sep=""),
  label=paste("perf_mat_dfa_i_o",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 87: eco2.rnw:3465-3480
###################################################
file = paste("z_dfa_wn_b_dist.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,2))
for (i in 2:length(L_vec))
{
  b_opt<-dfa_ms(L_vec[i],rep(1,length(Gamma)),Lag,Gamma)$b
  ts.plot(b_sim[[1]][[i]],ylim=c(ymin[i],ymax[i]),
  main=paste("Filter coefficients L = ",L_vec[i],sep=""),xlab="lag",
  ylab="coefficients")
  for (j in 2:anzsim)
    lines(b_sim[[j]][[i]])
  lines(b_opt,col="red",lwd=2)
}
dev.off()


###################################################
### code chunk number 88: z_dfa_wn_b_dist.pdf
###################################################
  file = paste("z_dfa_wn_b_dist.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Empirical distribution of filter coefficients as a function of the filter length L: DGP is white noise (a1=0)", sep = "")
  cat("\\label{z_dfa_wn_b_dist}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 89: eco2.rnw:3503-3518
###################################################
file = paste("z_dfa_wn_amp_dist.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
par(mfrow=c(2,2))
for (i in 2:length(L_vec)) #i<-2
{
  trffkt_opt<-dfa_ms(L_vec[i],rep(1,length(Gamma)),Lag,Gamma)$trffkt
  ts.plot(abs(trffkt[,i,1]),ylim=c(min(abs(trffkt[,i,])),max(abs(trffkt[,i,]))),
  main=paste("Amplitude:  L = ",
  L_vec[i],sep=""),xlab="lag",ylab="coefficients")
  for (j in 2:anzsim)
    lines(abs(trffkt[,i,j]))
    lines(abs(trffkt_opt),col="red",lwd=2)
}
dev.off()


###################################################
### code chunk number 90: z_dfa_wn_amp_dist.pdf
###################################################
  file = paste("z_dfa_wn_amp_dist.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Empirical distribution of amplitude functions as a function of the filter length L: DGP is white noise (a1=0)", sep = "")
  cat("\\label{z_dfa_wn_amp_dist}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 91: eco2.rnw:3534-3550
###################################################
file = paste("z_dfa_wn_shift_dist.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special",
width = 6, height = 6)
omega_k<-pi*0:(len/2)/(len/2)
par(mfrow=c(2,2))
for (i in 2:length(L_vec)) #i<-2
{
  trffkt_opt<-dfa_ms(L_vec[i],rep(1,length(Gamma)),Lag,Gamma)$trffkt
  ts.plot(Arg(trffkt[,i,1])/omega_k,ylim=c(min(na.exclude(Arg(trffkt[,i,])/omega_k)),max(na.exclude(Arg(trffkt[,i,])/omega_k))),
  main=paste("Shift:  L = ",L_vec[i],sep=""),xlab="lag",
  ylab="coefficients")
  for (j in 2:anzsim)
    lines(Arg(trffkt[,i,j])/omega_k)
    lines(Arg(trffkt_opt)/omega_k,col="red",lwd=2)
}
dev.off()


###################################################
### code chunk number 92: z_dfa_wn_shift_dist.pdf
###################################################
  file = paste("z_dfa_wn_shift_dist.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Empirical distribution of time-shift functions as a function of the filter length L: DGP is white noise (a1=0)", sep = "")
  cat("\\label{z_dfa_wn_shift_dist}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 93: eco2.rnw:3593-3595
###################################################
US_GDP<-read.csv(paste(path.dat,"US_GDP.csv",sep=""),header=T)
US_GDP_wp<-read.csv(paste(path.dat,"US_GDP_wp.csv",sep=""),header=T,sep=";")


###################################################
### code chunk number 94: US_GDP
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(US_GDP,
  paste("US-GDP: yearly vintages starting in Q1 2009 and ending in Q1 2013",sep=""),
  label=paste("US_GDP",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 95: eco2.rnw:3733-3750
###################################################
vintage<-array(dim=c(len,3,len))
dim(vintage)
# For each of the three AR(1)-processes We compute the vintage series
for (i in 1:3)
{
  for (j in L:len)#j<-L
  {
    vintage[(j-as.integer(L/2)):j,i,j]<-yhat_Lag[j,i,(as.integer(L/2)+1):1]
    vintage[1:(j-as.integer(L/2)-1),i,j]<-
    yhat_Lag[(as.integer(L/2)+1):(j-1),i,as.integer(L/2)+1]
  }
}
# We select the third DGP with a1=-0.9
i<-3
vintage_triangle<-vintage[,i,]
dimnames(vintage_triangle)[[2]]<-paste("Publ. ",1:len,sep="")
dimnames(vintage_triangle)[[1]]<-paste("Target ",1:len,sep="")


###################################################
### code chunk number 96: vintage_triangle
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(vintage_triangle[(len-5):len,(len-5):len], dec = 1,digits=rep(3,dim(vintage_triangle[(len-5):len,(len-5):len])[2]+1),
  paste("Last 5 vintages of finite sample filter for the AR(1)-process with a1=-0.9: columns correspond to vintages and are indexed
  by corresponding publication dates; rows correspond to revisions of estimates for a fixed historical target date",sep=""),
  label=paste("vintage_triangle",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 97: eco2.rnw:3786-3804
###################################################
colo<-rainbow(len)
file = paste("z_vintages.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,1))
for (i in 1:3)
{
  ymin<-min(vintage[,i,],na.rm=T)
  ymax<-max(vintage[,i,],na.rm=T)
  ts.plot(vintage[,i,L],col=colo[1],ylim=c(ymin,ymax),
  main=paste("Tentacle plot: vintages and full revision sequence,
  a1 = ",a_vec[i],sep=""),ylab="Vintages")
  for (j in (L+1):len)
  {
    lines(vintage[,i,j],col=colo[j])
  }
  lines(vintage[,i,len],col="red",lwd=2)
}
dev.off()


###################################################
### code chunk number 98: z_vintages.pdf
###################################################
  file = paste("z_vintages.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Tentacle plot: full historical revision sequence for
  a1=0.9 (top), a1=0 (middle) and a1=-0.9 (bottom). Final release is emphasized in bold red", sep = "")
  cat("\\label{z_vintages}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 99: eco2.rnw:3829-3854
###################################################
file = paste("z_vintages_2.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
i<-2
ymin<-min(vintage[,i,],na.rm=T)
ymax<-max(vintage[,i,],na.rm=T)
ts.plot(vintage[,i,L],col=colo[1],ylim=c(ymin,ymax),
main="Vintages: full revision sequence and final release (black)",ylab="Vintages")
for (j in (L+1):len)
{
  lines(vintage[,i,j],col=colo[j])
}
lines(vintage[,i,len],col="black",lwd=2)
i<-2
ymin<-min(vintage[,i,],na.rm=T)
ymax<-max(vintage[,i,],na.rm=T)
ts.plot(vintage[,i,L],col=colo[1],ylim=c(ymin,ymax),
main="Vintages: full revision sequence and real-time initial release (black)",
ylab="Vintages")
for (j in (L+1):len)
{
  lines(vintage[,i,j],col=colo[j])
}
lines(yhat_Lag[,i,1],col="black",lty=1)
dev.off()


###################################################
### code chunk number 100: z_vintages_2.pdf
###################################################
  file = paste("z_vintages_2.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Vintages: full historical revision sequence in the case of the white noise process", sep = "")
  cat("\\label{z_vintages_2}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 101: eco2.rnw:4148-4218
###################################################
# This function computes analytic DFA-solutions
# L is the length of the MA-filter,
# weight_func is the periodogram,
# lambda emphasizes phase artifacts in the customized criterion,
# eta emphasizes noise-suppression/smoothness
# Gamma is the transferfunction of the symmetric filter (target) and
# Lag is the lag-parameter: Lag=0 implies real-time filtering, Lag=L/2
#     implies symmetric filter
# i1 and i2 allow for filter restrictions in frequency zero
# The function returns the weights of the MA-Filter as well as its transferfunction
#
#
dfa_analytic<-function(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2)
{
  K<-length(weight_func)-1
# Define the amplitude weighting function weight_h (W(omega_k))
  omega_Gamma<-as.integer(cutoff*K/pi)
  if ((K-omega_Gamma+1)>0)
  {
    weight_h<-weight_func*(c(rep(1,omega_Gamma),(1:(K-omega_Gamma+1))^(eta)))
  } else
  {
    weight_h<-weight_func*rep(1,K+1)
  }
# First order filter restriction: assigning a `large' weight to frequency zero
  if (i1)
    weight_h[1]<-max(1.e+10,weight_h[1])

  X<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)*sqrt(weight_h)
  X_y<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)
  if (i2)
  {
# Second order restriction: time shift in frequency zero vanishes
    for (l in 2:(L-1))
    {
      X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
      cos((L-1-Lag)*pi*(0:(K))/(K))+
      sqrt(1+Gamma*lambda)*1.i*(sin((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
      sin((L-1-Lag)*pi*(0:(K))/(K))))*sqrt(weight_h))
      X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
      cos((L-1-Lag)*pi*(0:(K))/(K))+
      1.i*(sin((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*sin((L-1-Lag)*pi*(0:(K))/(K)))))
    }
    xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
  # MA-Filterweights
    b<-as.vector(solve(xtx)%*%(t(Re(X_y))%*%(Gamma*weight_h)))
# the last weight is a function of the previous ones through the second order restriction
    b<-c(b,-sum(b*(0:(length(b)-1)))/(length(b)))
  } else
  {
    for (l in 2:L)
    {
      X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))+
      sqrt(1+Gamma*lambda)*1.i*sin((l-1-Lag)*pi*(0:(K))/(K)))*sqrt(weight_h))
      X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))+
      1.i*sin((l-1-Lag)*pi*(0:(K))/(K))))
    }
    xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
  # MA-Filterweights
    b<-as.vector(solve(xtx)%*%(t(Re(X_y))%*%(Gamma*weight_h)))
  }
  # Transferfunction
  trffkt<-1:(K+1)
  trffkt[1]<-sum(b)
  for (k in 1:(K))#k<-1
  {
    trffkt[k+1]<-(b%*%exp(1.i*k*(0:(length(b)-1))*pi/(K)))
  }
  return(list(b=b,trffkt=trffkt))
}


###################################################
### code chunk number 102: eco2.rnw:4288-4309
###################################################
# Simulation of the AR(1)+cos series
len<-600
set.seed(10)
eps<-rnorm(len)
ar1<-0.9
z<-eps
for (i in 2:len)
  z[i]<-ar1*z[i-1]+eps[i]
x<-z+cos((1:len)*pi/6)
plot_T<-F

file = paste("z_seas_a.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,1))
ts.plot(x)
acf(x)
plot(per(x,plot_T)$per,type="l",axes=F,ylab="Periodogram")
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 103: z_seas_a.pdf
###################################################
  file = paste("z_seas_a.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Seasonal series (top), acf (middle) and periodogram (bottom)", sep = "")
  cat("\\label{z_seas_a}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 104: eco2.rnw:4338-4361
###################################################
# assigning the periodogram
weight_func<-per(x,plot_T)$per
K<-length(weight_func)-1
# Define target Gamma
Gamma<-rep(1,K+1)
Gamma[1+(length(Gamma)/6):(length(Gamma)/6)]<-0
# Compute amplitude of seasonal difference
omega_k<-(0:K)*pi/K
b1<-1
b12<--1
trffkt<-1+b12*complex(arg=12*omega_k)
amplitude<-abs(trffkt)

file = paste("z_seas_a_a.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,2),ylab="",xlab="")
lines(amplitude,lty=1,col="red")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
mtext("Amplitude seasonal differences", side = 3, line = -2,at=K/2,col="red")
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 105: z_seas_a_a.pdf
###################################################
  file = paste("z_seas_a_a.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Seasonal adjustment: target (blue) vs. traditional seasonal adjustment filter (black)", sep = "")
  cat("\\label{z_seas_a_a}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 106: eco2.rnw:4391-4437
###################################################
# Length of Moving-Average:
# Trick: selecting L=24 to be a multiple of 6 (the quotient in the
#    frequency pi/6 of the cosine) allows
#    the filter to damp/eliminate the component more easily!
L<-24
# Customization: lambda=eta=0 implies minimal mean-square revisions
lambda<-0
eta<-0
# cutoff frequency: this parameter is inactive if eta=0
cutoff<-pi/6
# Real-time filter
Lag<-0
# We do not impose filter restrictions: see exercises below
i1<-F
i2<-F
# Compute DFA
#
dfa_mse<-dfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2)
#
b_ms<-dfa_mse$b
trffkt_ms<-dfa_mse$trffkt

# Amplitude and time shift
amp_analytic_ms<-abs(trffkt_ms)
pha_analytic_ms<-Arg(trffkt_ms)/(pi*(0:(K))/(K))

file = paste("z_seas_a_a_at.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,max(amp_analytic_ms)),ylab="",
xlab="",main="Amplitude")
lines(amp_analytic_ms,lty=1,col="red")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
mtext("Real-time amplitude MSE-filter", side = 3, line = -2,at=K/2,col="red")
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(rep(0,K+1),type="l",axes=F,col="black",ylim=c(min(na.exclude(pha_analytic_ms)),
max(na.exclude(pha_analytic_ms))),ylab="",xlab="",main="Time-shift")
lines(pha_analytic_ms,lty=1,col="red")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
mtext("Real-time amplitude MSE-filter", side = 3, line = -2,at=K/2,col="red")
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 107: z_seas_a_a_at.pdf
###################################################
  file = paste("z_seas_a_a_at.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Seasonal adjustment: real-time MSE-filter (red) vs. target (black): amplitude (top) and
  time-shift (bottom)", sep = "")
  cat("\\label{z_seas_a_a_at}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 108: eco2.rnw:4453-4477
###################################################
# Filtering time series
xf_ms<-x
for (i in L:length(x))
{
  xf_ms[i]<-b_ms%*%x[i:(i-L+1)]
}

file = paste("z_seas_a_a_afbe.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
ts.plot(xf_ms,col="red",ylab="")
lines(x,col="blue")
mtext("Input", side = 3, line = -1,at=K,col="blue")
mtext("Real-time output", side = 3, line = -2,at=K,col="red")
plot((per(xf_ms[L:len],plot_T)$per),type="l",axes=F,col="red",
ylab="Periodograms",xlab="Frequency")
lines((per(x[L:len],plot_T)$per),lty=2,col="blue")
mtext("Input", side = 3, line = -1,at=K/2,col="blue")
mtext("Real-time output", side = 3, line = -2,at=K/2,col="red")
axis(1,at=1+0:6*(len-L)/12,labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 109: z_seas_a_a_afbe.pdf
###################################################
  file = paste("z_seas_a_a_afbe.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Filter effect: original (blue) and filtered (red) series (top graph) and
  periodograms (bottom graph)", sep = "")
  cat("\\label{z_seas_a_a_afbe}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 110: eco2.rnw:4516-4566
###################################################
weight_func_c<-weight_func
# Here We customize
weight_func_c[(length(weight_func)/6+1)]<-10^6*weight_func[length(weight_func)/6+1]

# Compute customized filter
dfa_c<-dfa_analytic(L,lambda,weight_func_c,Lag,Gamma,eta,cutoff,i1,i2)
b_c<-dfa_c$b
trffkt_c<-dfa_c$trffkt
# Amplitude and time shift
amp_analytic_c<-abs(trffkt_c)
pha_analytic_c<-Arg(trffkt_c)/(pi*(0:(K))/(K))

# Filter series
xf_c<-x
for (i in L:length(x))
{
  xf_c[i]<-b_c%*%x[i:(i-L+1)]
}

file = paste("z_seas_per_ms_c.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)

par(mfrow=c(2,1))
plot((per(xf_ms[L:len],plot_T)$per),type="l",axes=F,col="red",
ylab="Periodogram",xlab="Frequency",
main="Comparisons of periodograms of input and output signals")
lines((per(x[L:len],plot_T)$per),col="blue")
lines((per(xf_c[L:len],plot_T)$per),col="green")
mtext("Input", side = 3, line = -1,at=K/2,col="blue")
mtext("Output mean-square", side = 3, line = -2,at=K/2,col="red")
mtext("Output customized", side = 3, line = -3,at=K/2,col="green")
axis(1,at=1+0:6*(len-L)/12,labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()

hl<-10
plot((per(xf_ms[L:len],plot_T)$per)[((len-L)/12)+(-hl:hl)],type="l",
axes=F,col="red",ylab="Periodogram",xlab="Frequency",
main="Magnifying glass on pi/6",
ylim=c(0,max((per(x[L:len],plot_T)$per)[((len-L)/12)+(-hl:hl)])))
lines((per(x[L:len],plot_T)$per)[((len-L)/12)+(-hl:hl)],col="blue")
lines((per(xf_c[L:len],plot_T)$per)[((len-L)/12)+(-hl:hl)],col="green")
mtext("Input", side = 3, line = -1,at=K/2,col="blue")
mtext("Output mean-square", side = 3, line = -2,at=K/2,col="red")
mtext("Output customized", side = 3, line = -3,at=K/2,col="green")
axis(1,at=hl+2,labels="pi/6")
axis(2)
box()
dev.off()


###################################################
### code chunk number 111: z_seas_per_ms_c.pdf
###################################################
  file = paste("z_seas_per_ms_c.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Periodograms of input (blue), output MSE (red) and output customized (green): full
  bandwith (top) and vicinity of pi/6 (bottom)", sep = "")
  cat("\\label{z_seas_per_ms_c}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 112: eco2.rnw:4583-4608
###################################################
file = paste("z_seas_a_ms_c.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)

par(mfrow=c(2,1))
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,max(amp_analytic_ms)),
ylab="",xlab="",main="Amplitude")
lines(amp_analytic_ms,lty=1,col="red")
lines(amp_analytic_c,lty=1,col="green")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
mtext("Real-time amplitude MSE", side = 3, line = -2,at=K/2,col="red")
mtext("Real-time amplitude customized", side = 3, line = -3,at=K/2,col="green")
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(rep(0,K+1),type="l",axes=F,col="black",ylim=c(min(na.exclude(pha_analytic_c)),
max(na.exclude(pha_analytic_c))),ylab="",xlab="",main="Time-shift")
lines(pha_analytic_ms,lty=1,col="red")
lines(pha_analytic_c,lty=1,col="green")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
mtext("Real-time amplitude MSE", side = 3, line = -2,at=K/2,col="red")
mtext("Real-time amplitude customized", side = 3, line = -3,at=K/2,col="green")
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 113: z_seas_a_ms_c.pdf
###################################################
  file = paste("z_seas_a_ms_c.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Amplitude (top) and time-shifts (bottom) of MSE (red) and customized filters (green)", sep = "")
  cat("\\label{z_seas_a_ms_c}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 114: eco2.rnw:4642-4696
###################################################
# Length of Moving-Averages:
L_vec<-c(60,120)
# Compute analytical DFA: MSE-filters (no customization)
for (i in 1:length(L_vec))
{
  assign(paste("dfa_mse_",L_vec[i],sep=""),
  dfa_analytic(L_vec[i],lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2))
  assign(paste("b_ms_",L_vec[i],sep=""),
  get(paste("dfa_mse_",L_vec[i],sep=""))$b)
  assign(paste("trffkt_ms_",L_vec[i],sep=""),
  get(paste("dfa_mse_",L_vec[i],sep=""))$trffkt)
  assign(paste("amp_analytic_ms_",L_vec[i],sep=""),
  abs(get(paste("trffkt_ms_",L_vec[i],sep=""))))
  assign(paste("pha_analytic_ms_",L_vec[i],sep=""),
  Arg(get(paste("trffkt_ms_",L_vec[i],sep="")))/(pi*(0:(K))/(K)))
}

file = paste("z_seas_a_mse_c_mse.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,max(amp_analytic_c)),
ylab="",xlab="",main="Amplitude")
lines(amp_analytic_ms,lty=1,col="red")
lines(amp_analytic_c,lty=1,col="green")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
mtext("MSE-filter: L=24", side = 3, line = -2,at=K/2,col="red")
mtext("Customized filter: L=24", side = 3, line = -3,at=K/2,col="green")
colo<-c("violet","orange")
for (i in 1:length(L_vec))
{
  lines(get(paste("amp_analytic_ms_",L_vec[i],sep="")),col=colo[i])
  mtext(paste("MSE-Filter: L=",L_vec[i],sep=""), side = 3, line = -3-i,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
hl<-30
plot(Gamma[(K/6)+(-hl:hl)],type="l",axes=F,col="black",
ylim=c(0,max(amp_analytic_c)),ylab="",xlab="",
main="Amplitude: magnifying-glass in pi/6")
lines(amp_analytic_ms[(K/6)+(-hl:hl)],lty=1,col="red")
lines(amp_analytic_c[(K/6)+(-hl:hl)],lty=1,col="green")
mtext("Target", side = 3, line = -1,at=hl,col="black")
mtext("MSE-filter: L=24", side = 3, line = -2,at=hl,col="red")
mtext("Customized filter: L=24", side = 3, line = -3,at=hl,col="green")
for (i in 1:length(L_vec))
{
  lines(get(paste("amp_analytic_ms_",L_vec[i],sep=""))[(K/6)+(-hl:hl)],col=colo[i])
  mtext(paste("MSE-Filter: L=",L_vec[i],sep=""), side = 3, line = -3-i,at=hl,col=colo[i])
}
axis(1,at=hl+2,labels="pi/6")
axis(2)
box()
dev.off()


###################################################
### code chunk number 115: z_seas_a_mse_c_mse.pdf
###################################################
  file = paste("z_seas_a_mse_c_mse.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Amplitude functions full bandwith (top) and magnifying glass in pi/6 (bottom): MSE L=24 (red),
  customized L=24 (green), MSE L=60 (violet) and MSE L=120 (orange)", sep = "")
  cat("\\label{z_seas_a_mse_c_mse}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 116: eco2.rnw:4719-4766
###################################################
# Filter series
xf<-matrix(ncol=length(L_vec),nrow=length(x))
for (j in 1:length(L_vec))
{
  for (i in L_vec[j]:length(x))
  {
    xf[i,j]<-get(paste("b_ms_",L_vec[j],sep=""))%*%x[i:(i-L_vec[j]+1)]
  }
}

file = paste("z_seas_per_ms_c_ms.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)

par(mfrow=c(2,1))
plot((per(xf_ms[max(L_vec):len],plot_T)$per),type="l",axes=F,col="red",
ylab="Periodogram",xlab="Frequency",
main="Comparisons of periodograms of input and output signals")
lines((per(xf_c[max(L_vec):len],plot_T)$per),col="green")
mtext("MSE-Filter: L=24", side = 3, line = -3,at=K/2,col="red")
mtext("Output customized", side = 3, line = -2,at=K/2,col="green")
for (i in 1:length(L_vec))
{
  lines((per(xf[max(L_vec):len,i],plot_T)$per),col=colo[j])
  mtext(paste("MSE-Filter: L=",L_vec[i],sep=""), side = 3, line = -3-i,
  at=K/2,col=colo[i])
}
axis(1,at=1+0:6*(len-max(L_vec))/12,labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()

hl<-30
plot((per(xf_ms[max(L_vec):len],plot_T)$per)[(len-max(L_vec))/12+(-hl:hl)],type="l",axes=F,
col="red",ylab="Periodogram",xlab="Frequency",
main="Magnifying glass on pi/6")
lines((per(xf_c[max(L_vec):len],plot_T)$per)[(len-max(L_vec))/12+(-hl:hl)],col="green")
mtext("MSE-Filter: L=24", side = 3, line = -3,at=hl,col="red")
mtext("Output customized", side = 3, line = -2,at=hl,col="green")
for (i in 1:length(L_vec))   #i<-2
{
  lines((per(xf[max(L_vec):len,i],plot_T)$per)[(len-max(L_vec))/12+(-hl:hl)],col=colo[i])
  mtext(paste("MSE-Filter: L=",L_vec[i],sep=""), side = 3, line = -3-i,at=hl,col=colo[i])
}
axis(1,at=hl+2,labels="pi/6")
axis(2)
box()
dev.off()


###################################################
### code chunk number 117: z_seas_per_ms_c_ms.pdf
###################################################
  file = paste("z_seas_per_ms_c_ms.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Periodograms of  MSE L=24 (red), output customized L=24 (green),
  output MSE L=60 (violet) and output MSE L=120 (orange) : full
  bandwith (top) and vicinity of pi/6 (bottom)", sep = "")
  cat("\\label{z_seas_per_ms_c_ms}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 118: eco2.rnw:4782-4808
###################################################
file = paste("z_seas_per_ms_120_i.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)

par(mfrow=c(2,1))
plot((per(x[max(L_vec):len],plot_T)$per),type="l",axes=F,col="blue",
ylab="Periodogram",xlab="Frequency",
main="Comparisons of periodograms of input and output signals")
lines(per(xf[max(L_vec):len,i],plot_T)$per,col="orange")
mtext("Input", side = 3, line = -1,at=K/2,col="blue")
mtext("Output MSE 120", side = 3, line = -2,at=K/2,col="orange")
axis(1,at=1+0:6*(len-max(L_vec))/12,labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()

hl<-30
plot((per(x[max(L_vec):len],plot_T)$per)[(len-max(L_vec))/12+(-hl:hl)],
type="l",axes=F,col="blue",ylab="Periodogram",xlab="Frequency",
main="Comparisons of periodograms of input and output signals")
lines((per(xf[max(L_vec):len,i],plot_T)$per)[(len-max(L_vec))/12+(-hl:hl)],col="orange")
mtext("Input", side = 3, line = -1,at=hl,col="blue")
mtext("Output MSE 120", side = 3, line = -2,at=hl,col="orange")
axis(1,at=hl+2,labels="pi/6")
axis(2)
box()
dev.off()


###################################################
### code chunk number 119: z_seas_per_ms_120_i.pdf
###################################################
  file = paste("z_seas_per_ms_120_i.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Periodograms of input (blue) and ouput MSE L=120 (orange): full
  bandwith (top) and vicinity of pi/6 (bottom)", sep = "")
  cat("\\label{z_seas_per_ms_120_i}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 120: eco2.rnw:4845-4917
###################################################
Lag<-0
anzsim<-100
len1<-360
len_sim<-180
K<-len_sim/2
set.seed(10)
L_vec<-c(24,60,120)
# Define target Gamma
Gamma_sim<-rep(1,K+1)
Gamma_sim[1+(length(Gamma_sim)/6):(length(Gamma_sim)/6)]<-0
b_sim<-as.list(1:anzsim)
b_sim<-as.list(1:anzsim)
amp_sim<-b_sim
shift_sim<-b_sim
b<-as.list(1:length(L_vec))
amp<-matrix(ncol=length(L_vec),nrow=K+1)
shift<-amp
mse_in_sample<-matrix(ncol=length(L_vec),nrow=anzsim)
mse_out_sample<-mse_in_sample
ymin<-rep(10^6,length(L_vec))
ymax<--ymin
# Start simulation
for (i in 1:simanz)
{
  z1<-arima.sim(list(ar=0.9),n=len1+max(L_vec))
  x1<-z1+cos((1:(len1+max(L_vec)))*pi/6)
  z_sim<-z1[max(L_vec)-1+1:len_sim]
  x_sim<-x1[max(L_vec)-1+1:len_sim]
  weight_func_sim<-per(x_sim,plot_T)$per
# Compute MSE-filters
  for (j in 1:length(L_vec))
  {
    dfa<-dfa_analytic(L_vec[j],lambda,weight_func_sim,Lag,Gamma_sim,eta,cutoff,i1,i2)
    b[[j]]<-dfa$b
    ymin[j]<-min(ymin[j],min(b[[j]]))
    ymax[j]<-max(ymax[j],max(b[[j]]))
    amp[,j]<-abs(dfa$trffkt)
    shift[,j]<-Arg(dfa$trffkt)/((0:K)*pi/K)
    xf<-x1
    for (k in 1:len1)
    {
      xf[max(L_vec)-1+k]<-b[[j]]%*%x1[max(L_vec)-1+k:(k-L_vec[j]+1)]
    }
    
    mse_in_sample[i,j]<-mean((z1-xf)[max(L_vec)-1+1:len_sim]^2)
    mse_out_sample[i,j]<-mean((z1-xf)[max(L_vec)-1+(len_sim+1):len1]^2)
  }
  b_sim[[i]]<-b
  amp_sim[[i]]<-amp
  shift_sim[[i]]<-shift
}
perf_mat<-rbind(apply(mse_in_sample,2,mean),apply(mse_out_sample,2,mean))
dimnames(perf_mat)[[2]]<-paste("L=",L_vec,sep="")
dimnames(perf_mat)[[1]]<-c("In-sample","Out-of-sample")

file = paste("z_idfa_dist_SA_coef.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,1))
for (j in 1:length(L_vec))
{
  b<-b_sim[[1]][[j]]
  ts.plot(b_sim[[1]][[j]],ylim=c(ymin[j],ymax[j]),
  main=paste("Empirical distribution of coefficients: L = ",L_vec[j],sep=""),
  xlab="lag",ylab="coefficients")
  for (i in 2:anzsim)
  {
    lines(b_sim[[i]][[j]])
    b<-b+b_sim[[i]][[j]]
  }
  lines(b/anzsim,col="red",lwd=2)
}
dev.off()


###################################################
### code chunk number 121: z_idfa_dist_SA_coef.pdf
###################################################
  file = paste("z_idfa_dist_SA_coef.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Empirical distribution of real-time filter coefficients of SA-filters of length 24, 60 and 120:
  mean of distribution highlighted in red", sep = "")
  cat("\\label{z_idfa_dist_SA_coef}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 122: eco2.rnw:4930-4948
###################################################
file = paste("z_idfa_dist_SA_amp.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,1))

for (j in 1:length(L_vec))
{
  amp<-amp_sim[[1]][,j]
  ts.plot(amp_sim[[1]][,j],ylim=c(0,1.8),
  main=paste("Empirical distribution of amplitude: L = ",L_vec[j],sep=""),
  xlab="lag",ylab="coefficients")
  for (i in 2:anzsim)
  {
    lines(amp_sim[[i]][,j])
    amp<-amp+amp_sim[[i]][,j]
  }
  lines(amp/anzsim,col="red",lwd=2)
}
dev.off()


###################################################
### code chunk number 123: z_idfa_dist_SA_amp.pdf
###################################################
  file = paste("z_idfa_dist_SA_amp.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Empirical distribution of real-time amplitude functions of SA-filters of length 24, 60 and 120: mean of
  distribution is highlighted in red", sep = "")
  cat("\\label{z_idfa_dist_SA_amp}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 124: eco2.rnw:4961-4979
###################################################
file = paste("z_idfa_dist_SA_shift.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,1))

for (j in 1:length(L_vec))
{
  shift<-shift_sim[[1]][,j]
  ts.plot(shift_sim[[1]][,j],ylim=c(-2,2),
  main=paste("Empirical distribution of time-shift: L = ",L_vec[j],sep=""),
  xlab="lag",ylab="coefficients")
  for (i in 2:anzsim)
  {
    lines(shift_sim[[i]][,j])
    shift<-shift+shift_sim[[i]][,j]
  }
  lines(shift/anzsim,col="red",lwd=2)
}
dev.off()


###################################################
### code chunk number 125: z_idfa_dist_SA_shift.pdf
###################################################
  file = paste("z_idfa_dist_SA_shift.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Empirical distribution of time-shift of SA-filters of length 24, 60 and 120: mean of
  distribution is highlighted in red", sep = "")
  cat("\\label{z_idfa_dist_SA_shift}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 126: perf_mat_idfa_i_o
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(perf_mat, dec = 1,digits=rep(3,dim(perf_mat)[2]+1),
  paste("In- and out-of-sample MSE-performances of real-time SA-filters of length 24, 60 and 120",sep=""),
  label=paste("perf_mat_idfa_i_o",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 127: eco2.rnw:5054-5071
###################################################
K<-length(weight_func)-1
L<-121
yhat_Lag<-array(dim=c(len,L))
trffkt<-array(dim=c(len/2+1,L))
b<-array(dim=c(L,L))
# Compute real-time filters for Lag=,...,L/2 and for the above three AR-processes
weight_func<-per(x,plot_T)$per
for (Lag in 0:(L-1))
{
# Optimize filters
  filt<-dfa_ms(L,weight_func,Lag,Gamma)
  trffkt[,Lag+1]<-filt$trffkt
  b[,Lag+1]<-filt$b
# Compute outputs
  for (j in L:len)
    yhat_Lag[j,Lag+1]<-filt$b%*%x[j:(j-L+1)]
}


###################################################
### code chunk number 128: eco2.rnw:5077-5134
###################################################
omega_k<-pi*0:(len/2)/(len/2)
colo<-rainbow(L)
file = paste("z_idfa_SA_Lag02L.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,1))
amp<-abs(trffkt)
shift<-Arg(trffkt)/omega_k
ymin<-min(amp,na.rm=T)
ymax<-max(amp,na.rm=T)
plot(amp[,1],type="l",main="Amplitude functions",
  axes=F,xlab="Frequency",ylab="Amplitude",col="blue",ylim=c(ymin,ymax))
  mtext("Lag=0", side = 3, line = -1,at=len/4,col="blue")
select_i<-c(as.integer(L/2+1),L)
for (j in 1:length(select_i))
{
  lines(amp[,select_i[j]],col=colo[select_i[j]])
  mtext(paste("Lag=",select_i[j]-1,sep=""), side = 3, line = -j-1,
  at=len/4,col=colo[select_i[j]])
}
axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
ymin<-min(shift,na.rm=T)
ymax<-max(shift,na.rm=T)
plot(shift[,1],type="l",main="Time-Shifts",
  axes=F,xlab="Frequency",ylab="Shift",col="blue",ylim=c(-1,ymax))
  mtext("Lag=0", side = 3, line = -1,at=len/4,col="blue")
for (j in 1:length(select_i))
{
# We here recompute the shift because the function Arg has troubles
#   when dealing with the periodicity of the phase
  shift_s<-Arg(trffkt[,select_i[j]]*exp(-1.i*select_i[j]*omega_k))/omega_k
  lines(shift_s+select_i[j],col=colo[select_i[j]])
  mtext(paste("Lag=",select_i[j]-1,sep=""), side = 3,
  line = -j-1,at=len/4,col=colo[select_i[j]])
}
axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
ymin<-min(b,na.rm=T)
ymax<-max(b,na.rm=T)
plot(b[,1],type="l",col="blue",ylim=c(ymin,ymax),main="Filter coefficients",
ylab="Output",xlab="lag",
axes=F)
  mtext("Lag=0", side = 3, line = -1,at=L/2,col="blue")
for (j in 1:length(select_i))
{
  lines(b[,select_i[j]],col=colo[select_i[j]])
  mtext(paste("Lag=",select_i[j]-1,sep=""), side = 3, line = -j-1,
  at=L/2,col=colo[select_i[j]])
}
axis(1,at=20*(1:(L/20)),labels=20*(1:(L/20))-1)
axis(2)
box()
dev.off()


###################################################
### code chunk number 129: z_idfa_SA_Lag02L.pdf
###################################################
  file = paste("z_idfa_SA_Lag02L.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Amplitude (top), time-shift (middle) and coefficients (bottom) of filters with Lags 0 (blue), 60 (cyan) and 120 (red)", sep = "")
  cat("\\label{z_idfa_SA_Lag02L}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 130: eco2.rnw:5165-5180
###################################################
vintage<-array(dim=c(len,len))
dim(vintage)
for (j in L:len)#j<-L
{
# first L/2 releases
  vintage[(j-as.integer(L/2)):j,j]<-yhat_Lag[j,(as.integer(L/2)+1):1]
# symmetric filter
  vintage[1:(j-as.integer(L/2)-1),j]<-yhat_Lag[(as.integer(L/2)+1):(j-1),
  as.integer(L/2)+1]
# last L/2 releases
  vintage[1:as.integer(L/2),j]<-yhat_Lag[L,L:(as.integer(L/2)+2)]
}
vintage_triangle<-vintage
dimnames(vintage_triangle)[[2]]<-paste("Publ. ",1:len,sep="")
dimnames(vintage_triangle)[[1]]<-paste("Target ",1:len,sep="")


###################################################
### code chunk number 131: eco2.rnw:5183-5209
###################################################
colo<-rainbow(len)
file = paste("z_vintages_sa.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
ymin<-min(vintage,na.rm=T)
ymax<-max(vintage,na.rm=T)
ts.plot(vintage[,L],col=colo[1],ylim=c(ymin,ymax),
main="Tentacle plot: vintages and full revision sequence",ylab="Vintages")
for (j in (L+1):len)
{
  lines(vintage[,j],col=colo[j])
}
abline(v=c(L,len-L+1),col="black",lty=2)
mtext("Start of symmetric (final) filter", side = 3, line = -1,at=L,col="black")
mtext("End of symmetric (final) filter", side = 3, line = -1,at=len-L+1,col="black")
ymin<-min(vintage,na.rm=T)
ymax<-max(vintage,na.rm=T)
ts.plot(yhat_Lag[,1],col="blue",ylim=c(ymin,ymax),
  main="Initial release (blue), last vintage (red) and true signal (green)",
  ylab="Vintages")
lines(vintage[,len],col="red")
lines(z,col="green")
abline(v=c(L,len-L+1),col="black",lty=2)
mtext("Start of symmetric (final) filter", side = 3, line = -1,at=L,col="black")
mtext("End of symmetric (final) filter", side = 3, line = -1,at=len-L+1,col="black")
dev.off()


###################################################
### code chunk number 132: z_vintages_sa.pdf
###################################################
  file = paste("z_vintages_sa.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Tentacle plot (top) and initial release vs. last vintage and true signal (bottom)", sep = "")
  cat("\\label{z_vintages_sa}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 133: eco2.rnw:5275-5303
###################################################
# Simulation of the `complex' seasonal series
len<-120
len1<-1000
set.seed(10)
ar1<-0.9
z<-arima.sim(list(ar=ar1),n=len1)
eps<-rnorm(3000)
yh<-eps
for (i in 5:length(eps))
  yh[i]<-yh[i-4]+eps[i]
y<-scale(yh[(length(eps)-len1+1):length(eps)])
a_2<--0.995^2
a_1<-2*cos(4*pi/6)*sqrt(-a_2)
w<-scale(arima.sim(list(ar=c(a_1,a_2)),n=len1))
xh<-z+cos((1:len1)*pi/6)+y+w
x<-xh[1:len]
plot_T<-F

file = paste("z_seas_a_g.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(3,1))
ts.plot(x)
acf(x)
plot(per(x,plot_T)$per,type="l",axes=F,ylab="Periodogram")
axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 134: z_seas_a_g.pdf
###################################################
  file = paste("z_seas_a_g.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Seasonal series (top), acf (middle) and periodogram (bottom)", sep = "")
  cat("\\label{z_seas_a_g}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 135: eco2.rnw:5338-5422
###################################################
Gamma_interface<-function(list_freq_seasonal,list_freq_calendar,spect)
{
  K<-length(spect)-1
# determine asymmetrical width of seasonal peaks
#    (calendar peaks are assumed to be line spectra)
# Two rules are used:
#   The monotonic decay on both sides of the peak (as long as its decaying
#     it's still part of the peak)
#   The absolute value of the peak (which must be larger than
#     med_times the medians on both sides)
  med_times<-4
  spec_inflate<-4
  med_spec_u<-1:(length(list_freq_seasonal)-1)
  med_spec_l<-1:(length(list_freq_seasonal)-1)
  len_u<-list_freq_seasonal
  len_l<-len_u

  for (i in 1:(length(list_freq_seasonal)))#i<-6
  {
    med_spec_l[i]<-2*ifelse(i==1,NA,median(spect[(list_freq_seasonal[i]-1):
    (list_freq_seasonal[i-1]+1)]))
    med_spec_u[i]<-2*ifelse(i==length(list_freq_seasonal),
    NA,median(spect[(list_freq_seasonal[i]+1):(list_freq_seasonal[i+1]-1)]))
    if (i<length(list_freq_seasonal))
    {
# median rule on the right half of the peak (peaks can be asymmetric)
      med_u<-which(!(spect[list_freq_seasonal[i]:list_freq_seasonal[i+1]]>
      med_times*med_spec_u[i]))[1]-1
      med_u<-ifelse(length(med_u)==0,list_freq_seasonal[i+1]-list_freq_seasonal[i],med_u)
# monotonic decay rule on the right half of the peak
      diff_u<-which(!diff(spect[list_freq_seasonal[i]:list_freq_seasonal[i+1]])<0)[1]-1
      diff_u<-ifelse(length(diff_u)==0,list_freq_seasonal[i+1]-list_freq_seasonal[i],diff_u)
      len_u[i]<-min(diff_u,med_u)
    } else
    {
      len_u[i]<-0
    }
    if (i>1)
    {
# same as above but on the left side of the peak: median and monotonic decay rules
      med_l<-which(!(spect[list_freq_seasonal[i]:list_freq_seasonal[i-1]]>
      med_times*med_spec_l[i]))[1]-1
      med_l<-ifelse(length(med_l)==0,list_freq_seasonal[i]-list_freq_seasonal[i-1],med_l)
      diff_l<-which(!diff(spect[list_freq_seasonal[i]:list_freq_seasonal[i-1]])<0)[1]-1
      diff_l<-ifelse(length(diff_l)==0,list_freq_seasonal[i]-list_freq_seasonal[i-1],diff_l)
      len_l[i]<-min(diff_l,med_l)
    } else
    {
      len_l[i]<-0
    }
  #  which(diff(which(spect[list_freq[i]:list_freq[i-1]]>med_spec[i-1]))>1)[1]
  }

# Formatting Gamma based on location, width and height of peaks as determined above
  Gamma<-rep(1,K+1)
  for (i in 1:(length(list_freq_seasonal))) #i<-6
  {
  # If the width of the peak is not null (location problem: where do I have peaks)
    if (len_l[i]+len_u[i]>0)
    {
# determining the width
      width_peak<-(list_freq_seasonal[i]-max(0,len_l[i]-1)):
      (list_freq_seasonal[i]+max(0,len_u[i]-1))
# determining the height (deepness of target trough/dip)
      Gamma[width_peak]<-sqrt(min(na.exclude(c(med_spec_l[i],med_spec_u[i]))))/
      (spec_inflate*sqrt(spect[width_peak]))
# We dot not allow for amplitudes larger than one
      Gamma[which(Gamma>1)]<-1
# We drill a zero in the seasonal frequency (if it is loaded by a peak)
      Gamma[list_freq_seasonal[i]]<-0
    }
  }
  if (length(list_freq_calendar)>0)
    Gamma[list_freq_calendar]<-0

#  plot_T<-F
#  plot(Gamma,type="l",axes=F,ylim=c(0,1))
#  lines(per(x_data,plot_T)$per/max(per(x_data,plot_T)$per),col="red")
#  axis(1,at=1+(0:6)*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
#  axis(2)
#  box()

  return(list(Gamma=Gamma))
}


###################################################
### code chunk number 136: eco2.rnw:5447-5465
###################################################
weight_func<-per(x,plot_T)$per
spect<-weight_func
K<-length(weight_func)-1
list_freq_seasonal<-round(1+K*c((1:6)/6))
list_freq_calendar<-NULL#round(1+K*(4/6+0.1/pi))

Gamma<-Gamma_interface(list_freq_seasonal,list_freq_calendar,weight_func)$Gamma

file = paste("z_seas_agi.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
plot_T<-F
plot(Gamma,type="l",axes=F,ylim=c(0,1),
main="Periodogram of input (red) and automatic target (black)")
lines(per(x,plot_T)$per/max(per(x,plot_T)$per),col="red")
axis(1,at=1+(0:6)*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 137: z_seas_agi.pdf
###################################################
  file = paste("z_seas_agi.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Automatic Gamma-interface", sep = "")
  cat("\\label{z_seas_agi}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 138: eco2.rnw:5483-5552
###################################################
# Length of Moving-Average:
# Trick: selecting L=24 to be a multiple of 6 (the quotient in the
#   frequency pi/6 of the cosine) allows
#   the filter to damp/eliminate the component more easily!
L_vec<-c(24,len/2)
b_list<-as.list(1:length(L_vec))
xf<-matrix(nrow=len1,ncol=length(L_vec))
trffkt<-matrix(nrow=K+1,ncol=length(L_vec))
amp<-trffkt
shift<-trffkt
# Customization: lambda=eta=0 implies minimal mean-square revisions
lambda<-0
eta<-0
# cutoff frequency: this parameter is inactive if eta=0
cutoff<-pi/6
# Lag: Lag=0: real-time filter;
Lag<-0
# We do not impose filter restrictions: see exercises below
i1<-F
i2<-F
# Compute analytical DFA
for (i in 1:length(L_vec))
{
  dfa_mse<-dfa_analytic(L_vec[i],lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2)
  b_list[[i]]<-dfa_mse$b
  trffkt[,i]<-dfa_mse$trffkt
  amp[,i]<-abs(trffkt[,i])
  shift[,i]<-Arg(trffkt[,i])/(pi*(0:(K))/(K))
  for (j in L_vec[i]:len1)
  {
    xf[j,i]<-b_list[[i]]%*%xh[j:(j-L_vec[i]+1)]
  }
}
# Final Symmetric filter
L<-61
Lag<-(L-1)/2
dfa_symmetric<-dfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2)
b_symmetric<-dfa_symmetric$b
xf_symmetric<-rep(NA,len1)
for (j in ((L-1)/2+1):(len1-(L-1)/2))
{
  xf_symmetric[j]<-b_symmetric%*%xh[(j+(L-1)/2):(j-(L-1)/2)]
}
colo<-rainbow(2*length(L_vec))
file = paste("z_seas_a_a_asag.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,max(amp)),ylab="",xlab="",main="Amplitude")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:length(L_vec))
{
  lines(amp[,i],lty=1,col=colo[i])
  mtext(paste("Amplitude: L=",L_vec[i],sep=""), side = 3, line = -1-i,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(rep(0,K+1),type="l",axes=F,col="black",ylim=c(min(na.exclude(shift)),
max(na.exclude(shift))),ylab="",xlab="",main="Time-shift")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:length(L_vec))
{
  lines(shift[,i],lty=1,col=colo[i])
  mtext(paste("Shift: L=",L_vec[i],sep=""), side = 3, line = -1-i,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 139: z_seas_a_a_asag.pdf
###################################################
  file = paste("z_seas_a_a_asag.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Real-time SA-filters: amplitude (top) and shift (bottom)", sep = "")
  cat("\\label{z_seas_a_a_asag}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 140: eco2.rnw:5574-5590
###################################################
file = paste("z_seas_per_oos_sa.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
plot(per(xh[max(L_vec):len1],plot_T)$per,type="l",axes=F,col="blue",
ylab="Periodogram",xlab="Frequency",
main="Comparisons of periodograms of input and output signals")
mtext("Input", side = 3, line = -1,at=len1/4,col="blue")
for (i in 1:length(L_vec))
{
  lines(per(xf[max(L_vec):len1,i],plot_T)$per,col=colo[i])
  mtext(paste("Output: L=",L_vec[i],sep=""), side = 3, line = -1-i,at=len1/4,col=colo[i])
}
axis(1,at=1+0:6*(len1-max(L_vec))/12,labels=c("0","pi/6","2pi/6","3pi/6",
"4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 141: z_seas_per_oos_sa.pdf
###################################################
  file = paste("z_seas_per_oos_sa.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Periodograms of input and real-time outputs: out-of-sample", sep = "")
  cat("\\label{z_seas_per_oos_sa}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 142: eco2.rnw:5615-5653
###################################################
file = paste("z_sym_rt_sa_oos.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
mplot<-cbind(xf_symmetric,xf)
ymin<-min(mplot,na.rm=T)
ymax<-max(mplot,na.rm=T)
ts.plot(mplot[,3],col="blue",ylim=c(ymin,ymax),
main="Final and real-time filter outputs (out-of-sample)",ylab="Vintages")
mtext("Final (symmetric) filter",col="blue", side = 3, line = -1,at=dim(mplot)[1]/2)
for (j in 1:length(L_vec))
{
  lines(mplot[,j],col=colo[j])
  mtext(paste("Real-time filter: L=",L_vec[j],sep=""),col=colo[j], side = 3,
  line = -1-j,at=dim(mplot)[1]/2)
}
abline(v=c(L,len1-L+1),col="black",lty=2)
mtext("Start of symmetric (final) filter", side = 3, line = -1,at=L,col="black")
mtext("End of symmetric (final) filter", side = 3, line = -1,
at=dim(mplot)[1]-L+1,col="black")
anf<-400
enf<-600
mplot<-cbind(xf_symmetric,xf)[400:600,]
ymin<-min(mplot,na.rm=T)
ymax<-max(mplot,na.rm=T)
plot(mplot[,3],type="l",col="blue",ylim=c(ymin,ymax),
main="Final and real-time filter outputs (out-of-sample)",
ylab="Outputs",xlab="", axes=F)
mtext("Final (symmetric) filter",col="blue", side = 3, line = -1,at=dim(mplot)[1]/2)
for (j in 1:length(L_vec))
{
  lines(mplot[,j],col=colo[j])
  mtext(paste("Real-time filter: L=",L_vec[j],sep=""),col=colo[j], side = 3,
  line = -1-j,at=dim(mplot)[1]/2)
}
axis(1,at=c(1,as.integer(((enf-anf)/4)*1:4)),labels=anf+c(0,as.integer(((enf-anf)/4)*1:4)))
axis(2)
box()
dev.off()


###################################################
### code chunk number 143: z_sym_rt_sa_oos.pdf
###################################################
  file = paste("z_sym_rt_sa_oos.pdf", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Comparison of real-time (green: L=60; red: L=24) and final (blue) out-of-sample outputs", sep = "")
  cat("\\label{z_sym_rt_sa_oos}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 144: eco2.rnw:5822-5982
###################################################
# This function computes:
# 1. real-time and symmetric filters,
# 2. in- and out-of-sample performances: Selectivity, Mean-Shift, Curvature
#     and Peak-Correlation
# 3. the ATS error components
#
#-----------
#
Performance_func<-function(lambda_vec,eta_vec,len1,len,x1,weight_func,Gamma,
cutoff,L,L_sym,a1,mba)
{
  K<-length(weight_func)-1
  omega_k<-(0:K)*pi/K
  #
  amp<-matrix(ncol=length(lambda_vec)+1,nrow=K+1)
  shift<-amp
  selectivity_insamp<-rep(NA,dim(amp)[2])
  mean_shift_insamp<-rep(NA,dim(amp)[2])
  curvature_insamp<-rep(NA,dim(amp)[2])
  selectivity_outsamp<-rep(NA,dim(amp)[2])
  mean_shift_outsamp<-rep(NA,dim(amp)[2])
  curvature_outsamp<-rep(NA,dim(amp)[2])

  b<-matrix(ncol=dim(amp)[2],nrow=L)
  # Determine passband
  omega_Gamma<-as.integer(cutoff*K/pi)+1
  passband<-1:omega_Gamma
  stopband<-(omega_Gamma+1):(K+1)
  Accuracy<-rep(NA,dim(amp)[2])
  Smoothness<-Accuracy
  Timeliness<-Accuracy
  Residual<-Accuracy
  lag_cor<-0:10
  peak_cor_insamp<-matrix(ncol=length(lag_cor),nrow=dim(amp)[2])
  peak_cor_outsamp<-matrix(ncol=length(lag_cor),nrow=dim(amp)[2])
  lag_max_insamp<-rep(NA,dim(amp)[2])
  lag_max_outsamp<-rep(NA,dim(amp)[2])
  xf<-matrix(nrow=len1,ncol=dim(amp)[2])
# Final symmetric MSE-filter
  Lag_sym<-(L_sym-1)/2
  lambda<-0
  eta<-0
  filt_sym<-dfa_analytic(L_sym,lambda,weight_func,Lag_sym,Gamma,eta,cutoff,i1,i2)
  trffkt_sym<-filt_sym$trffkt
  amp_sym<-abs(trffkt_sym)
  shift_sym<-Arg(trffkt_sym)/omega_k
  b_sym<-filt_sym$b
# Compute outputs on long sample
  xf_sym<-rep(NA,len1)
  for (j in (1+(L_sym-1)/2):(len1-(L_sym-1)/2))
    xf_sym[j]<-b_sym%*%x1[(L_sym-1)/2+j:(j-L_sym+1)]
  weight_func_best<-rep(1,K+1)/abs(1-a1*exp(1.i*omega_k))^2               #a1<-0
#
# Estimation loop for all eta/lambda values
  for (i in 0:(dim(amp)[2]-1))         #i<-1
  {
    if (i==0)
    {
      if (mba)
      {
        dfa<-dfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2)
      } else
      {
        dfa<-dfa_analytic(L,lambda,weight_func_best,Lag,Gamma,eta,cutoff,i1,i2)
      }
    } else
    {
      dfa<-dfa_analytic(L,lambda_vec[i],weight_func,Lag,Gamma,eta_vec[i],cutoff,i1,i2)
    }

    b[,i+1]<-dfa$b                   #
    amp[,i+1]<-abs(dfa$trffkt)
# Selectivity: immunized against scaling effects
    selectivity_insamp[i+1]<-(amp[1:(K*cutoff/pi),i+1]^2%*%weight_func[1:(K*cutoff/pi)])/
    (amp[(K*cutoff/pi+1):(K+1),i+1]^2%*%(weight_func[(K*cutoff/pi+1):(K+1)]*
    (1+(0:(K-K*cutoff/pi)))^2))
    shift[1+1:K,i+1]<-Arg(dfa$trffkt)[1+1:K]/((pi*1:K)/K)
# Shift in frequency zero
    shift[1,i+1]<-b[,i+1]%*%(0:(L-1))/amp[1,i+1]
# Mean-shift statistic: immunized against scaling effects
    mean_shift_insamp[i+1]<-mean(shift[1:(K*cutoff/pi),i+1])
    amp_error<-((Gamma-amp[,i+1]))^2*weight_func
    shift_error<-4*Gamma*amp[,i+1]*sin(shift[,i+1]*omega_k/2)^2*weight_func
    Accuracy[i+1]<-sum(amp_error[passband])
    Smoothness[i+1]<-sum(amp_error[stopband])
    Timeliness[i+1]<-sum(shift_error[passband])
    Residual[i+1]<-sum(shift_error[stopband])
    yhat<-rep(NA,len1)
    for (j in L:len1)
      yhat[j]<-b[,i+1]%*%x1[j:(j-L+1)]
    xf[,i+1]<-yhat
# Relative Curvature (immunized against scaling effects)
    curvature_insamp[i+1]<-mean(diff(yhat[1:len],diff=2)^2,na.rm=T)/
    var(yhat[1:len],na.rm=T)
    curvature_outsamp[i+1]<-mean(diff(yhat[(len+1):len1],diff=2)^2,na.rm=T)/
    var(yhat[(len+1):len1],na.rm=T)
    for (j in lag_cor)   #j<-0
    {
      if (!mba)
         peak_cor_insamp[i+1,j+1]<-cor(xf_sym[(1+(L_sym-1)/2):len],
         yhat[j+((L_sym-1)/2+1):len])
       peak_cor_outsamp[i+1,j+1]<-cor(xf_sym[(len+1):(len1-L_sym)],
       yhat[j+(len+1):(len1-L_sym)])
    }
    if (!mba)
      lag_max_insamp[i+1]<-which(peak_cor_insamp[i+1,]==max(peak_cor_insamp[i+1,]))-1
    lag_max_outsamp[i+1]<-which(peak_cor_outsamp[i+1,]==max(peak_cor_outsamp[i+1,]))-1
  }
# Rownames
  if (!mba)
  {
    dim_names<-c("Theoretical best MSE",paste("Lambda=",lambda_vec,", eta=",eta_vec,sep=""))
    if (prod(lambda_vec*eta_vec)==0)
    {
      dim_names[1+which(lambda_vec==0&eta_vec==0)]<-"DFA-MSE"
    }
  } else
  {
    dim_names<-c("MBA-MSE",paste("Lambda=",lambda_vec,", eta=",eta_vec,sep=""))
  }

  dimnames(xf)[[2]]<-dim_names
  dimnames(amp)[[2]]<-dim_names
  dimnames(shift)[[2]]<-dim_names
  dimnames(b)[[2]]<-dim_names
#
# We collect all frequency-domain performance statistics in corresponding tables
  if (!mba)
  {
    amp_shift_mat_insamp<-cbind(selectivity_insamp,mean_shift_insamp,
    curvature_insamp,lag_max_insamp)
    dimnames(amp_shift_mat_insamp)[[2]]<-c("Selectivity","Mean-shift",
    "Curvature","Peak-Correlation lag")
    dimnames(amp_shift_mat_insamp)[[1]]<-dim_names
  } else
  {
    amp_shift_mat_insamp<-cbind(selectivity_insamp,mean_shift_insamp)
    dimnames(amp_shift_mat_insamp)[[2]]<-c("Selectivity","Mean-shift")
    dimnames(amp_shift_mat_insamp)[[1]]<-dim_names
  }
  amp_shift_mat_outsamp<-cbind(curvature_outsamp,lag_max_outsamp)
  dimnames(amp_shift_mat_outsamp)[[1]]<-dim_names
  if (mba)
  {
    dimnames(amp_shift_mat_outsamp)[[2]]<-c("Curvature","Peak-Correlation")
  } else
  {
    dimnames(amp_shift_mat_outsamp)[[2]]<-c("Curvature (out-of-sample)",
    "Peak-Correlation (out-of-sample)")
    dimnames(amp_shift_mat_outsamp)[[1]]<-dim_names
  }

  ats_mat<-cbind(Accuracy,Timeliness,Smoothness,Residual,
  Accuracy+Timeliness+Smoothness+Residual)
  dimnames(ats_mat)[[2]][dim(ats_mat)[2]]<-"Total MSE"
  dimnames(ats_mat)[[1]]<-dim_names
  return(list(xf=xf,amp_shift_mat_insamp=amp_shift_mat_insamp,
  amp_shift_mat_outsamp=amp_shift_mat_outsamp,
  ats_mat=ats_mat,amp=amp,shift=shift,b=b,xf_sym=xf_sym))
}


###################################################
### code chunk number 145: eco2.rnw:6019-6053
###################################################
# Generate series
gen_ser<-function(setseed,len,len1,cutoff_period)
{
  set.seed(setseed)

  # The longer sample is used for implementing the symmetric filter and
  # for computing peak-correlations
  x1<-rnorm(len1)
  # The shorter sample of length 120 is used for estimating the real-time filter
  x<-x1[1:len]
  plot_T<-F
  yhat<-x
  Gamma<-c(1,(1:(len/2))<len/(2*cutoff_period))
  # Compute periodogram (in-sample)
  weight_func<-per(x,plot_T)$per
  K<-length(weight_func)-1
  omega_k<-(0:K)*pi/K
  return(list(omega_k=omega_k,x1=x1,x=x,Gamma=Gamma,weight_func=weight_func))
}
#
setseed<-10
len<-120
len1<-1000
cutoff_period<-12
cutoff<-pi/cutoff_period
#
gen_obj<-gen_ser(setseed,len,len1,cutoff_period)
#
x1<-gen_obj$x1
x<-gen_obj$x
omega_k<-gen_obj$omega_k
Gamma<-gen_obj$Gamma
weight_func<-gen_obj$weight_func
K<-length(weight_func)-1


###################################################
### code chunk number 146: eco2.rnw:6069-6083
###################################################
# Generate series
L_sym<-61
L<-24
# Specify filter settings
Lag<-0
# No filter restrictions
i1<-i2<-F
# Optimize unconstrained real-time MSE-filter
eta_vec<-0
lambda_vec<-0
# True DGP: white noise spectrum
a1<-0
# mba<-T is used when replicating model-based approaches later on
mba<-F


###################################################
### code chunk number 147: eco2.rnw:6088-6099
###################################################
# Start processing
#
perf<-Performance_func(lambda_vec,eta_vec,len1,len,x1,weight_func,
Gamma,cutoff,L,L_sym,a1,mba)

xf<-perf$xf
amp_shift_mat_insamp<-perf$amp_shift_mat_insamp
amp_shift_mat_outsamp<-perf$amp_shift_mat_outsamp
ats_mat<-perf$ats_mat
amp<-perf$amp
shift<-perf$shift


###################################################
### code chunk number 148: eco2.rnw:6103-6145
###################################################
colo<-rainbow(dim(xf)[2])
file = paste("z_dfa_ar1_output_e.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,2))
ts.plot(x,main="DFA (red) and theoretical MSE (cyan) vs. input (black)",col="black")
mtext("Input", side = 3, line = -1,at=len/2,col="black",ylab="series")
for (i in 1:dim(xf)[2])
{
  lines(xf[,i],col=colo[i])
  mtext(dimnames(xf)[[2]][i], side = 3, line = -i-1,at=len/2,col=colo[i])
}
plot(per(x[L:len],plot_T)$per,type="l",axes=F,col="black",ylim=c(0,max(amp)),
ylab="Periodogram",xlab="",main="Periodograms")
mtext("Input", side = 3, line = -1,at=(K-L/2)/2,col="black")
lines(per(xf[L:len,2],plot_T)$per,lty=1,col="red")
mtext("MSE-Output", side = 3, line = -2,at=(K-L/2)/2,col="red")
axis(1,at=1+0:6*(K-L/2)/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,max(amp,1)),ylab="",xlab="",
main=paste("Amplitude MSE-Filter: L=",L,sep=""))
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:dim(xf)[2])
{
  lines(amp[,i],col=colo[i])
  mtext(dimnames(amp)[[2]][i], side = 3, line = -i-1,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(rep(0,K+1),type="l",axes=F,col="black",ylim=c(0,max(na.exclude(shift[2:(K+1),]))),
ylab="",xlab="",main=paste("Shift MSE-Filter: L=",L,sep=""))
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:dim(xf)[2])
{
  lines(shift[,i],col=colo[i])
  mtext(dimnames(shift)[[2]][i], side = 3, line = -i-1,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 149: z_dfa_ar1_output_e.pdf
###################################################
  file = paste("z_dfa_ar1_output_e", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Inputs (blue) and real-time lowpass MSE output (red) (top left); periodograms (top right);
  amplitude functions (bottom left) and time-shifts (bottom right): the theoretically best filter (blue)
  assumes knowledge of the DGP (white noise)", sep = "")
  cat("\\label{z_dfa_ar1_output_e}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 150: eco2.rnw:6185-6210
###################################################

# i1 restriction only (not i2)
i1<-T
i2<-F
# Compute filter and performances
perf<-Performance_func(lambda_vec,eta_vec,len1,len,x1,weight_func,Gamma,cutoff,L,L_sym,a1,mba)
xf_i1<-perf$xf
amp_i1<-perf$amp
shift_i1<-perf$shift
#
# i2 restriction only (not (i1)
i1<-F
i2<-T
perf<-Performance_func(lambda_vec,eta_vec,len1,len,x1,weight_func,Gamma,cutoff,L,L_sym,a1,mba)
xf_i2<-perf$xf
amp_i2<-perf$amp
shift_i2<-perf$shift
#
# i1 and i2 restrictions
i1<-T
i2<-T
perf<-Performance_func(lambda_vec,eta_vec,len1,len,x1,weight_func,Gamma,cutoff,L,L_sym,a1,mba)
xf_i1i2<-perf$xf
amp_i1i2<-perf$amp
shift_i1i2<-perf$shift


###################################################
### code chunk number 151: eco2.rnw:6217-6251
###################################################
#
# Generate graphs
file = paste("z_dfa_ar1_output_ee.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(1,2))
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,max(amp_i1i2,1)),
ylab="",xlab="",main="Amplitude")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
lines(amp[,2],lty=1,col="red")
lines(amp_i1[,2],lty=1,col="orange")
lines(amp_i2[,2],lty=1,col="green")
lines(amp_i1i2[,2],lty=1,col="violet")
mtext("Amplitude i1=i2=F", side = 3, line = -2,at=K/2,col="red")
mtext("Amplitude i1=T,i2=F", side = 3, line = -3,at=K/2,col="orange")
mtext("Amplitude i1=F,i2=T", side = 3, line = -4,at=K/2,col="green")
mtext("Amplitude i1=i2=T", side = 3, line = -5,at=K/2,col="violet")
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(rep(0,K+1),type="l",axes=F,col="black",ylim=c(0,max(na.exclude(shift_i1))),
ylab="",xlab="",main="Shift")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
lines(shift[,2],lty=1,col="red")
lines(shift_i1[,2],lty=1,col="orange")
lines(shift_i2[,2],lty=1,col="green")
lines(shift_i1i2[,2],lty=1,col="violet")
mtext("Shift i1=i2=F", side = 3, line = -2,at=K/2,col="red")
mtext("Shift i1=T,i2=F", side = 3, line = -3,at=K/2,col="orange")
mtext("Shift i1=F,i2=T", side = 3, line = -4,at=K/2,col="green")
mtext("Shift i1=i2=T", side = 3, line = -5,at=K/2,col="violet")
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 152: z_dfa_ar1_output_ee.pdf
###################################################
  file = paste("z_dfa_ar1_output_ee", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Effect of Filter restrictions in frequency zero on
  amplitude and time-shifts. Input: white noise.", sep = "")
  cat("\\label{z_dfa_ar1_output_ee}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 153: eco2.rnw:6304-6319
###################################################

eta_vec<-(0:6)/4
lambda_vec<-rep(0,length(eta_vec))
i1<-F
i2<-F
# Compute in/out-of-sample performances

perf<-Performance_func(lambda_vec,eta_vec,len1,len,x1,weight_func,Gamma,cutoff,L,L_sym,a1,mba)

xf<-perf$xf
amp_shift_mat_insamp<-perf$amp_shift_mat_insamp
amp_shift_mat_outsamp<-perf$amp_shift_mat_outsamp
ats_mat<-perf$ats_mat
amp<-perf$amp
shift<-perf$shift


###################################################
### code chunk number 154: eco2.rnw:6323-6351
###################################################
#
# Plots
colo<-rainbow(2*dim(amp)[2])
file = paste("z_dfa_cust_amp.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(1,2))
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,1),ylab="",xlab="",main="Amplitude")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:dim(amp)[2])
{
  lines(amp[,i],lty=1,col=colo[i])
  mtext(dimnames(amp)[[2]][i], side = 3, line = -1-i,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(rep(0,K+1),type="l",axes=F,col="black",ylim=c(0,max(na.exclude(shift))),
ylab="",xlab="",main="Shift")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:dim(amp)[2])
{
  lines(shift[,i],lty=1,col=colo[i])
  mtext(dimnames(shift)[[2]][i], side = 3, line = -1-i,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 155: z_dfa_cust_amp.pdf
###################################################
  file = paste("z_dfa_cust_amp", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Customized filters: the effect of eta (strength of noise supression/smoothness) on amplitude and time-shifts (delays)", sep = "")
  cat("\\label{z_dfa_cust_amp}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 156: ats_mat_1
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(ats_mat, dec = 1,digits=rep(3,dim(ats_mat)[2]+1),
  paste("ATS-components: emphasizing Smoothness only",sep=""),
  label=paste("ats_mat_1",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 157: amp_shift_mat_insamp_1
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(amp_shift_mat_insamp, dec = 1,digits=rep(3,dim(amp_shift_mat_insamp)[2]+1),
  paste("In-sample performance measures: Selectivity vs. Mean-Shift and Curvature vs. Peak-Correlation: emphasizing Smoothness only",sep=""),
  label=paste("amp_shift_mat_insamp_1",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 158: amp_shift_mat_outsamp_1
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(amp_shift_mat_outsamp, dec = 1,digits=rep(3,dim(amp_shift_mat_outsamp)[2]+1),
  paste("Out-of-sample performance measures: Selectivity vs. Mean-Shift and Curvature vs. Peak-Correlation: emphasizing Smoothness only",sep=""),
  label=paste("amp_shift_mat_outsamp_1",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 159: eco2.rnw:6422-6438
###################################################

#-------------------
# Compare MSE and heavily customized filter
xf0<-xf[L+1:len,1]
xf1<-xf[L+1:len,dim(amp)[2]]
file = paste("z_dfa_cust_amp_out_iut.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
ts.plot(as.ts(xf0),type="l",axes=F,col=colo[1],ylim=c(min(na.exclude(xf0)),
max(na.exclude(xf0))),
ylab="",xlab="",main=paste("Filter outputs: MSE (eta=0) vs. customized (eta=",
eta_vec[length(eta_vec)],")",sep=""))
mtext("MSE (eta=0)", side = 3, line = -1,at=len/2,col=colo[1])
lines(as.ts(xf1),col=colo[dim(amp)[2]])
mtext(paste("Customized (eta=",eta_vec[length(eta_vec)],")",sep=""),
side = 3, line = -2,at=len/2,col=colo[dim(amp)[2]])
dev.off()


###################################################
### code chunk number 160: z_dfa_cust_amp_out_iut.pdf
###################################################
  file = paste("z_dfa_cust_amp_out_iut", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Filter outputs: MSE vs. customized", sep = "")
  cat("\\label{z_dfa_cust_amp_out_iut}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 161: eco2.rnw:6467-6482
###################################################

lambda_vec<-c(0,2^(0:4))
eta_vec<-rep(0,length(lambda_vec))
i1<-F
i2<-F
# Compute in/out-of-sample performances

perf<-Performance_func(lambda_vec,eta_vec,len1,len,x1,weight_func,Gamma,cutoff,L,L_sym,a1,mba)

xf<-perf$xf
amp_shift_mat_insamp<-perf$amp_shift_mat_insamp
amp_shift_mat_outsamp<-perf$amp_shift_mat_outsamp
ats_mat<-perf$ats_mat
amp<-perf$amp
shift<-perf$shift


###################################################
### code chunk number 162: eco2.rnw:6487-6517
###################################################
#
# Plots
colo<-rainbow(2*dim(amp)[2])
file = paste("z_dfa_cust_shift.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(1,2))
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,1),ylab="",xlab="",main="Amplitude")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:dim(amp)[2])
{
  lines(amp[,i],lty=1,col=colo[i])
  mtext(dimnames(amp)[[2]][i], side = 3, line = -1-i,
  at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(rep(0,K+1),type="l",axes=F,col="black",ylim=c(min(shift),max(na.exclude(shift))),
ylab="",xlab="",main="Shift")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:dim(amp)[2])
{
  lines(shift[,i],lty=1,col=colo[i])
  mtext(dimnames(shift)[[2]][i], side = 3,
  line = -1-i,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 163: z_dfa_cust_shift.pdf
###################################################
  file = paste("z_dfa_cust_shift", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Customized filters: emphasize timeliness", sep = "")
  cat("\\label{z_dfa_cust_shift}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 164: ats_mat_2
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(ats_mat, dec = 1,digits=rep(3,dim(ats_mat)[2]+1),
  paste("ATS-components: emphasizing Timeliness only",sep=""),
  label=paste("ats_mat_2",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 165: amp_shift_mat_insamp_2
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(amp_shift_mat_insamp, dec = 1,digits=rep(3,dim(amp_shift_mat_insamp)[2]+1),
  paste("In-sample performances measures: Selectivity vs. Mean-Shift and Curvature vs. Peak-Correlation: emphasizing Timeliness only",sep=""),
  label=paste("amp_shift_mat_insamp_2",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 166: amp_shift_mat_outsamp_2
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(amp_shift_mat_outsamp, dec = 1,digits=rep(3,dim(amp_shift_mat_outsamp)[2]+1),
  paste("Out-of-sample performances measures: Selectivity vs. Mean-Shift and Curvature vs. Peak-Correlation: emphasizing Timeliness only",sep=""),
  label=paste("amp_shift_mat_outsamp_2",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 167: eco2.rnw:6589-6602
###################################################

#-------------------
# Compare MSE and heavily customized filter
xf0<-xf[L+1:len,1]
xf1<-xf[L+1:len,dim(amp)[2]]
file = paste("z_dfa_cust_amp_out_puti.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
ts.plot(as.ts(xf0),type="l",axes=F,col=colo[1],ylim=c(min(na.exclude(xf1)),max(na.exclude(xf1))),ylab="",xlab="",
main="Filter outputs: MSE (lambda=0) vs. customized (lambda=16)")
mtext("MSE (lambda=0)", side = 3, line = -1,at=len/2,col=colo[1])
lines(as.ts(xf1),col=colo[dim(amp)[2]])
mtext("Customized (lambda=16)", side = 3, line = -2,at=len/2,col=colo[dim(amp)[2]])
dev.off()


###################################################
### code chunk number 168: z_dfa_cust_amp_out_puti.pdf
###################################################
  file = paste("z_dfa_cust_amp_out_puti", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=4in]{", file, "}\n",sep = "")
  cat("\\caption{Filter outputs: MSE vs. customized", sep = "")
  cat("\\label{z_dfa_cust_amp_out_puti}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 169: eco2.rnw:6641-6661
###################################################
# Filter length
L<-24
# No filter restrictions
i1<-F
i2<-F
# Designs: MSE and two customized filters
lambda_vec<-c(0,32,100,8)
eta_vec<-c(0,1,0.7,1.5)
mba<-F
# Compute in/out-of-sample performances

perf<-Performance_func(lambda_vec,eta_vec,len1,len,x1,weight_func,Gamma,cutoff,L,L_sym,a1,mba)

xf<-perf$xf
amp_shift_mat_insamp<-perf$amp_shift_mat_insamp
amp_shift_mat_outsamp<-perf$amp_shift_mat_outsamp
ats_mat<-perf$ats_mat
amp<-perf$amp
shift<-perf$shift
xf_sym<-perf$xf_sym


###################################################
### code chunk number 170: eco2.rnw:6664-6694
###################################################
#
# Plots
colo<-rainbow(dim(amp)[2])
file = paste("z_dfa_cust_ats.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(1,2))
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,1),ylab="",xlab="",main="Amplitude")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:dim(amp)[2])
{
  lines(amp[,i],lty=1,col=colo[i])
  mtext(dimnames(amp)[[2]][i], side = 3,
  line = -1-i,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(rep(0,K+1),type="l",axes=F,col="black",ylim=c(min(shift),max(na.exclude(shift))),
ylab="",xlab="",main="Shift")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:dim(amp)[2])
{
  lines(shift[,i],lty=1,col=colo[i])
  mtext(dimnames(shift)[[2]][i], side = 3,
  line = -1-i,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 171: z_dfa_cust_ats.pdf
###################################################
  file = paste("z_dfa_cust_ats", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{MSE vs. customized filters: emphasize timeliness and smoothness simultaneously", sep = "")
  cat("\\label{z_dfa_cust_ats}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 172: ats_mat_2
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(ats_mat, dec = 1,digits=rep(3,dim(ats_mat)[2]+1),
  paste("ATS-components: emphasizing Timeliness only",sep=""),
  label=paste("ats_mat_3",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 173: amp_shift_mat_insamp_3
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(amp_shift_mat_insamp, dec = 1,digits=rep(3,dim(amp_shift_mat_insamp)[2]+1),
  paste("In-sample performances measures: Selectivity vs. Mean-Shift and Curvature vs. Peak-Correlation: emphasizing
  Smoothness and Timeliness",sep=""),
  label=paste("amp_shift_mat_insamp_3",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 174: amp_shift_mat_outsamp_3
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(amp_shift_mat_outsamp, dec = 1,digits=rep(3,dim(amp_shift_mat_outsamp)[2]+1),
  paste("Out-of-sample performances measures: Selectivity vs. Mean-Shift and Curvature vs. Peak-Correlation:
  emphasizing Smoothness and Timeliness",sep=""),
  label=paste("amp_shift_mat_outsamp_3",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 175: eco2.rnw:6769-6791
###################################################
colo<-rainbow(dim(amp)[2]+1)
file = paste("z_dfa_cust_amp_out_ats.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
# Select out-of-sample period
anf<-4*len#len+1
enf<-5*len
sel<-1:dim(amp)[2]
mplot<-scale(cbind(xf_sym,xf[,sel])[anf:enf,])
plot(as.ts(mplot[,1]),type="l",axes=F,col="black",ylim=c(min(na.exclude(mplot)),
max(na.exclude(mplot))),ylab="",xlab="",
main="Out-of-sample final symmetric vs. real-time MSE and customized",lwd=2)
mtext("Final symmetric", side = 3, line = -1,at=(enf-anf)/2,col="black")
for (i in 1:length(sel))
{
  lines(as.ts(mplot[,i+1]),col=colo[sel[i]],lwd=2)
  mtext(dimnames(amp)[[2]][i], side = 3, line = -1-i,at=(enf-anf)/2,col=colo[sel[i]])
}
axis(1,at=c(1,rep(0,6))+as.integer((0:6)*(enf-anf)/6),
labels=as.integer(anf+(0:6)*(enf-anf)/6))
axis(2)
box()
dev.off()


###################################################
### code chunk number 176: z_dfa_cust_amp_out_ats.pdf
###################################################
  file = paste("z_dfa_cust_amp_out_ats", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Standardized filter outputs out-of-sample: final symmetric vs. real-time MSE and customized designs", sep = "")
  cat("\\label{z_dfa_cust_amp_out_ats}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 177: eco2.rnw:6823-6874
###################################################

# Filter length
L<-24
# No filter restrictions
i1<-F
i2<-F
# Designs: MSE and two customized filters
lambda_vec<-c(0,32,100,8)
eta_vec<-c(0,1,0.7,1.5)
# Compute in/out-of-sample performances
setseed<-10
len<-120
len1<-1000
cutoff_period<-12
cutoff<-pi/cutoff_period
anzsim<-200
amp_shift_mat_insamp_array<-array(dim=c(dim(amp_shift_mat_insamp),anzsim))
amp_shift_mat_outsamp_array<-array(dim=c(dim(amp_shift_mat_outsamp),anzsim))
ats_mat_array<-array(dim=c(dim(ats_mat),anzsim))
amp_array<-array(dim=c(dim(amp),anzsim))
shift_array<-array(dim=c(dim(shift),anzsim))
dimnames(amp_shift_mat_insamp_array)<-list(dimnames(amp_shift_mat_insamp)[[1]],
dimnames(amp_shift_mat_insamp)[[2]],1:anzsim)
dimnames(amp_shift_mat_outsamp_array)<-list(dimnames(amp_shift_mat_outsamp)[[1]],
dimnames(amp_shift_mat_outsamp)[[2]],1:anzsim)
dimnames(ats_mat_array)<-list(dimnames(ats_mat)[[1]],
dimnames(ats_mat)[[2]],1:anzsim)
for (i in 1:anzsim)      #i<-1
{
  setseed<-setseed+1
# generate series and in-sample periodogram
  gen_obj<-gen_ser(setseed,len,len1,cutoff_period)
#
  x1<-gen_obj$x1
  x<-gen_obj$x
  omega_k<-gen_obj$omega_k
  Gamma<-gen_obj$Gamma
  weight_func<-gen_obj$weight_func
#
# Compute real-time filters
  perf<-Performance_func(lambda_vec,eta_vec,len1,len,x1,weight_func,
  Gamma,cutoff,L,L_sym,a1,mba)

  xf<-perf$xf
  amp_shift_mat_insamp_array[,,i]<-perf$amp_shift_mat_insamp
  amp_shift_mat_outsamp_array[,,i]<-perf$amp_shift_mat_outsamp
  ats_mat_array[,,i]<-perf$ats_mat
  amp_array[,,i]<-perf$amp
  shift_array[,,i]<-perf$shift
}
  


###################################################
### code chunk number 178: eco2.rnw:6877-6892
###################################################
colo<-rainbow(dim(amp)[2]+1)[1:dim(amp)[2]]
file = paste("z_box_plot_in.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
boxplot(list(amp_shift_mat_insamp_array[1,3,],amp_shift_mat_insamp_array[2,3,],
amp_shift_mat_insamp_array[3,3,],amp_shift_mat_insamp_array[4,3,],
amp_shift_mat_insamp_array[5,3,]),outline=F,
names=c("Best MSE",paste("DFA(",lambda_vec,",",eta_vec,")",sep="")),
main="Curvature in-sample",cex.axis=0.8,col=colo)
boxplot(list(amp_shift_mat_insamp_array[1,4,],amp_shift_mat_insamp_array[2,4,],
amp_shift_mat_insamp_array[3,4,],amp_shift_mat_insamp_array[4,4,],
amp_shift_mat_insamp_array[5,4,]),outline=F,
names=c("Best MSE",paste("DFA(",lambda_vec,",",eta_vec,")",sep="")),
main="Peak-correlation in-sample",cex.axis=0.8,col=colo)
dev.off()


###################################################
### code chunk number 179: z_box_plot_in.pdf
###################################################
  file = paste("z_box_plot_in", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Box-plots of in-of-sample empirical distributions of Curvature (top) and Peak-Correlation (bottom)", sep = "")
  cat("\\label{z_box_plot_in}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 180: eco2.rnw:6905-6919
###################################################
file = paste("z_box_plot_out.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
boxplot(list(amp_shift_mat_outsamp_array[1,1,],amp_shift_mat_outsamp_array[2,1,],
amp_shift_mat_outsamp_array[3,1,],amp_shift_mat_outsamp_array[4,1,],
amp_shift_mat_outsamp_array[5,1,]),outline=F,names=c("Best MSE",paste("DFA(",
lambda_vec,",",eta_vec,")",sep="")),#dimnames(amp_shift_mat_outsamp_array)[[1]],
main="Curvature out-of-sample",cex.axis=0.8,col=colo)
boxplot(list(amp_shift_mat_outsamp_array[1,2,],amp_shift_mat_outsamp_array[2,2,],
amp_shift_mat_outsamp_array[3,2,],amp_shift_mat_outsamp_array[4,2,],
amp_shift_mat_outsamp_array[5,2,]),outline=T,names=c("Best MSE",paste("DFA(",
lambda_vec,",",eta_vec,")",sep="")),#dimnames(amp_shift_mat_outsamp_array)[[1]],
main="Peak-correlation out-of-sample",cex.axis=0.8,col=colo)
dev.off()


###################################################
### code chunk number 181: z_box_plot_out.pdf
###################################################
  file = paste("z_box_plot_out", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Box-plots of out-of-sample empirical distributions of Curvature (top) and Peak-Correlation (bottom)", sep = "")
  cat("\\label{z_box_plot_out}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 182: eco2.rnw:6938-6959
###################################################
file = paste("z_box_plot_agg.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
boxplot(list(amp_shift_mat_outsamp_array[1,1,]*amp_shift_mat_outsamp_array[1,2,],
amp_shift_mat_outsamp_array[2,1,]*amp_shift_mat_outsamp_array[2,2,],
amp_shift_mat_outsamp_array[3,1,]*amp_shift_mat_outsamp_array[3,2,],
amp_shift_mat_outsamp_array[4,1,]*amp_shift_mat_outsamp_array[4,2,],
amp_shift_mat_outsamp_array[5,1,]*amp_shift_mat_outsamp_array[5,2,]),
outline=F,names=c("Best MSE",paste("DFA(",
lambda_vec,",",eta_vec,")",sep="")),#dimnames(amp_shift_mat_outsamp_array)[[1]],
main="Aggregate performance: Curvature*Peakcorr out-of-sample",cex.axis=0.8,col=colo)
boxplot(list(10*amp_shift_mat_outsamp_array[1,1,]+amp_shift_mat_outsamp_array[1,2,],
10*amp_shift_mat_outsamp_array[2,1,]+amp_shift_mat_outsamp_array[2,2,],
10*amp_shift_mat_outsamp_array[3,1,]+amp_shift_mat_outsamp_array[3,2,],
10*amp_shift_mat_outsamp_array[4,1,]+amp_shift_mat_outsamp_array[4,2,],
10*amp_shift_mat_outsamp_array[5,1,]+amp_shift_mat_outsamp_array[5,2,]),
outline=F,names=c("Best MSE",paste("DFA(",
lambda_vec,",",eta_vec,")",sep="")),#dimnames(amp_shift_mat_outsamp_array)[[1]],
main="Aggregate performance: 10*Curvature+Peakcorr out-of-sample",
cex.axis=0.8,col=colo)
dev.off()


###################################################
### code chunk number 183: z_box_plot_agg.pdf
###################################################
  file = paste("z_box_plot_agg", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Box-plots of product (top) and sum (bottom) of Curvature and Peak-Correlation: out-of-sample", sep = "")
  cat("\\label{z_box_plot_agg}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 184: eco2.rnw:7117-7130
###################################################
# Frequency resolution
K<-1200
# Hypothetical sample length
len<-120
# Symmetric lowpass target
cutoff<-pi/12
# Order of approximation : 10, 100, 1000, 10000
ord<-len
# Compute coefficients gamma
gamma_k<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))
#sum(gamma_k)+sum(gamma_k[2:ord])
# Target in frequency-domain
Gamma<-(0:K)<as.integer(cutoff*K/pi)+1


###################################################
### code chunk number 185: eco2.rnw:7133-7147
###################################################
omega_k<-(0:K)*pi/K
# AR(1)-coefficient
a1<-0.5
# Model-based filter
gamma_0<-gamma_k%*%a1^(0:ord)
gamma_mba_rt<-c(gamma_0,gamma_k[2:(ord)])
trffkt_mba<-rep(NA,K+1)
weight_func<-trffkt_mba
for (i in 0:K)
{
  trffkt_mba[i+1]<-gamma_mba_rt%*%exp(1.i*omega_k[i+1]*(0:(ord-1)))

}
amp_mba<-abs(trffkt_mba)


###################################################
### code chunk number 186: eco2.rnw:7152-7175
###################################################
L<-len
# model-based spectrum
weight_func<-1/abs(1-a1*exp(1.i*omega_k))^2
# MSE-filter
lambda<-0
eta<-0
# Real-time design
Lag<-0
# Unconstrained filter
i1<-F
i2<-F

dfa_ar1_long<-dfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2)

# Much shorter sample length
L<-len/10
lambda<-0
eta<-0
Lag<-0
i1<-F
i2<-F

dfa_ar1_short<-dfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2)


###################################################
### code chunk number 187: eco2.rnw:7178-7201
###################################################
file = paste("z_mbaedfa.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(2,1))
plot(dfa_ar1_long$b,type="l",axes=F,col="blue",ylim=c(min(dfa_ar1_long$b),
max(dfa_ar1_long$b)),ylab="",xlab="",main=paste("Filter coefficients of real-time
model-based (red) and DFA-MSE (blue): a1=",a1,", T=120",sep=""),lwd=2)
mtext("DFA", side = 3, line = -1,at=len/2,col="blue")
lines(gamma_mba_rt,col="red",lwd=1)
mtext("Model-based", side = 3, line = -2,at=len/2,col="red")
axis(1,at=c(1,rep(0,6))+as.integer((0:6)*len/6),
labels=as.integer(0+(0:6)*(len)/6))
axis(2)
box()
plot(dfa_ar1_short$b,type="l",axes=F,col="blue",ylim=c(min(dfa_ar1_short$b),
max(dfa_ar1_short$b)),ylab="",xlab="",main=paste("Filter coefficients of real-time
model-based (red) and DFA-MSE (blue): a1=",a1,", T=12",sep=""),lwd=2)
mtext("DFA", side = 3, line = -1,at=len/20,col="blue")
lines(gamma_mba_rt[1:12],col="red",lwd=1)
mtext("Model-based", side = 3, line = -2,at=len/20,col="red")
axis(1,at=1:12,labels=0:11)
axis(2)
box()
dev.off()


###################################################
### code chunk number 188: z_mbaedfa.pdf
###################################################
  file = paste("z_mbaedfa", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Filter coefficients of real-time model-based (red) and DFA-MSE (blue)", sep = "")
  cat("\\label{z_mbaedfa}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 189: eco2.rnw:7247-7258
###################################################
set.seed(10)
len<-120
len1<-1000
# The longer sample is used for implementing the symmetric filter and
# for computing peak-correlations
x1<-arima.sim(list(ar=a1),n=len1)
# In previous sections the shorter sample of length 120 was used for
# estimating the real-time DFA-filters
# We do not use in-sample data here because we know the DGP
#   (but we still have to define x in the head of the function)
x<-x1[1:len]


###################################################
### code chunk number 190: eco2.rnw:7263-7286
###################################################
# Length of symmetric target filter
L_sym<-61
# Length of model-based filters
L<-120
# Real-time
Lag<-0
# No filter restrictions
i1<-i2<-F
# Designs: three customized filters (the MSE is automatically computed)
lambda_vec<-c(32,100,500)
eta_vec<-c(0.5,1,0.3)
# Use model-based spectrum
mba<-T

# Compute in/out-of-sample performances
perf<-Performance_func(lambda_vec,eta_vec,len1,len,x1,weight_func,Gamma,cutoff,L,L_sym,a1,mba)

xf<-perf$xf
amp_shift_mat<-cbind(perf$amp_shift_mat_insamp,perf$amp_shift_mat_outsamp)
ats_mat<-perf$ats_mat
amp<-perf$amp
shift<-perf$shift
xf_sym<-perf$xf_sym


###################################################
### code chunk number 191: eco2.rnw:7290-7320
###################################################
#
# Plots
colo<-rainbow(dim(amp)[2]+1)
file = paste("z_dfa_cust_ats_mba.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(1,2))
plot(Gamma,type="l",axes=F,col="black",ylim=c(0,1),ylab="",xlab="",main="Amplitude")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:dim(amp)[2])
{
  lines(amp[,i],lty=1,col=colo[i])
  mtext(dimnames(amp)[[2]][i], side = 3,
  line = -1-i,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
plot(rep(0,K+1),type="l",axes=F,col="black",ylim=c(min(shift),max(na.exclude(shift))),
ylab="",xlab="",main="Shift")
mtext("Target", side = 3, line = -1,at=K/2,col="black")
for (i in 1:dim(amp)[2])
{
  lines(shift[,i],lty=1,col=colo[i])
  mtext(dimnames(shift)[[2]][i], side = 3,
  line = -1-i,at=K/2,col=colo[i])
}
axis(1,at=1+0:6*K/6,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
axis(2)
box()
dev.off()


###################################################
### code chunk number 192: z_dfa_cust_ats_mba.pdf
###################################################
  file = paste("z_dfa_cust_ats_mba", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=4in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{MBA vs. customized filters: amplitude and time-shifts", sep = "")
  cat("\\label{z_dfa_cust_ats_mba}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


###################################################
### code chunk number 193: amp_shift_mat_mba
###################################################
  library(Hmisc)
  require(xtable)
  #latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
  #center = "centering", file = "", floating = FALSE)
  xtable(amp_shift_mat, dec = 1,digits=rep(6,dim(amp_shift_mat)[2]+1),
  paste("In-sample performances measures: MBA-MSE vs. MBA-customized",sep=""),
  label=paste("amp_shift_mat_mba",sep=""),
  center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 194: eco2.rnw:7363-7384
###################################################
file = paste("z_dfa_cust_amp_out_ats_mba.pdf", sep = "")
pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
# Select out-of-sample period
anf<-3*len
enf<-5*len
sel<-1:dim(amp)[2]
mplot<-scale(cbind(xf_sym,xf)[anf:enf,])
plot(as.ts(mplot[,1]),type="l",axes=F,col="black",ylim=c(min(na.exclude(mplot)),
max(na.exclude(mplot))),ylab="",xlab="",
main="Final symmetric vs. real-time MBA: MSE and customized",lwd=2)
mtext("Final symmetric", side = 3, line = -1,at=(enf-anf)/2,col="black")
for (i in 1:length(sel))
{
  lines(as.ts(mplot[,i+1]),col=colo[sel[i]],lwd=2)
  mtext(dimnames(amp)[[2]][i], side = 3, line = -1-i,at=(enf-anf)/2,col=colo[sel[i]])
}
axis(1,at=c(1,rep(0,6))+as.integer((0:6)*(enf-anf)/6),
labels=as.integer(anf+(0:6)*(enf-anf)/6))
axis(2)
box()
dev.off()


###################################################
### code chunk number 195: z_dfa_cust_amp_out_ats_mba.pdf
###################################################
  file = paste("z_dfa_cust_amp_out_ats_mba", sep = "")
  cat("\\begin{figure}[H]")
  cat("\\begin{center}")
  cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
  cat("\\caption{Filter outputs out-of-sample: final symmetric vs. real-time MSE and customized designs", sep = "")
  cat("\\label{z_dfa_cust_amp_out_ats_mba}}", sep = "")
  cat("\\end{center}")
  cat("\\end{figure}")


