disk_id<-"C"

# Set paths
path.main<-paste(disk_id,":\\wia_desktop\\Projekte\\2017\\MDFA-Legacy\\Sweave\\",sep="")
path.pgm<-paste(path.main,"Rnw\\",sep="")
path.out<-paste(path.main,"Latex\\",sep="")


script <- paste(path.pgm,"MDFA_Legacy",sep="")

## enforce par(ask=FALSE)
options(device.ask.default=FALSE)

## create a LaTeX file
Sweave(script,output=paste(path.out,"MDFA_Legacy.tex",sep=""))


# The following piece of code is useful when extracting the R-code


script_e <- paste(path.main,"Rnw\\MDFA_Legacy.Rnw",sep="")

## create an R source file from the code chunks
# The newly generated file MDFA_Legay.r will be written to wia_desktop\\Projekte\\2014\\MDFA-Legacy\\

Stangle(script_e)

