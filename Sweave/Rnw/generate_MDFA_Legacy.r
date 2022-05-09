

rm(list=ls())


library(xts)
# State-space models (will be replicated by MDFA) 
library(dlm)
# Numerical package 
library(numDeriv)
# Graphical package for recession-shading (empirical examples based on US-GDP)
library(tis) 
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
devtools::install_github("wiaidp/MDFA")
# MDFA package
library(MDFA) 
library(RCurl)    # For getURL() and curl handler / cookie / google login
library(stringr)  # For str_trim() to trip whitespace from strings
library(Quandl)
require (Quandl)
#Quandl.api_key("ivVdJGV57TXA1RX5jgvp")
#  Tucker's libraries
library(expm)
library(R.utils)
library(xtable)

# set directory to GitHub/MDFA-Legacy
setwd("C:/Users/neide/OneDrive/Documents/GitHub/MDFA-Legacy")
path.main <- paste(getwd(),"/Sweave/",sep="")
path.pgm <- paste(path.main,"Rnw/",sep="")
path.out <- paste(path.main,"Latex/",sep="")
 
## Load MDFA R codes #
sourceDirectory(path=paste(path.main,"RcodeTSM/",sep=""), modifiedOnly=FALSE)

script <- paste(path.pgm,"MDFA_Legacy",sep="")

## enforce par(ask=FALSE)
options(device.ask.default=FALSE)

## create a LaTeX file
Sweave(script,output=paste(path.out,"MDFA_Legacy.tex",sep=""))


# The following piece of code is useful when extracting the R-code


script_e <- paste(path.main,"Rnw/MDFA_Legacy.Rnw",sep="")

## create an R source file from the code chunks
# The newly generated file MDFA_Legay.r will be written to wia_desktop\\Projekte\\2014\\MDFA-Legacy\\

Stangle(script_e)

