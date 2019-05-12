

# Latest R-studio?
install.packages("rstudioapi")
rstudioapi::isAvailable("0.99.149")

install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
install.packages("devtools")

devtools::install_github("hadley/devtools")
library(devtools)
library(roxygen2)
library(testthat)
library(knitr)
has_devel()




devtools::build_github_devtools()

#### Restart R before continuing ####
install.packages("devtools.zip", repos = NULL)

# Remove the package after installation
unlink("devtools.zip")


#---------------------
# Make a git repository (will highlight the git symbol in the source panel of rstudio (with commit/pull/push))

# Go to the projects menu in the top-right corner of RStudio.
# Click and select project options
# Click on Git/Svn
# In the version control system checkbox select git (instead of none)

#----------------------------------


.libPaths()
.Library
.Library.site
.libPaths(new=.Library)
Sys.getenv("R_LIBS_USER")
Sys.getenv("R_LIBS_SITE")
Sys.getenv("R_HOME")

install.packages("pacman")
library(pacman)
p_unlock()

Sys.glob(path.expand(new))
Sys.glob(path.expand("~/foo"))
path.expand("~/foo")
paths <- unique(normalizePath(c(new, .Library.site, .Library), "/"))

ret<-1:6
write.table(ret,file="C:/PROGRA~1/R/R-33~1.3/library/ret.txt")

write.table(ret,file="C:/Program Files/R/R-3.3.3/library/ret.txt")

read.table("C:/Users/Marc/Documents/ret.txt")


# Generate help files (folder man)
devtools::document()
# run checks
devtools::check()


# In the environment panel check that upper right entry is on "master": otherwise pull/push are not possible
#

# Change file, save it, commit (stage all relevant files), push
