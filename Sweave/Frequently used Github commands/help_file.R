install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
install.packages("devtools")

devtools::install_github("hadley/devtools")


devtools::document()



# Generate help files (folder man)
devtools::document()
# run checks
devtools::check()


# In the environment panel check that upper right entry is on "master": otherwise pull/push are not possible
#

# Change file, save it, commit (stage all relevant files), push

# To export a file (add to package on github)
#   Supply the following line in the description of the functions
#' @export



