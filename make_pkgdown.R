#'--------------------------------
#'
#' A script to build a package website
#'
#'--------------------------------

#set the working directory here
setwd(here::here())

#build the site
pkgdown::build_site()
