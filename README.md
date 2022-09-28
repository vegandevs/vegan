# vegan: an R package for community ecologists

Ordination methods, diversity analysis and other functions for community and
vegetation ecologists.

<!-- badges: start -->
[![R build
status](https://github.com/vegandevs/vegan/workflows/R-CMD-check/badge.svg)](https://github.com/vegandevs/vegan/actions)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/vegan)](https://cran.r-project.org/package=vegan)
[![](http://cranlogs.r-pkg.org/badges/grand-total/vegan)](http://cran.rstudio.com/web/packages/vegan/index.html)
<!-- badges: end -->

Website for the development version of the **vegan** package.

# Installation

To install the development version of **vegan** you can use the usual `git` and `R CMD build -> R CMD INSTALL` dance on the cloned repo (or downloaded sources). You'll need to be able to install packages from source for that to work; if you don't have the relevant developer tools, you won't be able to install **vegan** this way.

## Using **remotes**

If you do have the developer tools installed but don't want the hassle of keeping a local source code tree up-to-date, use the **remotes** package:

```r
install.packages("remotes")
remotes::install_github("vegandevs/vegan")
```

## Installing binaries from R Universe

If you just want to install a binary version of the packages, just as you would from CRAN, you can install from our R Universe repository. Run the following in your R session:

```r
# Enable repository from vegandevs
options(repos = c(
  vegandevs = 'https://vegandevs.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))
# Download and install vegan in R
install.packages('vegan')
```

To make this permanent, you'll need to include the `options()` part of that in your person R settings file, typically a `.Rprofile` profile in your home drive. See `?Startup`.
