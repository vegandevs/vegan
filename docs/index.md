# vegan: an R package for community ecologists

Ordination methods, diversity analysis and other functions for community
and vegetation ecologists.

Website for the development version of the **vegan** package.

Vignettes are available on
[R-universe](https://vegandevs.r-universe.dev/vegan)

- [Introduction to ordination in
  vegan](https://vegandevs.r-universe.dev/vegan/doc/intro-vegan.pdf)
- [Partition of
  Variation](https://vegandevs.r-universe.dev/vegan/doc/partitioning.pdf)
- [Diversity analysis in
  vegan](https://vegandevs.r-universe.dev/vegan/doc/diversity-vegan.pdf)
- [Design decisions and
  implementation](https://vegandevs.r-universe.dev/vegan/doc/decision-vegan.pdf)
- [vegan FAQ](https://vegandevs.github.io/vegan/articles/FAQ-vegan.html)

# Installation

To install the development version of **vegan** you can use the usual
`git` and `R CMD build -> R CMD INSTALL` dance on the cloned repo (or
downloaded sources). You’ll need to be able to install packages from
source for that to work; if you don’t have the relevant developer tools,
you won’t be able to install **vegan** this way.

## Using **remotes**

If you do have the developer tools installed but don’t want the hassle
of keeping a local source code tree up-to-date, use the **remotes**
package:

``` r
install.packages("remotes")
remotes::install_github("vegandevs/vegan")
```

## Installing binaries from R Universe

If you just want to install a binary version of the packages, just as
you would from CRAN, you can install from our R Universe repository. Run
the following in your R session:

``` r
install.packages('vegan',
    repos = c('https://vegandevs.r-universe.dev','https://cloud.r-project.org'))
```
