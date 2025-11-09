# Community Ecology Package: Ordination, Diversity and Dissimilarities

The vegan package provides tools for descriptive community ecology. It
has most basic functions of diversity analysis, community ordination and
dissimilarity analysis. Most of its multivariate tools can be used for
other data types as well.

## Details

The functions in the vegan package contain tools for diversity analysis,
ordination methods and tools for the analysis of dissimilarities.
Together with the labdsv package, the vegan package provides most
standard tools of descriptive community analysis. Package ade4 provides
an alternative comprehensive package, and several other packages
complement vegan and provide tools for deeper analysis in specific
fields. Package BiodiversityR provides a GUI for a large subset of vegan
functionality.

The vegan package is developed at GitHub
(<https://github.com/vegandevs/vegan/>). GitHub provides up-to-date
information and forums for bug reports.

Most important changes in vegan documents can be read with
`news(package="vegan")` and vignettes can be browsed with
`browseVignettes("vegan")`. The vignettes include a vegan FAQ,
discussion on design decisions, short introduction to ordination and
discussion on diversity methods.

To see the preferable citation of the package, type `citation("vegan")`.

## Author

The vegan development team is Jari Oksanen, F. Guillaume Blanchet,
Roeland Kindt, Pierre Legendre, Peter R. Minchin, R. B. O'Hara, Gavin L.
Simpson, Peter Solymos, M. Henry H. Stevens, Helene Wagner. Many other
people have contributed to individual functions: see credits in function
help pages.

## Examples

``` r
### Example 1: Unconstrained ordination
## NMDS
data(varespec, varechem)
ord <- metaMDS(varespec)
#> Square root transformation
#> Wisconsin double standardization
#> Run 0 stress 0.1843196 
#> Run 1 stress 0.2067744 
#> Run 2 stress 0.1976151 
#> Run 3 stress 0.2444586 
#> Run 4 stress 0.2095882 
#> Run 5 stress 0.233192 
#> Run 6 stress 0.1825658 
#> ... New best solution
#> ... Procrustes: rmse 0.04162469  max resid 0.1517992 
#> Run 7 stress 0.2295042 
#> Run 8 stress 0.2092456 
#> Run 9 stress 0.1825658 
#> ... Procrustes: rmse 3.746352e-06  max resid 1.292456e-05 
#> ... Similar to previous best
#> Run 10 stress 0.2092456 
#> Run 11 stress 0.2290222 
#> Run 12 stress 0.237339 
#> Run 13 stress 0.1869637 
#> Run 14 stress 0.2807394 
#> Run 15 stress 0.1955836 
#> Run 16 stress 0.1969805 
#> Run 17 stress 0.2066027 
#> Run 18 stress 0.2088293 
#> Run 19 stress 0.1869637 
#> Run 20 stress 0.1948413 
#> *** Best solution repeated 1 times
plot(ord, optimize = TRUE, type = "t")
## Fit environmental variables
ef <- envfit(ord, varechem)
ef
#> 
#> ***VECTORS
#> 
#>             NMDS1    NMDS2     r2 Pr(>r)    
#> N        -0.05730 -0.99836 0.2536  0.042 *  
#> P         0.61972  0.78482 0.1938  0.111    
#> K         0.76645  0.64230 0.1809  0.139    
#> Ca        0.68519  0.72836 0.4119  0.004 ** 
#> Mg        0.63252  0.77455 0.4270  0.003 ** 
#> S         0.19138  0.98152 0.1752  0.124    
#> Al       -0.87160  0.49021 0.5269  0.001 ***
#> Fe       -0.93602  0.35194 0.4450  0.005 ** 
#> Mn        0.79871 -0.60172 0.5231  0.001 ***
#> Zn        0.61755  0.78653 0.1879  0.106    
#> Mo       -0.90308  0.42948 0.0609  0.511    
#> Baresoil  0.92489 -0.38024 0.2508  0.044 *  
#> Humdepth  0.93283 -0.36032 0.5200  0.001 ***
#> pH       -0.64798  0.76166 0.2308  0.069 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Permutation: free
#> Number of permutations: 999
#> 
#> 
plot(ef, p.max = 0.05)

### Example 2: Constrained ordination (RDA)
## The example uses formula interface to define the model
data(dune, dune.env)
## No constraints: PCA
mod0 <- rda(dune ~ 1, dune.env)
mod0
#> 
#> Call: rda(formula = dune ~ 1, data = dune.env)
#> 
#>               Inertia Rank
#> Total           84.12     
#> Unconstrained   84.12   19
#> 
#> Inertia is variance
#> 
#> Eigenvalues for unconstrained axes:
#>    PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
#> 24.795 18.147  7.629  7.153  5.695  4.333  3.199  2.782 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
plot(mod0, spe.par = list(arrows = TRUE))

## All environmental variables: Full model
mod1 <- rda(dune ~ ., dune.env)
#> 
#> Some constraints or conditions were aliased because they were redundant. This
#> can happen if terms are constant or linearly dependent (collinear): ‘Manure^4’
mod1
#> 
#> Call: rda(formula = dune ~ A1 + Moisture + Management + Use + Manure, data
#> = dune.env)
#> 
#>               Inertia Proportion Rank
#> Total         84.1237     1.0000     
#> Constrained   63.2062     0.7513   12
#> Unconstrained 20.9175     0.2487    7
#> 
#> Inertia is variance
#> 
#> -- NOTE:
#> Some constraints or conditions were aliased because they were redundant.
#> This can happen if terms are constant or linearly dependent (collinear):
#> ‘Manure^4’
#> 
#> Eigenvalues for constrained axes:
#>   RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7   RDA8   RDA9  RDA10  RDA11 
#> 22.396 16.208  7.039  4.038  3.760  2.609  2.167  1.803  1.404  0.917  0.582 
#>  RDA12 
#>  0.284 
#> 
#> Eigenvalues for unconstrained axes:
#>   PC1   PC2   PC3   PC4   PC5   PC6   PC7 
#> 6.627 4.309 3.549 2.546 2.340 0.934 0.612 
#> 
plot(mod1)

## Automatic selection of variables by permutation P-values
mod <- ordistep(mod0, scope=formula(mod1))
#> 
#> Start: dune ~ 1 
#> 
#>              Df    AIC      F Pr(>F)   
#> + Moisture    3 87.707 2.5883  0.005 **
#> + Management  3 87.082 2.8400  0.010 **
#> + Manure      4 89.232 1.9539  0.015 * 
#> + A1          1 89.591 1.9217  0.050 * 
#> + Use         2 91.032 1.1741  0.245   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Step: dune ~ Moisture 
#> 
#>            Df   AIC      F Pr(>F)   
#> - Moisture  3 89.62 2.5883  0.005 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>              Df    AIC      F Pr(>F)   
#> + Management  3 85.567 2.1769   0.01 **
#> + Manure      4 86.060 1.8598   0.01 **
#> + Use         2 88.003 1.4245   0.13   
#> + A1          1 88.886 0.6286   0.77   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Step: dune ~ Moisture + Management 
#> 
#>              Df    AIC      F Pr(>F)  
#> - Moisture    3 87.082 1.9764  0.020 *
#> - Management  3 87.707 2.1769  0.015 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>          Df    AIC      F Pr(>F)
#> + Manure  3 85.762 1.1225  0.340
#> + A1      1 86.220 0.8359  0.610
#> + Use     2 86.842 0.8027  0.645
#> 
mod
#> 
#> Call: rda(formula = dune ~ Moisture + Management, data = dune.env)
#> 
#>               Inertia Proportion Rank
#> Total         84.1237     1.0000     
#> Constrained   46.4249     0.5519    6
#> Unconstrained 37.6988     0.4481   13
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>   RDA1   RDA2   RDA3   RDA4   RDA5   RDA6 
#> 21.588 14.075  4.123  3.163  2.369  1.107 
#> 
#> Eigenvalues for unconstrained axes:
#>   PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8   PC9  PC10  PC11  PC12  PC13 
#> 8.241 7.138 5.355 4.409 3.143 2.770 1.878 1.741 0.952 0.909 0.627 0.311 0.227 
#> 
plot(mod, spe.par = list(optimize = TRUE))

## Permutation test for all variables
anova(mod)
#> Permutation test for rda under reduced model
#> Permutation: free
#> Number of permutations: 999
#> 
#> Model: rda(formula = dune ~ Moisture + Management, data = dune.env)
#>          Df Variance      F Pr(>F)    
#> Model     6   46.425 2.6682  0.001 ***
#> Residual 13   37.699                  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## Permutation test of "type III" effects, or significance when a term
## is added to the model after all other terms
anova(mod, by = "margin")
#> Permutation test for rda under reduced model
#> Marginal effects of terms
#> Permutation: free
#> Number of permutations: 999
#> 
#> Model: rda(formula = dune ~ Moisture + Management, data = dune.env)
#>            Df Variance      F Pr(>F)   
#> Moisture    3   17.194 1.9764  0.009 **
#> Management  3   18.938 2.1769  0.004 **
#> Residual   13   37.699                 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## Plot only sample plots, use different symbols and draw SD ellipses 
## for Managemenet classes
plot(mod, display = "sites", type = "n")
with(dune.env, points(mod, disp = "si", pch = as.numeric(Management)))
with(dune.env, legend("topleft", levels(Management), pch = 1:4,
  title = "Management"))
with(dune.env, ordiellipse(mod, Management, label = TRUE))
## add fitted surface of diversity to the model
ordisurf(mod, diversity(dune), add = TRUE)

#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> y ~ s(x1, x2, k = 10, bs = "tp", fx = FALSE)
#> 
#> Estimated degrees of freedom:
#> 1.28  total = 2.28 
#> 
#> REML score: 3.00623     
### Example 3: analysis of dissimilarites a.k.a. non-parametric
### permutational anova
adonis2(dune ~ ., dune.env, by = "margin")
#> Permutation test for adonis under reduced model
#> Marginal effects of terms
#> Permutation: free
#> Number of permutations: 999
#> 
#> adonis2(formula = dune ~ ., data = dune.env, by = "margin")
#>            Df SumOfSqs      R2      F Pr(>F)
#> A1          1   0.1283 0.02983 0.9231  0.466
#> Moisture    3   0.6596 0.15343 1.5826  0.115
#> Management  2   0.1959 0.04556 0.7050  0.724
#> Use         2   0.1305 0.03036 0.4697  0.900
#> Manure      3   0.4208 0.09787 1.0096  0.460
#> Residual    7   0.9725 0.22621              
#> Total      19   4.2990 1.00000              
adonis2(dune ~ Management + Moisture, dune.env, by = "term")
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Permutation: free
#> Number of permutations: 999
#> 
#> adonis2(formula = dune ~ Management + Moisture, data = dune.env, by = "term")
#>            Df SumOfSqs      R2      F Pr(>F)    
#> Management  3   1.4686 0.34161 3.7907  0.002 ** 
#> Moisture    3   1.1516 0.26788 2.9726  0.001 ***
#> Residual   13   1.6788 0.39051                  
#> Total      19   4.2990 1.00000                  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
