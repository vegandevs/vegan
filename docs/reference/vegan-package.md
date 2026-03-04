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
data(varespec)
data(varechem)
ord <- metaMDS(varespec)
#> Square root transformation
#> Wisconsin double standardization
#> Run 0 stress 0.1843196 
#> Run 1 stress 0.2291372 
#> Run 2 stress 0.2265716 
#> Run 3 stress 0.18584 
#> Run 4 stress 0.196245 
#> Run 5 stress 0.241719 
#> Run 6 stress 0.2088293 
#> Run 7 stress 0.1852397 
#> Run 8 stress 0.2104573 
#> Run 9 stress 0.2058803 
#> Run 10 stress 0.1948413 
#> Run 11 stress 0.1982376 
#> Run 12 stress 0.2079059 
#> Run 13 stress 0.2251281 
#> Run 14 stress 0.2365003 
#> Run 15 stress 0.2382236 
#> Run 16 stress 0.195049 
#> Run 17 stress 0.2143611 
#> Run 18 stress 0.1948413 
#> Run 19 stress 0.2069728 
#> Run 20 stress 0.1955837 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     20: stress ratio > sratmax
plot(ord, type = "t")
## Fit environmental variables
ef <- envfit(ord, varechem)
ef
#> 
#> ***VECTORS
#> 
#>             NMDS1    NMDS2     r2 Pr(>r)    
#> N        -0.05039 -0.99873 0.2081  0.091 .  
#> P         0.68707  0.72659 0.1755  0.139    
#> K         0.82733  0.56171 0.1657  0.151    
#> Ca        0.75017  0.66124 0.2811  0.032 *  
#> Mg        0.69680  0.71727 0.3493  0.013 *  
#> S         0.27631  0.96107 0.1774  0.127    
#> Al       -0.83771  0.54611 0.5155  0.002 ** 
#> Fe       -0.86190  0.50708 0.4000  0.003 ** 
#> Mn        0.80229 -0.59693 0.5322  0.001 ***
#> Zn        0.66523  0.74664 0.1779  0.124    
#> Mo       -0.84873  0.52882 0.0517  0.590    
#> Baresoil  0.87204 -0.48943 0.2494  0.057 .  
#> Humdepth  0.92632 -0.37673 0.5590  0.001 ***
#> pH       -0.79905  0.60127 0.2624  0.044 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Permutation: free
#> Number of permutations: 999
#> 
#> 
plot(ef, p.max = 0.05)

### Example 2: Constrained ordination (RDA)
## The example uses formula interface to define the model
data(dune)
data(dune.env)
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
plot(mod0)

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
#> + Management  3 87.082 2.8400  0.005 **
#> + Moisture    3 87.707 2.5883  0.005 **
#> + Manure      4 89.232 1.9539  0.005 **
#> + A1          1 89.591 1.9217  0.055 . 
#> + Use         2 91.032 1.1741  0.310   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Step: dune ~ Management 
#> 
#>              Df   AIC    F Pr(>F)   
#> - Management  3 89.62 2.84  0.005 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>            Df    AIC      F Pr(>F)   
#> + Moisture  3 85.567 1.9764  0.005 **
#> + Manure    3 87.517 1.3902  0.155   
#> + A1        1 87.424 1.2965  0.205   
#> + Use       2 88.284 1.0510  0.430   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Step: dune ~ Management + Moisture 
#> 
#>              Df    AIC      F Pr(>F)   
#> - Moisture    3 87.082 1.9764   0.02 * 
#> - Management  3 87.707 2.1769   0.01 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>          Df    AIC      F Pr(>F)
#> + Manure  3 85.762 1.1225  0.290
#> + A1      1 86.220 0.8359  0.540
#> + Use     2 86.842 0.8027  0.665
#> 
mod
#> 
#> Call: rda(formula = dune ~ Management + Moisture, data = dune.env)
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
plot(mod)

## Permutation test for all variables
anova(mod)
#> Permutation test for rda under reduced model
#> Permutation: free
#> Number of permutations: 999
#> 
#> Model: rda(formula = dune ~ Management + Moisture, data = dune.env)
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
#> Model: rda(formula = dune ~ Management + Moisture, data = dune.env)
#>            Df Variance      F Pr(>F)   
#> Management  3   18.938 2.1769  0.006 **
#> Moisture    3   17.194 1.9764  0.003 **
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
adonis2(dune ~ ., dune.env)
#> Permutation test for adonis under reduced model
#> Permutation: free
#> Number of permutations: 999
#> 
#> adonis2(formula = dune ~ ., data = dune.env)
#>          Df SumOfSqs      R2      F Pr(>F)  
#> Model    12   3.3265 0.77379 1.9954  0.011 *
#> Residual  7   0.9725 0.22621                
#> Total    19   4.2990 1.00000                
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis2(dune ~ Management + Moisture, dune.env)
#> Permutation test for adonis under reduced model
#> Permutation: free
#> Number of permutations: 999
#> 
#> adonis2(formula = dune ~ Management + Moisture, data = dune.env)
#>          Df SumOfSqs      R2      F Pr(>F)    
#> Model     6   2.6202 0.60949 3.3816  0.001 ***
#> Residual 13   1.6788 0.39051                  
#> Total    19   4.2990 1.00000                  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
