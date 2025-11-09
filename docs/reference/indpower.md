# Indicator Power of Species

Indicator power calculation of Halme et al. (2009) or the congruence
between indicator and target species.

## Usage

``` r
indpower(x, type = 0)
```

## Arguments

- x:

  Community data frame or matrix.

- type:

  The type of statistic to be returned. See Details for explanation.

## Details

Halme et al. (2009) described an index of indicator power defined as
\\IP_I = \sqrt{a \times b}\\, where \\a = S / O_I\\ and \\b = 1 - (O_T -
S) / (N - O_I)\\. \\N\\ is the number of sites, \\S\\ is the number of
shared occurrences of the indicator (\\I\\) and the target (\\T\\)
species. \\O_I\\ and \\O_T\\ are number of occurrences of the indicator
and target species. The `type` argument in the function call enables to
choose which statistic to return. `type = 0` returns \\IP_I\\,
`type = 1` returns \\a\\, `type = 2` returns \\b\\. Total indicator
power (TIP) of an indicator species is the column mean (without its own
value, see examples). Halme et al. (2009) explain how to calculate
confidence intervals for these statistics, see Examples.

## Value

A matrix with indicator species as rows and target species as columns
(this is indicated by the first letters of the row/column names).

## References

Halme, P., Mönkkönen, M., Kotiaho, J. S, Ylisirniö, A-L. 2009.
Quantifying the indicator power of an indicator species. *Conservation
Biology* 23: 1008–1016.

## Author

Peter Solymos

## Examples

``` r
data(dune)
## IP values
ip <- indpower(dune)
## and TIP values
diag(ip) <- NA
(TIP <- rowMeans(ip, na.rm=TRUE))
#> i.Achimill i.Agrostol i.Airaprae i.Alopgeni i.Anthodor i.Bellpere i.Bromhord 
#>  0.3186250  0.3342800  0.2168133  0.3416198  0.3567884  0.3432281  0.3665632 
#> i.Chenalbu i.Cirsarve i.Comapalu i.Eleopalu i.Elymrepe i.Empenigr i.Hyporadi 
#>  0.2095044  0.2781640  0.1713273  0.2414787  0.3263516  0.2016196  0.2378197 
#> i.Juncarti i.Juncbufo i.Lolipere i.Planlanc  i.Poaprat  i.Poatriv i.Ranuflam 
#>  0.2915850  0.3331330  0.3998442  0.3426064  0.4094319  0.3929520  0.2663080 
#> i.Rumeacet i.Sagiproc i.Salirepe i.Scorautu i.Trifprat i.Trifrepe i.Vicilath 
#>  0.3484684  0.3788905  0.2898512  0.4362493  0.3145854  0.4503764  0.2605349 
#> i.Bracruta i.Callcusp 
#>  0.4252676  0.2070766 

## p value calculation for a species
## from Halme et al. 2009
## i is ID for the species
i <- 1
fun <- function(x, i) indpower(x)[i,-i]
## 'c0' randomizes species occurrences
os <- oecosimu(dune, fun, "c0", i=i, nsimul=99)
#> Warning: nullmodel transformed 'comm' to binary data
## get z values from oecosimu output
z <- os$oecosimu$z
## p-value
(p <- sum(z) / sqrt(length(z)))
#> [1] -1.626625
## 'heterogeneity' measure
(chi2 <- sum((z - mean(z))^2))
#> [1] 92.62267
pchisq(chi2, df=length(z)-1)
#> [1] 1
## Halme et al.'s suggested output
out <- c(TIP=TIP[i], 
    significance=p,
    heterogeneity=chi2,
    minIP=min(fun(dune, i=i)),
    varIP=sd(fun(dune, i=i)^2))
out
#> TIP.i.Achimill   significance  heterogeneity          minIP          varIP 
#>      0.3186250     -1.6266254     92.6226720      0.0000000      0.2142097 
```
