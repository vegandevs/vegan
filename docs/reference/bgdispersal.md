# Coefficients of Biogeographical Dispersal Direction

This function computes coefficients of dispersal direction between
geographically connected areas, as defined by Legendre and Legendre
(1984), and also described in Legendre and Legendre (2012, section
13.3.4).

## Usage

``` r
bgdispersal(mat, PAonly = FALSE, abc = FALSE)
```

## Arguments

- mat:

  Data frame or matrix containing a community composition data table
  (species presence-absence or abundance data).

- PAonly:

  `FALSE` if the four types of coefficients, DD1 to DD4, are requested;
  `TRUE` if `DD1` and `DD2` only are sought (see Details).

- abc:

  If `TRUE`, return tables `a`, `b` and `c` used in `DD1` and `DD2`.

## Details

The signs of the DD coefficients indicate the direction of dispersal,
provided that the asymmetry is significant. A positive sign indicates
dispersal from the first (row in DD tables) to the second region
(column); a negative sign indicates the opposite. A McNemar test of
asymmetry is computed from the presence-absence data to test the
hypothesis of a significant asymmetry between the two areas under
comparison.

In the input data table, the rows are sites or areas, the columns are
taxa. Most often, the taxa are species, but the coefficients can be
computed from genera or families as well. DD1 and DD2 only are computed
for presence-absence data. The four types of coefficients are computed
for quantitative data, which are converted to presence-absence for the
computation of DD1 and DD2. `PAonly = FALSE` indicates that the four
types of coefficients are requested. `PAonly = TRUE` if DD1 and DD2 only
are sought.

## Value

Function `bgdispersal` returns a list containing the following matrices:

- DD1 :

  \\DD1\_{j,k} = (a(b - c))/((a + b + c)^2)\\

- DD2 :

  \\DD2\_{j,k} = (2 a (b - c))/((2a + b + c) (a + b + c))\\ where \\a\\,
  \\b\\, and \\c\\ have the same meaning as in the computation of binary
  similarity coefficients.

- DD3 :

  \\DD3\_{j,k} = {W(A-B) / (A+B-W)^2} \\

- DD4 :

  \\DD4\_{j,k} = 2W(A-B) / ((A+B)(A+B-W))\\ where
  `W = sum(pmin(vector1, vector2))`, `A = sum(vector1)`,
  `B = sum(vector2)`

- McNemar :

  McNemar chi-square statistic of asymmetry (Sokal and Rohlf 1995):
  \\2(b \log(b) + c \log(c) - (b+c) \log((b+c)/2)) / q\\, where \\q =
  1 + 1/(2(b+c))\\ (Williams correction for continuity)

- prob.McNemar :

  probabilities associated with McNemar statistics, chi-square test. H0:
  no asymmetry in \\(b-c)\\.

## References

Legendre, P. and V. Legendre. 1984. Postglacial dispersal of freshwater
fishes in the Québec peninsula. *Can. J. Fish. Aquat. Sci.* **41**:
1781-1802.

Legendre, P. and L. Legendre. 2012. *Numerical ecology*, 3rd English
edition. Elsevier Science BV, Amsterdam.

Sokal, R. R. and F. J. Rohlf. 1995. *Biometry. The principles and
practice of statistics in biological research.* 3rd edn. W. H. Freeman,
New York.

## Author

Pierre Legendre, Departement de Sciences Biologiques, Universite de
Montreal

## Note

The function uses a more powerful alternative for the McNemar test than
the classical formula. The classical formula was constructed in the
spirit of Pearson's Chi-square, but the formula in this function was
constructed in the spirit of Wilks Chi-square or the \\G\\ statistic.
Function [`mcnemar.test`](https://rdrr.io/r/stats/mcnemar.test.html)
uses the classical formula. The new formula was introduced in vegan
version 1.10-11, and the older implementations of `bgdispersal` used the
classical formula.

## Examples

``` r
mat <- matrix(c(32,15,14,10,70,30,100,4,10,30,25,0,18,0,40,
  0,0,20,0,0,0,0,4,0,30,20,0,0,0,0,25,74,42,1,45,89,5,16,16,20),
  4, 10, byrow=TRUE)
bgdispersal(mat)
#> $DD1
#>       [,1]  [,2] [,3]  [,4]
#> [1,]  0.00  0.24 0.21  0.00
#> [2,] -0.24  0.00 0.08 -0.24
#> [3,] -0.21 -0.08 0.00 -0.21
#> [4,]  0.00  0.24 0.21  0.00
#> 
#> $DD2
#>            [,1]       [,2]      [,3]       [,4]
#> [1,]  0.0000000  0.3428571 0.3230769  0.0000000
#> [2,] -0.3428571  0.0000000 0.1142857 -0.3428571
#> [3,] -0.3230769 -0.1142857 0.0000000 -0.3230769
#> [4,]  0.0000000  0.3428571 0.3230769  0.0000000
#> 
#> $DD3
#>             [,1]       [,2]      [,3]        [,4]
#> [1,]  0.00000000  0.1567922 0.1420408 -0.01325831
#> [2,] -0.15679216  0.0000000 0.1101196 -0.20049485
#> [3,] -0.14204082 -0.1101196 0.0000000 -0.13586560
#> [4,]  0.01325831  0.2004949 0.1358656  0.00000000
#> 
#> $DD4
#>             [,1]       [,2]      [,3]        [,4]
#> [1,]  0.00000000  0.2513176 0.2425087 -0.01960102
#> [2,] -0.25131757  0.0000000 0.1725441 -0.30993929
#> [3,] -0.24250871 -0.1725441 0.0000000 -0.23381521
#> [4,]  0.01960102  0.3099393 0.2338152  0.00000000
#> 
#> $McNemar
#>      [,1]     [,2]      [,3]     [,4]
#> [1,]   NA 7.677938 9.0571232 0.000000
#> [2,]   NA       NA 0.2912555 7.677938
#> [3,]   NA       NA        NA 9.057123
#> [4,]   NA       NA        NA       NA
#> 
#> $prob.McNemar
#>      [,1]        [,2]        [,3]        [,4]
#> [1,]   NA 0.005590001 0.002616734 1.000000000
#> [2,]   NA          NA 0.589417103 0.005590001
#> [3,]   NA          NA          NA 0.002616734
#> [4,]   NA          NA          NA          NA
#> 
```
