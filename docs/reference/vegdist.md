# Dissimilarity Indices for Community Ecologists

The function computes dissimilarity indices that are useful for or
popular with community ecologists. All indices use quantitative data,
although they would be named by the corresponding binary index, but you
can calculate the binary index using an appropriate argument. If you do
not find your favourite index here, you can see if it can be implemented
using
[`designdist`](https://vegandevs.github.io/vegan/reference/designdist.md).
Gower, Bray–Curtis, Jaccard and Kulczynski indices are good in detecting
underlying ecological gradients (Faith et al. 1987). Morisita,
Horn–Morisita, Binomial, Cao and Chao indices should be able to handle
different sample sizes (Wolda 1981, Krebs 1999, Anderson & Millar 2004),
and Mountford (1962) and Raup-Crick indices for presence–absence data
should be able to handle unknown (and variable) sample sizes. Most of
these indices are discussed by Krebs (1999) and Legendre & Legendre
(2012), and their properties further compared by Wolda (1981) and
Legendre & De Cáceres (2012). Aitchison (1986) distance is equivalent to
Euclidean distance between CLR-transformed samples (`"clr"`) and deals
with positive compositional data. Robust Aitchison distance by Martino
et al. (2019) uses robust CLR (`"rlcr"`), making it applicable to
non-negative data including zeroes (unlike the standard Aitchison).

## Usage

``` r
vegdist(x, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
        na.rm = FALSE, ...)
```

## Arguments

- x:

  Community data matrix.

- method:

  Dissimilarity index, partial match to `"manhattan"`, `"euclidean"`,
  `"canberra"`, `"clark"`, `"bray"`, `"kulczynski"`, `"jaccard"`,
  `"gower"`, `"altGower"`, `"morisita"`, `"horn"`, `"mountford"`,
  `"raup"`, `"binomial"`, `"chao"`, `"cao"`, `"mahalanobis"`, `"chisq"`,
  `"chord"`, `"hellinger"`, `"aitchison"`, or `"robust.aitchison"`.

- binary:

  Perform presence/absence standardization before analysis using
  [`decostand`](https://vegandevs.github.io/vegan/reference/decostand.md).

- diag:

  Compute diagonals.

- upper:

  Return only the upper diagonal.

- na.rm:

  Pairwise deletion of missing observations when computing
  dissimilarities (but some dissimilarities may still be `NA`, although
  calculation is handled).

- ...:

  Other parameters. These are ignored, except in `method ="gower"` which
  accepts `range.global` parameter of
  [`decostand`](https://vegandevs.github.io/vegan/reference/decostand.md),
  and in `method="aitchison"`, which accepts `pseudocount` parameter of
  [`decostand`](https://vegandevs.github.io/vegan/reference/decostand.md)
  used in the `clr` transformation.

## Details

Jaccard (`"jaccard"`), Mountford (`"mountford"`), Raup–Crick (`"raup"`),
Binomial and Chao indices are discussed later in this section. The
function also finds indices for presence/ absence data by setting
`binary = TRUE`. The following overview gives first the quantitative
version, where \\x\_{ij}\\ \\x\_{ik}\\ refer to the quantity on species
(column) \\i\\ and sites (rows) \\j\\ and \\k\\. In binary versions
\\A\\ and \\B\\ are the numbers of species on compared sites, and \\J\\
is the number of species that occur on both compared sites similarly as
in
[`designdist`](https://vegandevs.github.io/vegan/reference/designdist.md)
(many indices produce identical binary versions):

|              |                                                                                                                                 |
|--------------|---------------------------------------------------------------------------------------------------------------------------------|
| `euclidean`  | \\d\_{jk} = \sqrt{\sum_i (x\_{ij}-x\_{ik})^2}\\                                                                                 |
|              | binary: \\\sqrt{A+B-2J}\\                                                                                                       |
| `manhattan`  | \\d\_{jk}=\sum_i \|x\_{ij}-x\_{ik}\|\\                                                                                          |
|              | binary: \\A+B-2J\\                                                                                                              |
| `gower`      | \\d\_{jk} = (1/M) \sum_i \frac{\|x\_{ij}-x\_{ik}\|}{\max x_i-\min x_i}\\                                                        |
|              | binary: \\(A+B-2J)/M\\                                                                                                          |
|              | where \\M\\ is the number of columns (excluding missing values)                                                                 |
| `altGower`   | \\d\_{jk} = (1/NZ) \sum_i \|x\_{ij} - x\_{ik}\|\\                                                                               |
|              | where \\NZ\\ is the number of non-zero columns excluding double-zeros (Anderson et al. 2006).                                   |
|              | binary: \\\frac{A+B-2J}{A+B-J}\\                                                                                                |
| `canberra`   | \\d\_{jk}=\frac{1}{NZ} \sum_i \frac{\|x\_{ij}-x\_{ik}\|}{\|x\_{ij}\|+\|x\_{ik}\|}\\                                             |
|              | where \\NZ\\ is the number of non-zero entries.                                                                                 |
|              | binary: \\\frac{A+B-2J}{A+B-J}\\                                                                                                |
| `clark`      | \\d\_{jk}=\sqrt{\frac{1}{NZ} \sum_i (\frac{x\_{ij}-x\_{ik}}{x\_{ij}+x\_{ik}})^2}\\                                              |
|              | where \\NZ\\ is the number of non-zero entries.                                                                                 |
|              | binary: \\\frac{A+B-2J}{A+B-J}\\                                                                                                |
| `bray`       | \\d\_{jk} = \frac{\sum_i \|x\_{ij}-x\_{ik}\|}{\sum_i (x\_{ij}+x\_{ik})}\\                                                       |
|              | binary: \\\frac{A+B-2J}{A+B}\\                                                                                                  |
| `kulczynski` | \\d\_{jk} = 1-0.5(\frac{\sum_i \min(x\_{ij},x\_{ik})}{\sum_i x\_{ij}} + \frac{\sum_i \min(x\_{ij},x\_{ik})}{\sum_i x\_{ik}} )\\ |
|              | binary: \\1-(J/A + J/B)/2\\                                                                                                     |
| `morisita`   | \\d\_{jk} = 1 - \frac{2 \sum_i x\_{ij} x\_{ik}}{(\lambda_j + \lambda_k) \sum_i x\_{ij} \sum_i x\_{ik}}\\, where                 |
|              | \\\lambda_j = \frac{\sum_i x\_{ij} (x\_{ij} - 1)}{\sum_i x\_{ij} \sum_i (x\_{ij} - 1)}\\                                        |
|              | binary: cannot be calculated                                                                                                    |
| `horn`       | Like `morisita`, but \\\lambda_j = \sum_i x\_{ij}^2/(\sum_i x\_{ij})^2\\                                                        |
|              | binary: \\\frac{A+B-2J}{A+B}\\                                                                                                  |
| `binomial`   | \\d\_{jk} = \sum_i \[x\_{ij} \log (\frac{x\_{ij}}{n_i}) + x\_{ik} \log (\frac{x\_{ik}}{n_i}) - n_i \log(\frac{1}{2})\]/n_i\\,   |
|              | where \\n_i = x\_{ij} + x\_{ik}\\                                                                                               |
|              | binary: \\\log(2) \times (A+B-2J)\\                                                                                             |
| `cao`        | \\d\_{jk} = \frac{1}{S} \sum_i \log \left(\frac{n_i}{2}\right) - (x\_{ij} \log(x\_{ik}) + x\_{ik} \log(x\_{ij}))/n_i\\,         |
|              | where \\S\\ is the number of species in compared sites and \\n_i = x\_{ij}+x\_{ik}\\                                            |

Jaccard index is computed as \\2B/(1+B)\\, where \\B\\ is Bray–Curtis
dissimilarity.

Binomial index is derived from Binomial deviance under null hypothesis
that the two compared communities are equal. It should be able to handle
variable sample sizes. The index does not have a fixed upper limit, but
can vary among sites with no shared species. For further discussion, see
Anderson & Millar (2004).

Cao index or CYd index (Cao et al. 1997) was suggested as a minimally
biased index for high beta diversity and variable sampling intensity.
Cao index does not have a fixed upper limit, but can vary among sites
with no shared species. The index is intended for count (integer) data,
and it is undefined for zero abundances; these are replaced with
arbitrary value \\0.1\\ following Cao et al. (1997). Cao et al. (1997)
used \\\log\_{10}\\, but the current function uses natural logarithms so
that the values are approximately \\2.30\\ times higher than with
10-based logarithms. Anderson & Thompson (2004) give an alternative
formulation of Cao index to highlight its relationship with Binomial
index (above).

Mountford index is defined as \\M = 1/\alpha\\ where \\\alpha\\ is the
parameter of Fisher's logseries assuming that the compared communities
are samples from the same community (cf.
[`fisherfit`](https://vegandevs.github.io/vegan/reference/fisherfit.md),
[`fisher.alpha`](https://vegandevs.github.io/vegan/reference/diversity.md)).
The index \\M\\ is found as the positive root of equation \\\exp(aM) +
\exp(bM) = 1 + \exp\[(a+b-j)M\]\\, where \\j\\ is the number of species
occurring in both communities, and \\a\\ and \\b\\ are the number of
species in each separate community (so the index uses presence–absence
information). Mountford index is usually misrepresented in the
literature: indeed Mountford (1962) suggested an approximation to be
used as starting value in iterations, but the proper index is defined as
the root of the equation above. The function `vegdist` solves \\M\\ with
the Newton method. Please note that if either \\a\\ or \\b\\ are equal
to \\j\\, one of the communities could be a subset of other, and the
dissimilarity is \\0\\ meaning that non-identical objects may be
regarded as similar and the index is non-metric. The Mountford index is
in the range \\0 \dots \log(2)\\.

Raup–Crick dissimilarity (`method = "raup"`) is a probabilistic index
based on presence/absence data. It is defined as \\1 - prob(j)\\, or
based on the probability of observing at least \\j\\ species in shared
in compared communities. The current function uses analytic result from
hypergeometric distribution
([`phyper`](https://rdrr.io/r/stats/Hypergeometric.html)) to find the
probabilities. This probability (and the index) is dependent on the
number of species missing in both sites, and adding all-zero species to
the data or removing missing species from the data will influence the
index. The probability (and the index) may be almost zero or almost one
for a wide range of parameter values. The index is nonmetric: two
communities with no shared species may have a dissimilarity slightly
below one, and two identical communities may have dissimilarity slightly
above zero. The index uses equal occurrence probabilities for all
species, but Raup and Crick originally suggested that sampling
probabilities should be proportional to species frequencies (Chase et
al. 2011). A simulation approach with unequal species sampling
probabilities is implemented in
[`raupcrick`](https://vegandevs.github.io/vegan/reference/raupcrick.md)
function following Chase et al. (2011). The index can be also used for
transposed data to give a probabilistic dissimilarity index of species
co-occurrence (identical to Veech 2013).

Chao index tries to take into account the number of unseen species
pairs, similarly as in `method = "chao"` in
[`specpool`](https://vegandevs.github.io/vegan/reference/specpool.md).
Function `vegdist` implements a Jaccard, index defined as \\1-\frac{U
\times V}{U + V - U \times V}\\; other types can be defined with
function
[`chaodist`](https://vegandevs.github.io/vegan/reference/designdist.md).
In Chao equation, \\U = C_j/N_j + (N_k - 1)/N_k \times a_1/(2 a_2)
\times S_j/N_j\\, and \\V\\ is similar except for site index \\k\\.
\\C_j\\ is the total number of individuals in the species of site \\j\\
that are shared with site \\k\\, \\N_j\\ is the total number of
individuals at site \\j\\, \\a_1\\ (and \\a_2\\) are the number of
species occurring in site \\j\\ that have only one (or two) individuals
in site \\k\\, and \\S_j\\ is the total number of individuals in the
species present at site \\j\\ that occur with only one individual in
site \\k\\ (Chao et al. 2005).

Morisita index can be used with genuine count data (integers) only. It
is based on the idea of resampling without replacement and should not be
used with presence/absence data, but may give meaningless results if
compared sampling units (rows) have largest integer 1. Its Horn–Morisita
variant is able to handle any abundance data.

Mahalanobis distances are Euclidean distances of a matrix where columns
are centred, have unit variance, and are uncorrelated. The index is not
commonly used for community data, but it is sometimes used for
environmental variables. The calculation is based on transforming data
matrix and then using Euclidean distances following Mardia et al.
(1979). The Mahalanobis transformation usually fails when the number of
columns is larger than the number of rows (sampling units). When the
transformation fails, the distances are nearly constant except for small
numeric noise. Users must check that the returned Mahalanobis distances
are meaningful.

Euclidean and Manhattan dissimilarities are not good in gradient
separation without proper standardization but are still included for
comparison and special needs.

Chi-square distances (`"chisq"`) are Euclidean distances of Chi-square
transformed data (see
[`decostand`](https://vegandevs.github.io/vegan/reference/decostand.md)).
This is the internal standardization used in correspondence analysis
([`cca`](https://vegandevs.github.io/vegan/reference/cca.md),
[`decorana`](https://vegandevs.github.io/vegan/reference/decorana.md)).
Weighted principal coordinates analysis of these distances with row sums
as weights is equal to correspondence analysis (see the Example in
[`wcmdscale`](https://vegandevs.github.io/vegan/reference/wcmdscale.md)).
Chi-square distance is intended for non-negative data, such as typical
community data. However, it can be calculated as long as all margin sums
are positive, but warning is issued on negative data entries.

Chord distances (`"chord"`) are Euclidean distance of a matrix where
rows are standardized to unit norm (their sums of squares are 1) using
[`decostand`](https://vegandevs.github.io/vegan/reference/decostand.md).
Geometrically this standardization moves row points to a surface of
multidimensional unit sphere, and distances are the chords across the
hypersphere. Hellinger distances (`"hellinger"`) are related to Chord
distances, but data are standardized to unit total (row sums are 1)
using
[`decostand`](https://vegandevs.github.io/vegan/reference/decostand.md),
and then square root transformed. These distances have upper limit of
\\\sqrt{2}\\.

Bray–Curtis and Jaccard indices are rank-order similar, and some other
indices become identical or rank-order similar after some
standardizations, especially with presence/absence transformation of
equalizing site totals with
[`decostand`](https://vegandevs.github.io/vegan/reference/decostand.md).
Jaccard index is metric, and probably should be preferred instead of the
default Bray-Curtis which is semimetric.

Aitchison distance (1986) and robust Aitchison distance (Martino et al.
2019) are metrics that deal with compositional data. Aitchison distance
has been said to outperform Jensen-Shannon divergence and Bray-Curtis
dissimilarity, due to a better stability to subsetting and aggregation,
and it being a proper distance (Aitchison et al., 2000).

The naming conventions vary. The one adopted here is traditional rather
than truthful to priority. The function finds either quantitative or
binary variants of the indices under the same name, which correctly may
refer only to one of these alternatives For instance, the Bray index is
known also as Steinhaus, Czekanowski and Sørensen index. The
quantitative version of Jaccard should probably called Ružička index.
The abbreviation `"horn"` for the Horn–Morisita index is misleading,
since there is a separate Horn index. The abbreviation will be changed
if that index is implemented in `vegan`.

## Value

Function is a drop-in replacement for
[`dist`](https://rdrr.io/r/stats/dist.html) function and returns a
distance object of the same type. The result object adds attribute
`maxdist` that gives the theoretical maximum of the index for sampling
units that share no species, or `NA` when there is no such maximum.

## References

Aitchison, J. The Statistical Analysis of Compositional Data (1986).
London, UK: Chapman & Hall.

Aitchison, J., Barceló-Vidal, C., Martín-Fernández, J.A.,
Pawlowsky-Glahn, V. (2000). Logratio analysis and compositional
distance. *Math. Geol.* **32**, 271–275.

Anderson, M.J. and Millar, R.B. (2004). Spatial variation and effects of
habitat on temperate reef fish assemblages in northeastern New Zealand.
*Journal of Experimental Marine Biology and Ecology* 305, 191–221.

Anderson, M.J., Ellingsen, K.E. & McArdle, B.H. (2006). Multivariate
dispersion as a measure of beta diversity. *Ecology Letters* 9, 683–693.

Anderson, M.J & Thompson, A.A. (2004). Multivariate control charts for
ecological and environmental monitoring. *Ecological Applications* 14,
1921–1935.

Cao, Y., Williams, W.P. & Bark, A.W. (1997). Similarity measure bias in
river benthic Auswuchs community analysis. *Water Environment Research*
69, 95–106.

Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T. (2005). A new
statistical approach for assessing similarity of species composition
with incidence and abundance data. *Ecology Letters* 8, 148–159.

Chase, J.M., Kraft, N.J.B., Smith, K.G., Vellend, M. and Inouye, B.D.
(2011). Using null models to disentangle variation in community
dissimilarity from variation in \\\alpha\\-diversity. *Ecosphere*
2:art24 [doi:10.1890/ES10-00117.1](https://doi.org/10.1890/ES10-00117.1)

Faith, D. P, Minchin, P. R. and Belbin, L. (1987). Compositional
dissimilarity as a robust measure of ecological distance. *Vegetatio*
69, 57–68.

Gower, J. C. (1971). A general coefficient of similarity and some of its
properties. *Biometrics* 27, 623–637.

Krebs, C. J. (1999). *Ecological Methodology.* Addison Wesley Longman.

Legendre, P. & De Cáceres, M. (2012). Beta diversity as the variance of
community data: dissimilarity coefficients and partitioning. *Ecology
Letters* 16, 951–963.
[doi:10.1111/ele.12141](https://doi.org/10.1111/ele.12141)

Legendre, P. and Legendre, L. (2012) *Numerical Ecology*. 3rd English
ed. Elsevier.

Mardia, K.V., Kent, J.T. and Bibby, J.M. (1979). *Multivariate
analysis*. Academic Press.

Martino, C., Morton, J.T., Marotz, C.A., Thompson, L.R., Tripathi, A.,
Knight, R. & Zengler, K. (2019) A novel sparse compositional technique
reveals microbial perturbations. *mSystems* **4**, 1.

Mountford, M. D. (1962). An index of similarity and its application to
classification problems. In: P.W.Murphy (ed.), *Progress in Soil
Zoology*, 43–50. Butterworths.

Veech, J. A. (2013). A probabilistic model for analysing species
co-occurrence. *Global Ecology and Biogeography* 22, 252–260.

Wolda, H. (1981). Similarity indices, sample size and diversity.
*Oecologia* 50, 296–302.

## Author

Jari Oksanen, with contributions from Tyler Smith (Gower index), Michael
Bedward (Raup–Crick index), and Leo Lahti (Aitchison and robust
Aitchison distance).

## Note

The function is an alternative to
[`dist`](https://rdrr.io/r/stats/dist.html) adding some ecologically
meaningful indices. Both methods should produce similar types of objects
which can be interchanged in any method accepting either. Manhattan and
Euclidean dissimilarities should be identical in both methods. Canberra
index is divided by the number of variables in `vegdist`, but not in
[`dist`](https://rdrr.io/r/stats/dist.html). So these differ by a
constant multiplier, and the alternative in `vegdist` is in range (0,1).
Function [`daisy`](https://rdrr.io/pkg/cluster/man/daisy.html) (package
cluster) provides alternative implementation of Gower index that also
can handle mixed data of numeric and class variables. There are two
versions of Gower distance (`"gower"`, `"altGower"`) which differ in
scaling: `"gower"` divides all distances by the number of observations
(rows) and scales each column to unit range, but `"altGower"` omits
double-zeros and divides by the number of pairs with at least one
above-zero value, and does not scale columns (Anderson et al. 2006). You
can use
[`decostand`](https://vegandevs.github.io/vegan/reference/decostand.md)
to add range standardization to `"altGower"` (see Examples). Gower
(1971) suggested omitting double zeros for presences, but it is often
taken as the general feature of the Gower distances. See Examples for
implementing the Anderson et al. (2006) variant of the Gower index.

Most dissimilarity indices in `vegdist` are designed for community data,
and they will give misleading values if there are negative data entries.
The results may also be misleading or `NA` or `NaN` if there are empty
sites. In principle, you cannot study species composition without
species and you should remove empty sites from community data.

## See also

Function
[`designdist`](https://vegandevs.github.io/vegan/reference/designdist.md)
can be used for defining your own dissimilarity index. Function
[`betadiver`](https://vegandevs.github.io/vegan/reference/betadiver.md)
provides indices intended for the analysis of beta diversity.

## Examples

``` r
data(varespec)
vare.dist <- vegdist(varespec)
# Orlóci's Chord distance: range 0 .. sqrt(2)
vare.dist <- vegdist(decostand(varespec, "norm"), "euclidean")
# Anderson et al.  (2006) version of Gower
vare.dist <- vegdist(decostand(varespec, "log"), "altGower")
#> Warning: non-integer data: divided by smallest positive value
# Range standardization with "altGower" (that excludes double-zeros)
vare.dist <- vegdist(decostand(varespec, "range"), "altGower")
# Robust Aitchison distance equals to Euclidean distance for rclr transformed data
vare.dist <- vegdist(decostand(varespec, "rclr"), method = "euclidean")
vare.dist <- vegdist(varespec, "robust.aitchison")
```
