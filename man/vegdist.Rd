\encoding{UTF-8}
\name{vegdist}
\alias{vegdist}
\title{Dissimilarity Indices for Community Ecologists }
\description{
  The function computes dissimilarity indices that are useful for or
  popular with community ecologists. All indices use quantitative data,
  although they would be named by the corresponding binary index, but
  you can calculate the binary index using an appropriate argument.  If
  you do not find your favourite index here, you can see if it can be
  implemented using \code{\link{designdist}}.  Gower, Bray--Curtis,
  Jaccard and Kulczynski indices are good in detecting underlying
  ecological gradients (Faith et al. 1987). Morisita, Horn--Morisita,
  Binomial, Cao and Chao indices should be able to handle different
  sample sizes (Wolda 1981, Krebs 1999, Anderson & Millar 2004), and
  Mountford (1962) and Raup-Crick indices for presence--absence data
  should be able to handle unknown (and variable) sample sizes. Most of
  these indices are discussed by Krebs (1999) and Legendre & Legendre
  (2012), and their properties further compared by Wolda (1981) and
  Legendre & De \enc{Cáceres}{Caceres} (2012).
}

\usage{vegdist(x, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
        na.rm = FALSE, ...) } 
\arguments{
  \item{x}{Community data matrix.}
  \item{method}{Dissimilarity index, partial match to
    \code{"manhattan"}, \code{"euclidean"}, \code{"canberra"},
    \code{"clark"}, \code{"bray"}, \code{"kulczynski"},
    \code{"jaccard"}, \code{"gower"}, \code{"altGower"},
    \code{"morisita"}, \code{"horn"}, \code{"mountford"}, \code{"raup"},
    \code{"binomial"}, \code{"chao"}, \code{"cao"} or \code{"mahalanobis"}.}
  \item{binary}{Perform presence/absence standardization before analysis
    using \code{\link{decostand}}.}
  \item{diag}{Compute diagonals. }
  \item{upper}{Return only the upper diagonal. }
  \item{na.rm}{Pairwise deletion of missing observations when
    computing dissimilarities.}
  \item{\dots}{Other parameters.  These are ignored, except in
    \code{method ="gower"} which accepts \code{range.global} parameter of
    \code{\link{decostand}}. .}
}

\details{Jaccard (\code{"jaccard"}), Mountford (\code{"mountford"}),
  Raup--Crick (\code{"raup"}), Binomial and Chao indices are discussed
  later in this section.  The function also finds indices for presence/
  absence data by setting \code{binary = TRUE}. The following overview
  gives first the quantitative version, where \eqn{x_{ij}}{x[ij]}
  \eqn{x_{ik}}{x[ik]} refer to the quantity on species (column) \eqn{i}
  and sites (rows) \eqn{j} and \eqn{k}. In binary versions \eqn{A} and
  \eqn{B} are the numbers of species on compared sites, and \eqn{J} is
  the number of species that occur on both compared sites similarly as
  in \code{\link{designdist}} (many indices produce identical binary
  versions):
  
  \tabular{ll}{
    \code{euclidean}
    \tab \eqn{d_{jk} = \sqrt{\sum_i (x_{ij}-x_{ik})^2}}{d[jk] = sqrt(sum(x[ij]-x[ik])^2)}
    \cr \tab binary: \eqn{\sqrt{A+B-2J}}{sqrt(A+B-2*J)}
    \cr
    \code{manhattan}
    \tab \eqn{d_{jk}=\sum_i |x_{ij}-x_{ik}|}{d[jk] = sum(abs(x[ij] - x[ik]))}
    \cr \tab binary: \eqn{A+B-2J}{A+B-2*J}
    \cr
    \code{gower}
    \tab \eqn{d_{jk} = (1/M) \sum_i \frac{|x_{ij}-x_{ik}|}{\max x_i-\min
	x_i}}{d[jk] = (1/M) sum(abs(x[ij]-x[ik])/(max(x[i])-min(x[i])))}
    \cr \tab binary: \eqn{(A+B-2J)/M}{(A+B-2*J)/M}
    \cr
    \tab where \eqn{M} is the number of columns (excluding missing
    values)
    \cr
    \code{altGower}
    \tab \eqn{d_{jk} = (1/NZ) \sum_i |x_{ij} - x_{ik}|}{d[jk] = (1/NZ) sum(abs(x[ij] - x[ik]))}
    \cr
    \tab where \eqn{NZ} is the number of non-zero columns excluding
    double-zeros (Anderson et al. 2006).
    \cr \tab binary: \eqn{\frac{A+B-2J}{A+B-J}}{(A+B-2*J)/(A+B-J)}
    \cr
    \code{canberra}
    \tab \eqn{d_{jk}=\frac{1}{NZ} \sum_i
      \frac{|x_{ij}-x_{ik}|}{|x_{ij}|+|x_{ik}|}}{d[jk] = (1/NZ) sum (abs(x[ij]-x[ik])/(abs(x[ij])+abs(x[ik])))}
    \cr
    \tab where \eqn{NZ} is the number of non-zero entries.
    \cr \tab binary: \eqn{\frac{A+B-2J}{A+B-J}}{(A+B-2*J)/(A+B-J)}
    \cr
    \code{clark}
    \tab \eqn{d_{jk}=\sqrt{\frac{1}{NZ} \sum_i
      (\frac{x_{ij}-x_{ik}}{x_{ij}+x_{ik}})^2}}{d[jk] = sqrt( (1/NZ) sum (((x[ij]-x[ik])/(x[ij]+x[ik]))^2))}
    \cr
    \tab where \eqn{NZ} is the number of non-zero entries.
    \cr \tab binary: \eqn{\frac{A+B-2J}{A+B-J}}{(A+B-2*J)/(A+B-J)}
    \cr
    \code{bray}
    \tab \eqn{d_{jk} = \frac{\sum_i |x_{ij}-x_{ik}|}{\sum_i (x_{ij}+x_{ik})}}{d[jk] = (sum abs(x[ij]-x[ik]))/(sum (x[ij]+x[ik]))}
    \cr \tab binary: \eqn{\frac{A+B-2J}{A+B}}{(A+B-2*J)/(A+B)}
    \cr
    \code{kulczynski}
    \tab \eqn{d_{jk} = 1-0.5(\frac{\sum_i \min(x_{ij},x_{ik})}{\sum_i x_{ij}} +
      \frac{\sum_i \min(x_{ij},x_{ik})}{\sum_i x_{ik}} )}{d[jk] 1 - 0.5*(sum(min(x[ij],x[ik]))/(sum x[ij]) + sum(
      min(x[ij],x[ik]))/(sum x[ik]))}
    \cr \tab binary: \eqn{1-(J/A + J/B)/2}{1-(J/A + J/B)/2}
    \cr
    \code{morisita}
    \tab \eqn{d_{jk} =  1 - \frac{2 \sum_i x_{ij} x_{ik}}{(\lambda_j +
	  \lambda_k) \sum_i x_{ij} \sum_i
	  x_{ik}}}{d[jk] =  1 - 2*sum(x[ij]*x[ik])/((lambda[j]+lambda[k]) *
	sum(x[ij])*sum(x[ik]))}, where  
    \cr
    \tab \eqn{\lambda_j = \frac{\sum_i x_{ij} (x_{ij} - 1)}{\sum_i
	x_{ij} \sum_i (x_{ij} - 1)}}{lambda[j] = sum(x[ij]*(x[ij]-1))/sum(x[ij])*sum(x[ij]-1)}
    \cr \tab binary: cannot be calculated
    \cr
    \code{horn}
    \tab Like \code{morisita}, but \eqn{\lambda_j = \sum_i
      x_{ij}^2/(\sum_i x_{ij})^2}{lambda[j] = sum(x[ij]^2)/(sum(x[ij])^2)}
    \cr \tab binary: \eqn{\frac{A+B-2J}{A+B}}{(A+B-2*J)/(A+B)}
    \cr
    \code{binomial}
    \tab \eqn{d_{jk} = \sum_i [x_{ij} \log (\frac{x_{ij}}{n_i}) + x_{ik} \log
      (\frac{x_{ik}}{n_i}) - n_i \log(\frac{1}{2})]/n_i}{d[jk] =
      sum(x[ij]*log(x[ij]/n[i]) + x[ik]*log(x[ik]/n[i]) -
      n[i]*log(1/2))/n[i]},
    \cr
    \tab where \eqn{n_i = x_{ij} + x_{ik}}{n[i] = x[ij] + x[ik]}
    \cr \tab binary: \eqn{\log(2) \times (A+B-2J)}{log(2)*(A+B-2*J)}
    \cr
    \code{cao}
    \tab \eqn{d_{jk} = \frac{1}{S} \sum_i \log
    \left(\frac{n_i}{2}\right) - (x_{ij} \log(x_{ik}) + x_{ik}
    \log(x_{ij}))/n_i}{d[jk] = (1/S) * sum(log(n[i]/2) -
    (x[ij]*log(x[ik]) + x[ik]*log(x[ij]))/n[i])},
  \cr
  \tab where \eqn{S} is the number of species in compared sites and
    \eqn{n_i = x_{ij}+x_{ik}}{n[i] = x[ij] + x[ik]}
  }

  Jaccard index is computed as \eqn{2B/(1+B)}, where \eqn{B} is
  Bray--Curtis dissimilarity.

  Binomial index is derived from Binomial deviance under null hypothesis
  that the two compared communities are equal. It should be able to
  handle variable sample sizes. The index does not have a fixed upper
  limit, but can vary among sites with no shared species. For further
  discussion, see Anderson & Millar (2004).

  Cao index or CYd index (Cao et al. 1997) was suggested as a minimally
  biased index for high beta diversity and variable sampling intensity.
  Cao index does not have a fixed upper limit, but can vary among sites
  with no shared species.  The index is intended for count (integer)
  data, and it is undefined for zero abundances; these are replaced with
  arbitrary value \eqn{0.1} following Cao et al. (1997).  Cao et
  al. (1997) used \eqn{\log_{10}}{log10}, but the current function uses
  natural logarithms so that the values are approximately \eqn{2.30}
  times higher than with 10-based logarithms. Anderson & Thompson (2004)
  give an alternative formulation of Cao index to highlight its
  relationship with Binomial index (above).
  
  Mountford index is defined as \eqn{M = 1/\alpha} where \eqn{\alpha}
  is the parameter of Fisher's logseries assuming that the compared
  communities are samples from the same community
  (cf. \code{\link{fisherfit}}, \code{\link{fisher.alpha}}). The index
  \eqn{M} is found as the positive root of equation \eqn{\exp(aM) +
  \exp(bM) = 1 + \exp[(a+b-j)M]}{exp(a*M) + exp(b*M) = 1 +
  exp((a+b-j)*M)}, where \eqn{j} is the number of species occurring in
  both communities, and \eqn{a} and \eqn{b} are the number of species
  in each separate community (so the index uses presence--absence
  information). Mountford index is usually misrepresented in the
  literature: indeed Mountford (1962) suggested an approximation to be
  used as starting value in iterations, but the proper index is
  defined as the root of the equation above. The function
  \code{vegdist} solves \eqn{M} with the Newton method. Please note
  that if either \eqn{a} or \eqn{b} are equal to \eqn{j}, one of the
  communities could be a subset of other, and the dissimilarity is
  \eqn{0} meaning that non-identical objects may be regarded as
  similar and the index is non-metric. The Mountford index is in the
  range \eqn{0 \dots \log(2)}{0 \dots log(2)}, but the dissimilarities
  are divided by \eqn{\log(2)}{log(2)} so that the results will be in
  the conventional range \eqn{0 \dots 1}.

  Raup--Crick dissimilarity (\code{method = "raup"}) is a probabilistic
  index based on presence/absence data.  It is defined as \eqn{1 -
  prob(j)}, or based on the probability of observing at least \eqn{j}
  species in shared in compared communities.  The current function uses
  analytic result from hypergeometric distribution
  (\code{\link{phyper}}) to find the probabilities.  This probability
  (and the index) is dependent on the number of species missing in both
  sites, and adding all-zero species to the data or removing missing
  species from the data will influence the index.  The probability (and
  the index) may be almost zero or almost one for a wide range of
  parameter values.  The index is nonmetric: two communities with no
  shared species may have a dissimilarity slightly below one, and two
  identical communities may have dissimilarity slightly above zero. The
  index uses equal occurrence probabilities for all species, but Raup
  and Crick originally suggested that sampling probabilities should be
  proportional to species frequencies (Chase et al. 2011). A simulation
  approach with unequal species sampling probabilities is implemented in
  \code{\link{raupcrick}} function following Chase et al. (2011).  The
  index can be also used for transposed data to give a probabilistic
  dissimilarity index of species co-occurrence (identical to Veech
  2013).
  
  Chao index tries to take into account the number of unseen species
  pairs, similarly as in \code{method = "chao"} in
  \code{\link{specpool}}. Function \code{vegdist} implements a
  Jaccard, index defined as
  \eqn{1-\frac{U \times V}{U + V - U \times V}}{1 - U*V/(U + V - U*V)};
  other types can be defined with function \code{\link{chaodist}}. In Chao
  equation, \eqn{U = C_j/N_j + (N_k - 1)/N_k \times a_1/(2 a_2) \times
  S_j/N_j}{U = C[j]/N[j] + (N[k] -1)/N[k] * a1/(2*a2) * S[j]/N[j]},
  and \eqn{V} is similar except for site index
  \eqn{k}. \eqn{C_j}{C[j]} is the total number of individuals in the
  species of site \eqn{j} that are shared with site \eqn{k},
  \eqn{N_j}{N[j]} is the total number of individuals at site \eqn{j},
  \eqn{a_1}{a1} (and \eqn{a_2}{a2}) are the number of species
  occurring in site \eqn{j} that have only one (or two) individuals in
  site \eqn{k}, and \eqn{S_j}{S[j]} is the total number of individuals
  in the species present at site \eqn{j} that occur with only one
  individual in site \eqn{k} (Chao et al. 2005).

  Morisita index can be used with genuine count data (integers) only. Its
  Horn--Morisita variant is able to handle any abundance data.

  Mahalanobis distances are Euclidean distances of a matrix where
  columns are centred, have unit variance, and are uncorrelated.  The
  index is not commonly used for community data, but it is sometimes
  used for environmental variables. The calculation is based on
  transforming data matrix and then using Euclidean distances
  following Mardia et al. (1979). The Mahalanobis transformation
  usually fails when the number of columns is larger than the number
  of rows (sampling units). When the transformation fails, the
  distances are nearly constant except for small numeric noise. Users
  must check that the returned Mahalanobis distances are meaningful.

  Euclidean and Manhattan dissimilarities are not good in gradient
  separation without proper standardization but are still included for
  comparison and special needs.

  Bray--Curtis and Jaccard indices are rank-order similar, and some
  other indices become identical or rank-order similar after some 
  standardizations, especially with presence/absence transformation of
  equalizing site totals with \code{\link{decostand}}. Jaccard index is
  metric, and probably should be preferred instead of the default
  Bray-Curtis which is semimetric. 

  The naming conventions vary. The one adopted here is traditional
  rather than truthful to priority. The function finds either
  quantitative or binary variants of the indices under the same name,
  which correctly may refer only to one of these alternatives For
  instance, the Bray
  index is known also as Steinhaus, Czekanowski and \enc{Sørensen}{Sorensen} index.
  The quantitative version of Jaccard should probably called
  \enc{Ružička}{Ruzicka} index.
  The abbreviation \code{"horn"} for the Horn--Morisita index is
  misleading, since there is a separate Horn index. The abbreviation
  will be changed if that index is implemented in \code{vegan}. 
}
\value{
  Should provide a drop-in replacement for \code{\link{dist}} and
  return a distance object of the same type. 
}
\references{
  Anderson, M.J. and Millar, R.B. (2004). Spatial variation and effects
  of habitat on temperate reef fish assemblages in northeastern New
  Zealand.  \emph{Journal of Experimental Marine Biology and Ecology}
  305, 191--221.

  Anderson, M.J., Ellingsen, K.E. & McArdle, B.H. (2006). Multivariate
  dispersion as a measure of beta diversity. \emph{Ecology Letters} 
  9, 683--693.

  Anderson, M.J & Thompson, A.A. (2004). Multivariate control charts for
  ecological and environmental monitoring. \emph{Ecological
    Applications} 14, 1921--1935.

  Cao, Y., Williams, W.P. & Bark, A.W. (1997). Similarity measure bias
  in river benthic Auswuchs community analysis. \emph{Water
  Environment Research} 69, 95--106.

  Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T. (2005). A new
  statistical approach for assessing similarity of species composition
  with incidence and abundance data. \emph{Ecology Letters} 8, 148--159.

  Chase, J.M., Kraft, N.J.B., Smith, K.G., Vellend, M. and Inouye,
  B.D. (2011). Using null models to disentangle variation in community
  dissimilarity from variation in \eqn{\alpha}{alpha}-diversity.
  \emph{Ecosphere} 2:art24 \doi{10.1890/ES10-00117.1}
   
  Faith, D. P, Minchin, P. R. and Belbin, L. (1987).
  Compositional dissimilarity as a robust measure of ecological
  distance. \emph{Vegetatio} 69, 57--68.

  Gower, J. C. (1971). A general coefficient of similarity and some
  of its properties. \emph{Biometrics} 27, 623--637.

  Krebs, C. J. (1999). \emph{Ecological Methodology.} Addison Wesley
  Longman.

  Legendre, P. & De \enc{Cáceres}{Caceres}, M. (2012). Beta diversity as
  the variance of community data: dissimilarity coefficients and
  partitioning. \emph{Ecology Letters} 16, 951--963.
  \doi{10.1111/ele.12141}

  Legendre, P. and Legendre, L. (2012) \emph{Numerical Ecology}. 3rd English
  ed. Elsevier.

  Mardia, K.V., Kent, J.T. and Bibby, J.M. (1979). \emph{Multivariate analysis}.
  Academic Press.

  Mountford, M. D. (1962). An index of similarity and its application to
  classification problems. In: P.W.Murphy (ed.),
  \emph{Progress in Soil Zoology}, 43--50. Butterworths.

  Veech, J. A. (2013). A probabilistic model for analysing species
  co-occurrence. \emph{Global Ecology and Biogeography} 22, 252--260. 

  Wolda, H. (1981). Similarity indices, sample size and
  diversity. \emph{Oecologia} 50, 296--302.
}

\author{ Jari Oksanen, with contributions from Tyler Smith (Gower index)
  and Michael Bedward (Raup--Crick index). }

\note{The function is an alternative to \code{\link{dist}} adding some
  ecologically meaningful indices.  Both methods should produce similar
  types of objects which can be interchanged in any method accepting
  either.  Manhattan and Euclidean dissimilarities should be identical
  in both methods. Canberra index is divided by the number of variables
  in \code{vegdist}, but not in \code{\link{dist}}.  So these differ by
  a constant multiplier, and the alternative in \code{vegdist} is in
  range (0,1).  Function \code{\link[cluster]{daisy}} (package
  \pkg{cluster}) provides alternative implementation of Gower index that
  also can handle mixed data of numeric and class variables.  There are
  two versions of Gower distance (\code{"gower"}, \code{"altGower"})
  which differ in scaling: \code{"gower"} divides all distances by the
  number of observations (rows) and scales each column to unit range,
  but \code{"altGower"} omits double-zeros and divides by the number of
  pairs with at least one above-zero value, and does not scale columns
  (Anderson et al. 2006).  You can use \code{\link{decostand}} to add
  range standardization to \code{"altGower"} (see Examples). Gower
  (1971) suggested omitting double zeros for presences, but it is often
  taken as the general feature of the Gower distances. See Examples for
  implementing the Anderson et al. (2006) variant of the Gower index.

  Most dissimilarity indices in \code{vegdist} are designed for
  community data, and they will give misleading values if there are
  negative data entries.  The results may also be misleading or
  \code{NA} or \code{NaN} if there are empty sites.  In principle, you
  cannot study species composition without species and you should remove
  empty sites from community data.
}

\seealso{
 Function \code{\link{designdist}} can be used for defining your own
 dissimilarity index. Alternative dissimilarity functions include
 \code{\link{dist}} in base \R,
  \code{\link[cluster]{daisy}} (package \pkg{cluster}), and
  \code{\link[labdsv]{dsvdis}} (package \pkg{labdsv}).  Function
  \code{\link{betadiver}} provides indices intended for the analysis of
  beta diversity.}

\examples{
data(varespec)
vare.dist <- vegdist(varespec)
# Orlóci's Chord distance: range 0 .. sqrt(2)
vare.dist <- vegdist(decostand(varespec, "norm"), "euclidean")
# Anderson et al.  (2006) version of Gower
vare.dist <- vegdist(decostand(varespec, "log"), "altGower")
# Range standardization with "altGower" (that excludes double-zeros)
vare.dist <- vegdist(decostand(varespec, "range"), "altGower")
}
\keyword{ multivariate }
