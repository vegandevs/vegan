\encoding{UTF-8}
\name{kendall.global}
\alias{kendall.global}
\alias{kendall.post}

\title{ Kendall coefficient of concordance }

\description{ 
   Function \code{kendall.global} computes and tests the coefficient of
  concordance among several judges (variables, species) through a
  permutation test.

  Function \code{kendall.post} carries out \emph{a posteriori} tests
  of the contributions of individual judges (variables, species) to
  the overall concordance of their group through permutation tests.

  If several groups of judges are identified in the data table,
  coefficients of concordance (\code{kendall.global}) or a posteriori
  tests (\code{kendall.post}) will be computed for each group
  separately. Use in ecology: to identify significant species
  associations.
}
\usage{
kendall.global(Y, group, nperm = 999, mult = "holm")
kendall.post(Y, group, nperm = 999, mult = "holm")
}
\arguments{
  \item{Y}{ Data file (data frame or matrix) containing quantitative or
    semiquantitative data. Rows are objects and columns are judges
    (variables). In community ecology, that table is often a
    site-by-species table. }
  \item{group}{ A vector defining how judges should be divided into
    groups. See example below. If groups are not explicitly defined,
    all judges in the data file will be considered as forming a single
    group. }
  \item{nperm}{ Number of permutations to be performed. Default is
    999. }
  \item{mult}{Correct P-values for multiple testing using the
    alternatives described in \code{\link{p.adjust}} and in addition
    \code{"sidak"} (see Details). The Bonferroni correction is overly
    conservative; it is not recommended. It is included to allow
    comparisons with the other methods.}
}
\details{
  \code{Y} must contain quantitative data. They will be transformed to
  ranks within each column before computation of the coefficient of
  concordance.

  The search for species associations described in Legendre (2005)
  proceeds in 3 steps:

  (1) Correlation analysis of the species. A possible method is to
  compute Ward's agglomerative clustering of a matrix of correlations
  among the species. In detail: (1.1) compute a Pearson or Spearman
  correlation matrix (\code{correl.matrix}) among the species; (1.2)
  turn it into a distance matrix: \code{mat.D =
  as.dist(1-correl.matrix)}; (1.3) carry out Ward's hierarchical
  clustering of that matrix using \code{hclust}: \code{clust.ward =
  hclust(mat.D, "ward")}; (1.4) plot the dendrogram:
  \code{plot(clust.ward, hang=-1)}; (1.5) cut the dendrogram in two
  groups, retrieve the vector of species membership: \code{group.2 =
  cutree(clust.ward, k=2)}. (1.6) After steps 2 and 3 below, you may
  have to come back and try divisions of the species into k = 3, 4, 5,
  \dots groups.

  (2) Compute global tests of significance of the 2 (or more) groups
  using the function \code{kendall.global} and the vector defining the
  groups. Groups that are not globally significant must be refined or
  abandoned.

  (3) Compute a posteriori tests of the contribution of individual
  species to the concordance of their group using the function
  \code{kendall.post} and the vector defining the groups. If some
  species have negative values for "Spearman.mean", this means that
  these species clearly do not belong to the group, hence that group
  is too inclusive. Go back to (1.5) and cut the dendrogram more
  finely. The left and right groups can be cut separately,
  independently of the levels along the dendrogram; write your own
  vector of group membership if \code{cutree} does not produce the
  desired groups.

  The corrections used for multiple testing are applied to the list of
  P-values (P); they take into account the number of tests (k) carried
  out simultaneously (number of groups in \code{kendall.global}, or
  number of species in \code{kendall.post}). The corrections are
  performed using function \code{\link{p.adjust}}; see that function
  for the description of the correction methods. In addition, there is
  \enc{Šidák}{Sidak} correction which defined as 
  \eqn{P_{corr} = 1 -(1 - P)^k}.
}
\value{

  A table containing the following information in rows. The columns
  correspond to the groups of "judges" defined in vector "group". When
  function \code{Kendall.post} is used, there are as many tables as
  the number of predefined groups.

  \item{W }{Kendall's coefficient of concordance, W. }

  \item{F }{F statistic. F = W*(m-1)/(1-W) where m is the number of
  judges. }

  \item{Prob.F }{Probability associated with the F statistic, computed
  from the F distribution with nu1 = n-1-(2/m) and nu2 = nu1*(m-1); n is
  the number of objects. }

  \item{Corrected prob.F }{Probabilities associated with F, corrected
  using the method selected in parameter \code{mult}. Shown only if
  there are more than one group. }

  \item{Chi2 }{Friedman's chi-square statistic (Friedman 1937) used in
  the permutation test of W. }

  \item{Prob.perm }{Permutational probabilities, uncorrected. }

  \item{Corrected prob.perm }{Permutational probabilities corrected
  using the method selected in parameter \code{mult}. Shown only if
  there are more than one group. }

  \item{Spearman.mean }{Mean of the Spearman correlations between the
  judge under test and all the other judges in the same group. }

  \item{W.per.species }{Contribution of the judge under test to the
  overall concordance statistic for that group. }
}
\references{ 

Friedman, M. 1937. The use of ranks to avoid the assumption of normality
implicit in the analysis of variance. Journal of the American
Statistical Association 32: 675-701.

Kendall, M. G. and B. Babington Smith. 1939. The problem of m
rankings. Annals of Mathematical Statistics 10: 275-287.

Legendre, P. 2005. Species associations: the Kendall coefficient of
concordance revisited. Journal of Agricultural, Biological, and
Environmental Statistics 10: 226-245.

Legendre, P. 2009. Coefficient of concordance. In: Encyclopedia of
Research Design. SAGE Publications (in press).

Siegel, S. and N. J. Castellan, Jr. 1988. Nonparametric statistics for
the behavioral sciences. 2nd edition. McGraw-Hill, New York.
}

\seealso{\code{\link{cor}}, \code{\link{friedman.test}},
\code{\link{hclust}}, \code{\link{cutree}}, \code{\link{kmeans}},
\code{\link[vegan]{cascadeKM}}, \code{\link[labdsv]{indval}}}

\author{ F. Guillaume Blanchet, University of Alberta, and Pierre
  Legendre, Université de Montréal }

\examples{
data(mite)
mite.hel <- decostand(mite, "hel")

# Reproduce the results shown in Table 2 of Legendre (2005), a single group
mite.small <- mite.hel[c(4,9,14,22,31,34,45,53,61,69),c(13:15,23)]
kendall.global(mite.small, nperm=99)
kendall.post(mite.small, mult="holm", nperm=99)

# Reproduce the results shown in Tables 3 and 4 of Legendre (2005), 2 groups
group <-c(1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,2,1,2,1,1,1,1,2,1,2,1,1,1,1,1,2,2,2,2,2)
kendall.global(mite.hel, group=group, nperm=99)
kendall.post(mite.hel, group=group, mult="holm", nperm=99)

# NOTE: 'nperm' argument usually needs to be larger than 99.
# It was set to this low value for demonstration purposes.
}

\keyword{ multivariate }
\keyword{ nonparametric }
