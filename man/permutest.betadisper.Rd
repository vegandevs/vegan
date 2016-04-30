\name{permutest.betadisper}
\alias{permutest.betadisper}

\title{Permutation test of multivariate homogeneity of groups dispersions
  (variances)}
\description{
  Implements a permutation-based test of multivariate homogeneity of
  group dispersions (variances) for the results of a call to
  \code{\link{betadisper}}.
}
\usage{
\method{permutest}{betadisper}(x, pairwise = FALSE,
          permutations = 999,
          parallel = getOption("mc.cores"),
          \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"betadisper"}, the result of a
    call to \code{betadisper}.}
  \item{pairwise}{logical; perform pairwise comparisons of group means?}
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}
  \item{parallel}{Number of parallel processes or a predefined socket
    cluster.  With \code{parallel = 1} uses ordinary, non-parallel
    processing.}
  \item{\dots}{Arguments passed to other methods.}
}
\details{
  To test if one or more groups is more variable than the others, ANOVA
  of the distances to group centroids can be performed and parametric
  theory used to interpret the significance of F. An alternative is to
  use a permutation test. \code{permutest.betadisper} permutes model
  residuals to generate a permutation distribution of F under the Null
  hypothesis of no difference in dispersion between groups.

  Pairwise comparisons of group mean dispersions can be performed by
  setting argument \code{pairwise} to \code{TRUE}. A classical t test
  is performed on the pairwise group dispersions. This is combined with a
  permutation test based on the t statistic calculated on pairwise group
  dispersions. An alternative to the classical comparison of group
  dispersions, is to calculate Tukey's Honest Significant Differences
  between groups, via \code{\link{TukeyHSD.betadisper}}.
}
\value{
  \code{permutest.betadisper} returns a list of class
  \code{"permutest.betadisper"} with the following components:

  \item{tab}{the ANOVA table which is an object inheriting from class
    \code{"data.frame"}.}
  \item{pairwise}{a list with components \code{observed} and
    \code{permuted} containing the observed and permuted p-values for
    pairwise comparisons of group mean distances (dispersions or variances).}
  \item{groups}{character; the levels of the grouping factor.}
  \item{control}{a list, the result of a call to
    \code{\link{how}}.}
}
\references{
  Anderson, M.J. (2006) Distance-based tests for homogeneity of
  multivariate dispersions. \emph{Biometrics} \strong{62(1)}, 245--253.

  Anderson, M.J., Ellingsen, K.E. & McArdle, B.H. (2006) Multivariate
  dispersion as a measure of beta diversity. \emph{Ecology Letters}
  \strong{9(6)}, 683--693.
}
\author{Gavin L. Simpson}
\seealso{For the main fitting function see \code{\link{betadisper}}. For
  an alternative approach to determining which groups are more variable,
  see \code{\link{TukeyHSD.betadisper}}.}
\examples{
data(varespec)

## Bray-Curtis distances between samples
dis <- vegdist(varespec)

## First 16 sites grazed, remaining 8 sites ungrazed
groups <- factor(c(rep(1,16), rep(2,8)), labels = c("grazed","ungrazed"))

## Calculate multivariate dispersions
mod <- betadisper(dis, groups)
mod

## Perform test
anova(mod)

## Permutation test for F
pmod <- permutest(mod, permutations = 99, pairwise = TRUE)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)

## Has permustats() method
pstat <- permustats(pmod)
densityplot(pstat, scales = list(x = list(relation = "free")))
qqmath(pstat, scales = list(relation = "free"))
}
\keyword{methods}
\keyword{multivariate}
