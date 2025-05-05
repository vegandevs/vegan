\encoding{UTF-8}
\name{betadisper}
\alias{betadisper}
\alias{scores.betadisper}
\alias{anova.betadisper}
\alias{plot.betadisper}
\alias{boxplot.betadisper}
\alias{TukeyHSD.betadisper}
\alias{eigenvals.betadisper}
\alias{betadistances}
\alias{print.betadisper}
\alias{ordimedian}

\title{Multivariate homogeneity of groups dispersions (variances)}
\description{
  Implements Marti Anderson's PERMDISP2 procedure for the analysis of
  multivariate homogeneity of group dispersions (variances).
  \code{betadisper} is a multivariate analogue of Levene's test for
  homogeneity of variances. Non-euclidean distances between objects and
  group centres (centroids or medians) are handled by reducing the
  original distances to principal coordinates. This procedure has
  latterly been used as a means of assessing beta diversity. There are
  \code{anova}, \code{scores}, \code{plot} and \code{boxplot} methods.

  \code{TukeyHSD.betadisper} creates a set of confidence intervals on
  the differences between the mean distance-to-centroid of the levels of
  the grouping factor with the specified family-wise probability of
  coverage.  The intervals are based on the Studentized range statistic,
  Tukey's 'Honest Significant Difference' method.
}
\usage{
betadisper(d, group, type = c("median","centroid"), bias.adjust = FALSE,
       sqrt.dist = FALSE, add = FALSE)

\method{anova}{betadisper}(object, \dots)

\method{scores}{betadisper}(x, display = c("sites", "centroids"),
       choices = c(1,2), \dots)

\method{eigenvals}{betadisper}(x, \dots)

\method{plot}{betadisper}(x, axes = c(1,2), cex = 0.7,
     pch = seq_len(ng), col = NULL, lty = "solid", lwd = 1, hull = TRUE,
     ellipse = FALSE, conf,
     segments = TRUE, seg.col = "grey", seg.lty = lty, seg.lwd = lwd,
     label = TRUE, label.cex = 1,
     ylab, xlab, main, sub, \dots)

\method{boxplot}{betadisper}(x, ylab = "Distance to centroid", ...)

\method{TukeyHSD}{betadisper}(x, which = "group", ordered = FALSE,
         conf.level = 0.95, \dots)

\method{print}{betadisper}(x, digits = max(3, getOption("digits") - 3),
                           neigen = 8, \dots)

betadistances(x, \dots)
}

\arguments{
  \item{d}{a distance structure such as that returned by 
    \code{\link[stats]{dist}}, \code{\link{betadiver}} or
    \code{\link{vegdist}}.}
  \item{group}{vector describing the group structure, usually a factor
    or an object that can be coerced to a factor using
    \code{\link{as.factor}}. Can consist of a factor with a single
    level (i.e., one group).}
  \item{type}{the type of analysis to perform. Use the spatial median or
    the group centroid? The spatial median is now the default.}
  \item{bias.adjust}{logical: adjust for small sample bias in beta
    diversity estimates?}
  \item{sqrt.dist}{Take square root of dissimilarities. This often
    euclidifies dissimilarities.}
  \item{add}{Add a constant to the non-diagonal dissimilarities such
    that all eigenvalues are non-negative in the underlying Principal
    Co-ordinates Analysis (see \code{\link{wcmdscale}} for
    details). Choice \code{"lingoes"} (or \code{TRUE}) use the
    recommended method of Legendre & Anderson (1999: \dQuote{method
    1}) and \code{"cailliez"} uses their \dQuote{method 2}.}
  \item{display}{character; partial match to access scores for
    \code{"sites"} or \code{"species"}.}
  \item{object, x}{an object of class \code{"betadisper"}, the result of a
    call to \code{betadisper}.}
  \item{choices, axes}{the principal coordinate axes wanted.}
  \item{hull}{logical; should the convex hull for each group be plotted?}
  \item{ellipse}{logical; should the standard deviation data ellipse for
    each group be plotted?}
  \item{conf}{Expected fractions of data coverage for data ellipses,
    e.g. 0.95. The default is to draw a 1 standard deviation data
    ellipse, but if supplied, \code{conf} is multiplied with the
    corresponding value found from the Chi-squared distribution with 2df
    to provide the requested coverage (probability contour).}
  \item{pch}{plot symbols for the groups, a vector of length equal to
    the number of groups.}
  \item{col}{colors for the plot symbols and centroid labels for the groups,
    a vector of length equal to the number of groups.}
  \item{lty, lwd}{linetype, linewidth for convex hulls and confidence
    ellipses.}
  \item{segments}{logical; should segments joining points to their
    centroid be drawn?}
  \item{seg.col}{colour to draw segments between points and their
    centroid. Can be a vector, in which case one colour per group.}
  \item{seg.lty, seg.lwd}{linetype and line width for segments.}
  \item{label}{logical; should the centroids by labelled with their
    respective factor label?}
  \item{label.cex}{numeric; character expansion for centroid labels.}
  \item{cex, ylab, xlab, main, sub}{graphical parameters. For details,
    see \code{\link{plot.default}}.}
  \item{which}{A character vector listing terms in the fitted model for
    which the intervals should be calculated. Defaults to the grouping
    factor.}
  \item{ordered}{logical; see \code{\link{TukeyHSD}}.}
  \item{conf.level}{A numeric value between zero and one giving the
    family-wise confidence level to use.}
  \item{digits, neigen}{numeric; for the \code{print} method, sets the
  number of digits to use (as per \code{\link{print.default}}) and the
  maximum number of axes to display eigenvalues for, repsectively.}
  \item{\dots}{arguments, including graphical parameters (for
    \code{plot.betadisper} and \code{boxplot.betadisper}), passed to
    other methods.}
}
\details{
  One measure of multivariate dispersion (variance) for a group of
  samples is to calculate the average distance of group members to the
  group centroid or spatial median (both referred to as 'centroid' from
  now on unless stated otherwise) in multivariate space. To test if the
  dispersions (variances) of one or more groups are different, the
  distances of group members to the group centroid are subject to
  ANOVA. This is a multivariate analogue of Levene's test for
  homogeneity of variances if the distances between group members and
  group centroids is the Euclidean distance.

  However, better measures of distance than the Euclidean distance are
  available for ecological data. These can be accommodated by reducing
  the distances produced using any dissimilarity coefficient to
  principal coordinates, which embeds them within a Euclidean space. The
  analysis then proceeds by calculating the Euclidean distances between
  group members and the group centroid on the basis of the principal
  coordinate axes rather than the original distances.
  
  Non-metric dissimilarity coefficients can produce principal coordinate
  axes that have negative Eigenvalues. These correspond to the
  imaginary, non-metric part of the distance between objects. If
  negative Eigenvalues are produced, we must correct for these imaginary
  distances.

  The distance to its centroid of a point is \deqn{z_{ij}^c =
  \sqrt{\Delta^2(u_{ij}^+, c_i^+) - \Delta^2(u_{ij}^-, c_i^-)},}{z[ij]^c
  = sqrt(Delta^2(u[ij]^+, c[i]^+) - Delta^2(u[ij]^-, c[i]^-)),} where
  \eqn{\Delta^2}{Delta^2} is the squared Euclidean distance between
  \eqn{u_{ij}}{u[ij]}, the principal coordinate for the \eqn{j}th
  point in the \eqn{i}th group, and \eqn{c_i}{c[i]}, the
  coordinate of the centroid for the \eqn{i}th group. The
  super-scripted \sQuote{\eqn{+}} and \sQuote{\eqn{-}} indicate the
  real and imaginary parts respectively. This is equation (3) in
  Anderson (2006). If the imaginary part is greater in magnitude than
  the real part, then we would be taking the square root of a negative
  value, resulting in NaN, and these cases are changed to zero distances
  (with a warning). This is in line with the behaviour of Marti Anderson's
  PERMDISP2 programme. Function \code{betadistances} returns distances
  from all points to all centroids. Moreover, it gives the original
  \code{group} and nearest group for each point.
  
  To test if one or more groups is more variable than the others, ANOVA
  of the distances to group centroids can be performed and parametric
  theory used to interpret the significance of \eqn{F}. An alternative is to
  use a permutation test. \code{\link{permutest.betadisper}} permutes model
  residuals to generate a permutation distribution of \eqn{F} under the Null
  hypothesis of no difference in dispersion between groups.

  Pairwise comparisons of group mean dispersions can also be performed
  using \code{\link{permutest.betadisper}}. An alternative to the classical
  comparison of group dispersions, is to calculate Tukey's Honest
  Significant Differences between groups, via
  \code{TukeyHSD.betadisper}. This is a simple wrapper to
  \code{\link{TukeyHSD}}. The user is directed to read the help file
  for \code{\link{TukeyHSD}} before using this function. In particular,
  note the statement about using the function with 
  unbalanced designs.

  The results of the analysis can be visualised using the \code{plot}
  and \code{boxplot} methods. The distances of points to all centroids
  (\code{group}) can be found with function \code{betadistances}.

  One additional use of these functions is in assessing beta diversity
  (Anderson \emph{et al} 2006). Function \code{\link{betadiver}}
  provides some popular dissimilarity measures for this purpose.

  As noted in passing by Anderson (2006) and in a related
  context by O'Neill (2000), estimates of dispersion around a
  central location (median or centroid) that is calculated from the same data
  will be biased downward. This bias matters most when comparing diversity
  among treatments with small, unequal numbers of samples.  Setting
  \code{bias.adjust=TRUE} when using \code{betadisper} imposes a 
  \eqn{\sqrt{n/(n-1)}}{sqrt(n/(n-1))} correction (Stier et al. 2013).
}
\value{
  The \code{anova} method returns an object of class \code{"anova"}
  inheriting from class \code{"data.frame"}.

  The \code{scores} method returns a list with one or both of the
  components \code{"sites"} and \code{"centroids"}.
  
  The \code{plot} function invisibly returns an object of class
  \code{"ordiplot"}, a plotting structure which can be used by
  \code{\link{identify.ordiplot}} (to identify the points) or other
  functions in the \code{\link{ordiplot}} family. 

  The \code{boxplot} function invisibly returns a list whose components
  are documented in \code{\link{boxplot}}.

  \code{eigenvals.betadisper} returns a named vector of eigenvalues.

  \code{TukeyHSD.betadisper} returns a list. See \code{\link{TukeyHSD}}
  for further details.

  \code{betadisper} returns a list of class \code{"betadisper"} with the
  following components:

  \item{eig}{numeric; the eigenvalues of the principal coordinates
    analysis.}
  \item{vectors}{matrix; the eigenvectors of the principal coordinates
    analysis.}
  \item{distances}{numeric; the Euclidean distances in principal
    coordinate space between the samples and their respective group
    centroid or median.}
  \item{group}{factor; vector describing the group structure}
  \item{centroids}{matrix; the locations of the group centroids or
    medians on the principal coordinates.}
  \item{group.distances}{numeric; the mean distance to each group
    centroid or median.}
  \item{call}{the matched function call.}
}
\note{
  If \code{group} consists of a single level or group, then the
  \code{anova} and \code{permutest} methods are not appropriate and if
  used on such data will stop with an error.

  Missing values in either \code{d} or \code{group} will be removed
  prior to performing the analysis.
}
\section{Warning}{
  Stewart Schultz noticed that the permutation test for
  \code{type="centroid"} had the wrong type I error and was
  anti-conservative. As such, the default for \code{type} has been
  changed to \code{"median"}, which uses the spatial median as the group
  centroid. Tests suggests that the permutation test for this type of
  analysis gives the correct error rates.
}
\references{
  Anderson, M.J. (2006) Distance-based tests for homogeneity of
  multivariate dispersions. \emph{Biometrics} \strong{62}, 245--253.

  Anderson, M.J., Ellingsen, K.E. & McArdle, B.H. (2006) Multivariate
  dispersion as a measure of beta diversity. \emph{Ecology Letters}
  \strong{9}, 683--693.

  O'Neill, M.E. (2000) A Weighted Least Squares Approach to Levene's 
  Test of Homogeneity of Variance. \emph{Australian & New Zealand Journal of 
  Statistics} \strong{42}, 81-–100.

  Stier, A.C., Geange, S.W., Hanson, K.M., & Bolker, B.M. (2013) Predator 
  density and timing of arrival affect reef fish community
  assembly. \emph{Ecology} \strong{94}, 1057--1068.
}
\author{Gavin L. Simpson; bias correction by Adrian Stier and Ben Bolker.}
\seealso{\code{\link{permutest.betadisper}},
  \code{\link[stats]{anova.lm}},
  \code{\link{scores}}, \code{\link{boxplot}},
  \code{\link{TukeyHSD}}. Further measure of beta diversity
  can be found in \code{\link{betadiver}}.}
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
permutest(mod, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)

## Plot  and show the groups and distances to centroids on the
## first two PCoA axes
plot(mod)
betadistances(mod)

## with data ellipses instead of hulls
plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse

# plot with manual colour specification
my_cols <- c("#1b9e77", "#7570b3")
plot(mod, col = my_cols, pch = c(16,17), cex = 1.1)

## can also specify which axes to plot, ordering respected
plot(mod, axes = c(3,1), seg.col = "forestgreen", seg.lty = "dashed")

## Draw a boxplot of the distances to centroid for each group
boxplot(mod)

## `scores` and `eigenvals` also work
scrs <- scores(mod)
str(scrs)
head(scores(mod, 1:4, display = "sites"))
# group centroids/medians 
scores(mod, 1:4, display = "centroids")
# eigenvalues from the underlying principal coordinates analysis
eigenvals(mod) 

## try out bias correction; compare with mod3
(mod3B <- betadisper(dis, groups, type = "median", bias.adjust=TRUE))
anova(mod3B)
permutest(mod3B, permutations = 99)

## should always work for a single group
group <- factor(rep("grazed", NROW(varespec)))
(tmp <- betadisper(dis, group, type = "median"))
(tmp <- betadisper(dis, group, type = "centroid"))

## simulate missing values in 'd' and 'group'
## using spatial medians
groups[c(2,20)] <- NA
dis[c(2, 20)] <- NA
mod2 <- betadisper(dis, groups) ## messages
mod2
permutest(mod2, permutations = 99)
anova(mod2)
plot(mod2)
boxplot(mod2)
plot(TukeyHSD(mod2))

## Using group centroids
mod3 <- betadisper(dis, groups, type = "centroid")
mod3
permutest(mod3, permutations = 99)
anova(mod3)
plot(mod3)
boxplot(mod3)
plot(TukeyHSD(mod3))

}
\keyword{methods}
\keyword{multivariate}
\keyword{hplot}
