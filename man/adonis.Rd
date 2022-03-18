\encoding{UTF-8}
\name{adonis}
\alias{adonis2}

\title{Permutational Multivariate Analysis of Variance Using Distance
  Matrices}

\description{Analysis of variance using distance matrices --- for
  partitioning distance matrices among sources of variation and fitting
  linear models (e.g., factors, polynomial regression) to distance 
  matrices; uses a permutation test with pseudo-\eqn{F} ratios.}

\usage{
adonis2(formula, data, permutations = 999, method = "bray",
    sqrt.dist = FALSE, add = FALSE, by = "terms",
    parallel = getOption("mc.cores"), na.action = na.fail,
    strata = NULL, ...)
}

\arguments{

  \item{formula}{Model formula. The left-hand side (LHS) of the formula
    must be either a community data matrix or a dissimilarity matrix,
    e.g., from \code{\link{vegdist}} or \code{\link{dist}}.  If the LHS
    is a data matrix, function \code{\link{vegdist}} will be used to
    find the dissimilarities. The right-hand side (RHS) of the formula
    defines the independent variables. These can be continuous variables
    or factors, they can be transformed within the formula, and they can
    have interactions as in a typical \code{\link{formula}}.}
  \item{data}{ the data frame for the independent variables.}
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}
  \item{method}{ the name of any method used in \code{\link{vegdist}} to
    calculate pairwise distances if the left hand side of the
    \code{formula} was a data frame or a matrix. }
  \item{sqrt.dist}{Take square root of dissimilarities. This often
    euclidifies dissimilarities.}
  \item{add}{Add a constant to the non-diagonal dissimilarities such
    that all eigenvalues are non-negative in the underlying Principal
    Co-ordinates Analysis (see \code{\link{wcmdscale}} for
    details). Choice \code{"lingoes"} (or \code{TRUE}) use the
    recommended method of Legendre & Anderson (1999: \dQuote{method
    1}) and \code{"cailliez"} uses their \dQuote{method 2}.}
  \item{by}{\code{by = "terms"} will assess significance for each term
    (sequentially from first to last), setting \code{by = "margin"}
    will assess the marginal effects of the terms (each marginal term
    analysed in a model with all other variables), and \code{by =
    NULL} will assess the overall significance of all terms
    together. The arguments is passed on to \code{\link{anova.cca}}.}
  \item{parallel}{Number of parallel processes or a predefined socket
    cluster.  With \code{parallel = 1} uses ordinary, non-parallel
    processing. The parallel processing is done with \pkg{parallel}
    package.}
  \item{na.action}{Handling of missing values on the right-hand-side
    of the formula (see \code{\link{na.fail}} for explanation and
    alternatives). Missing values are not allowed on the
    left-hand-side. NB, argument \code{subset} is not implemented.}

  \item{strata}{Groups within which to constrain permutations. The
    traditional non-movable strata are set as Blocks in the
    \CRANpkg{permute} package, but some more flexible alternatives may
    be more appropriate.}

  \item{\dots}{Other arguments passed to \code{\link{vegdist}}.}
}

\details{

\code{adonis2} is a function for the analysis and partitioning sums of
squares using dissimilarities. The function is based on the principles
of McArdle & Anderson (2001) and can perform sequential, marginal and
overall tests. The function also allows using additive constants or
squareroot of dissimilarities to avoid negative eigenvalues, but can
also handle semimetric indices (such as Bray-Curtis) that produce
negative eigenvalues. The \code{adonis2} tests are identical to
\code{\link{anova.cca}} of \code{\link{dbrda}}. With Euclidean
distances, the tests are also identical to \code{\link{anova.cca}} of
\code{\link{rda}}.

The function partitions sums of squares of a multivariate data set, and
they are directly analogous to MANOVA (multivariate analysis of
variance). McArdle and Anderson (2001) and Anderson (2001) refer to the
method as \dQuote{permutational MANOVA} (formerly \dQuote{nonparametric
MANOVA}). Further, as the inputs are linear predictors, and a response
matrix of an arbitrary number of columns, they are a robust alternative
to both parametric MANOVA and to ordination methods for describing how
variation is attributed to different experimental treatments or
uncontrolled covariates. The method is also analogous to distance-based
redundancy analysis in functions \code{\link{dbrda}} and
\code{\link{capscale}} (Legendre and Anderson 1999), and provides
an alternative to AMOVA (nested analysis of molecular variance,
Excoffier, Smouse, and Quattro, 1992; \code{amova} in the
\pkg{ade4} package) for both crossed and nested factors.

}

\value{

  The function returns an \code{\link{anova.cca}} result object with a
  new column for partial \eqn{R^2}{R-squared}: This is the proportion
  of sum of squares from the total, and in marginal models
  (\code{by = "margin"}) the \eqn{R^2}{R-squared} terms do not add up to
  1.

}

\note{Anderson (2001, Fig. 4) warns that the method may confound
  location and dispersion effects: significant differences may be caused
  by different within-group variation (dispersion) instead of different
  mean values of the groups (see Warton et al. 2012 for a general
  analysis). However, it seems that \code{adonis2} is less sensitive to
  dispersion effects than some of its alternatives (\code{\link{anosim}},
  \code{\link{mrpp}}). Function \code{\link{betadisper}} is a sister
  function to \code{adonis2} to study the differences in dispersion
  within the same geometric framework.
}

\references{
Anderson, M.J. 2001. A new method for non-parametric multivariate
analysis of variance. \emph{Austral Ecology}, \strong{26}: 32--46.

Excoffier, L., P.E. Smouse, and J.M. Quattro. 1992. Analysis of
molecular variance inferred from metric distances among DNA haplotypes:
Application to human mitochondrial DNA restriction data. \emph{Genetics},
\strong{131}:479--491.

Legendre, P. and M.J. Anderson. 1999. Distance-based redundancy
analysis: Testing multispecies responses in multifactorial ecological
experiments. \emph{Ecological Monographs}, \strong{69}:1--24.

McArdle, B.H.  and M.J. Anderson. 2001. Fitting multivariate models to
community data: A comment on distance-based redundancy
analysis. \emph{Ecology}, \strong{82}: 290--297.

Warton, D.I., Wright, T.W., Wang, Y. 2012. Distance-based multivariate
analyses confound location and dispersion effects. \emph{Methods in
Ecology and Evolution}, 3, 89--101.
}

\author{Martin Henry H. Stevens and Jari Oksanen.}

\seealso{ \code{\link{mrpp}}, \code{\link{anosim}},
  \code{\link{mantel}}, \code{\link{varpart}}. }
\examples{
data(dune)
data(dune.env)
## default test by terms
adonis2(dune ~ Management*A1, data = dune.env)
## overall tests
adonis2(dune ~ Management*A1, data = dune.env, by = NULL)

### Example of use with strata, for nested (e.g., block) designs.
dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
dat
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
total <- Agropyron + Schizachyrium
dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field,
        type=c('p','a'), xlab="NO3", auto.key=list(columns=3, lines=TRUE) )

Y <- data.frame(Agropyron, Schizachyrium)
mod <- metaMDS(Y, trace = FALSE)
plot(mod)
### Ellipsoid hulls show treatment
with(dat, ordiellipse(mod, field, kind = "ehull", label = TRUE))
### Spider shows fields
with(dat, ordispider(mod, field, lty=3, col="red"))

### Incorrect (no strata)
adonis2(Y ~ NO3, data = dat, permutations = 199)
## Correct with strata
with(dat, adonis2(Y ~ NO3, data = dat, permutations = 199, strata = field))
}

\keyword{multivariate }
\keyword{nonparametric }
