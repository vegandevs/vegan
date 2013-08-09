\name{cca}
\alias{cca}
\alias{cca.default}
\alias{cca.formula}
\alias{rda}
\alias{rda.default}
\alias{rda.formula}

\title{ [Partial] [Constrained] Correspondence Analysis and Redundancy
  Analysis } 
\description{
  Function \code{cca} performs correspondence analysis, or optionally
  constrained correspondence analysis (a.k.a. canonical correspondence
  analysis), or optionally partial constrained correspondence
  analysis. Function \code{rda} performs redundancy analysis, or
  optionally principal components analysis.
  These are all very popular ordination techniques in community ecology.
}
\usage{
\method{cca}{formula}(formula, data, na.action = na.fail, subset = NULL,
  ...)
\method{cca}{default}(X, Y, Z, ...)
\method{rda}{formula}(formula, data, scale=FALSE, na.action = na.fail,
  subset = NULL, ...)
\method{rda}{default}(X, Y, Z, scale=FALSE, ...)
}

\arguments{
  \item{formula}{Model formula, where the left hand side gives the
    community data matrix, right hand side gives the constraining variables,
    and conditioning variables can be given within a special function
    \code{Condition}.}
  \item{data}{Data frame containing the variables on the right hand side
    of the model formula.}
  \item{X}{ Community data matrix. }
  \item{Y}{ Constraining matrix, typically of environmental variables.
    Can be missing. }
  \item{Z}{ Conditioning matrix, the effect of which is removed
    (`partialled out') before next step. Can be missing.}
  \item{scale}{Scale species to unit variance (like correlations).}
  \item{na.action}{Handling of missing values in constraints or
    conditions. The default (\code{\link{na.fail}}) is to stop with
    missing value. Choice \code{\link{na.omit}} removes all rows with
    missing values. Choice \code{\link{na.exclude}} keeps all
    observations but gives \code{NA} for results that cannot be
    calculated. The WA scores of rows may be found also for missing
    values in constraints. Missing values are never allowed in
    dependent community data. }
  \item{subset}{Subset of data rows. This can be a logical vector which
    is \code{TRUE} for kept observations, or a logical expression which
    can contain variables in the working environment, \code{data} or
    species names of the community data.}
  \item{...}{Other arguments for \code{print} or \code{plot} functions
    (ignored in other functions).}
}
\details{
  Since their introduction (ter Braak 1986), constrained, or canonical,
  correspondence analysis and its spin-off, redundancy analysis, have
  been the most popular ordination methods in community ecology.
  Functions \code{cca} and \code{rda} are  similar to popular
  proprietary software \code{Canoco}, although the implementation is
  completely different.  The functions are based on Legendre &
  Legendre's (2012) algorithm: in \code{cca}
  Chi-square transformed data matrix is subjected to weighted linear
  regression on constraining variables, and the fitted values are
  submitted to correspondence analysis performed via singular value
  decomposition (\code{\link{svd}}). Function \code{rda} is similar, but uses
  ordinary, unweighted linear regression and unweighted SVD. Legendre &
  Legendre (2012), Table 11.5 (p. 650) give a skeleton of the RDA
  algorithm of \pkg{vegan}. The algorithm of CCA is similar, but
  involves standardization by row and column weights.

  The functions can be called either with matrix-like entries for community
  data and constraints, or with formula interface.  In general, the
  formula interface is preferred, because it allows a better control of
  the model and allows factor constraints.

  In the following sections, \code{X}, \code{Y} and \code{Z}, although
  referred to as matrices, are more commonly data frames.

  In the matrix interface, the
  community data matrix \code{X} must be given, but the other data
  matrices may be omitted, and the corresponding stage of analysis is
  skipped.  If matrix \code{Z} is supplied, its effects are removed from
  the community matrix, and the residual matrix is submitted to the next
  stage.  This is called `partial' correspondence or redundancy
  analysis.  If matrix
  \code{Y} is supplied, it is used to constrain the ordination,
  resulting in constrained or canonical correspondence analysis, or
  redundancy analysis.
  Finally, the residual is submitted to ordinary correspondence
  analysis (or principal components analysis).  If both matrices
  \code{Z} and \code{Y} are missing, the
  data matrix is analysed by ordinary correspondence analysis (or
  principal components analysis).

  Instead of separate matrices, the model can be defined using a model
  \code{\link{formula}}.  The left hand side must be the
  community data matrix (\code{X}).  The right hand side defines the
  constraining model.
  The constraints can contain ordered or unordered factors,
  interactions among variables and functions of variables.  The defined
  \code{\link{contrasts}} are honoured in \code{\link{factor}}
  variables.  The constraints can also be matrices (but not data
  frames).
  The formula can include a special term \code{Condition}
  for conditioning variables (``covariables'') ``partialled out'' before
  analysis.  So the following commands are equivalent:
  \code{cca(X, Y, Z)},  \code{cca(X ~ Y + Condition(Z))}, where \code{Y}
  and \code{Z} refer to constraints and conditions matrices respectively.

  Constrained correspondence analysis is indeed a constrained method:
  CCA does not try to display all variation in the
  data, but only the part that can be explained by the used constraints.
  Consequently, the results are strongly dependent on the set of
  constraints and their transformations or interactions among the
  constraints.  The shotgun method is to use all environmental variables
  as constraints.  However, such exploratory problems are better
  analysed with
  unconstrained methods such as correspondence analysis
  (\code{\link{decorana}}, \code{\link[MASS]{corresp}}) or non-metric
  multidimensional scaling (\code{\link{metaMDS}}) and
  environmental interpretation after analysis
  (\code{\link{envfit}}, \code{\link{ordisurf}}).
  CCA is a good choice if the user has
  clear and strong \emph{a priori} hypotheses on constraints and is not
  interested in the major structure in the data set.  

  CCA is able to correct the curve artefact commonly found in
  correspondence analysis by forcing the configuration into linear
  constraints.  However, the curve artefact can be avoided only with a
  low number of constraints that do not have a curvilinear relation with
  each other.  The curve can reappear even with two badly chosen
  constraints or a single factor.  Although the formula interface makes
  easy to include polynomial or interaction terms, such terms often
  produce curved artefacts (that are difficult to interpret), these
  should probably be avoided.

  According to folklore, \code{rda} should be used with ``short
  gradients'' rather than \code{cca}. However, this is not based
  on research which finds methods based on Euclidean metric as uniformly
  weaker than those based on Chi-squared metric.  However, standardized
  Euclidean distance may be an appropriate measures (see Hellinger
  standardization in \code{\link{decostand}} in particular).
  
  Partial CCA (pCCA; or alternatively partial RDA) can be used to remove
  the effect of some
  conditioning or ``background'' or ``random'' variables or
  ``covariables'' before CCA proper.  In fact, pCCA compares models
  \code{cca(X ~ Z)} and \code{cca(X ~ Y + Z)} and attributes their
  difference to the effect of \code{Y} cleansed of the effect of
  \code{Z}.  Some people have used the method for extracting
  ``components of variance'' in CCA.  However, if the effect of
  variables together is stronger than sum of both separately, this can
  increase total Chi-square after ``partialling out'' some
  variation, and give negative ``components of variance''.  In general,
  such components of ``variance'' are not to be trusted due to
  interactions between two sets of variables.
  
  The functions have \code{summary} and \code{plot} methods which are
  documented separately (see \code{\link{plot.cca}}, \code{\link{summary.cca}}).

}
\value{
  Function \code{cca} returns a huge object of class \code{cca}, which
  is described separately in \code{\link{cca.object}}.

  Function \code{rda} returns an object of class \code{rda} which
  inherits from class \code{cca} and is described in \code{\link{cca.object}}.
  The scaling used in \code{rda} scores is described in a separate
  vignette with this package.
}
\references{ The original method was by ter Braak, but the current
  implementations follows Legendre and Legendre.

  Legendre, P. and Legendre, L. (2012) \emph{Numerical Ecology}. 3rd English
  ed. Elsevier.

  McCune, B. (1997) Influence of noisy environmental data on canonical
  correspondence analysis. \emph{Ecology} \strong{78}, 2617-2623.
  
  Palmer, M. W. (1993) Putting things in even better order: The
  advantages of canonical correspondence analysis.  \emph{Ecology}
  \strong{74},2215-2230. 
  
  Ter Braak, C. J. F. (1986) Canonical Correspondence Analysis: a new
  eigenvector technique for multivariate direct gradient
  analysis. \emph{Ecology} \strong{67}, 1167-1179.
  
}
\author{
  The responsible author was Jari Oksanen, but the code borrows heavily
  from Dave Roberts (\url{http://labdsv.nr.usu.edu/}).
}

\seealso{
  
  This help page describes two constrained ordination functions,
  \code{cca} and \code{rda}.  A related method, distance-based
  redundancy analysis (dbRDA) is described separately
  (\code{\link{capscale}}). All these functions return similar objects
  (described in \code{\link{cca.object}}). There are numerous support
  functions that can be used to access the result object. In the list
  below, functions of type \code{cca} will handle all three constrained
  ordination objects, and functions of \code{rda} only handle \code{rda}
  and \code{\link{capscale}} results.

  The main plotting functions are \code{\link{plot.cca}} for all
  methods, and \code{\link{biplot.rda}} for RDA and dbRDA.  However,
  generic \pkg{vegan} plotting functions can also handle the results.
  The scores can be accessed and scaled with \code{\link{scores.cca}},
  and summarized with \code{\link{summary.cca}}. The eigenvalues can
  be accessed with \code{\link{eigenvals.cca}} and the regression
  coefficients for constraints with \code{\link{coef.cca}}.  The
  eigenvalues can be plotted with \code{\link{screeplot.cca}}, and the
  (adjusted) \eqn{R^2}{R-squared} can be found with
  \code{\link{RsquareAdj.rda}}. The scores can be also calculated for
  new data sets with \code{\link{predict.cca}} which allows adding
  points to ordinations.  The values of constraints can be inferred
  from ordination and community composition with
  \code{\link{calibrate.cca}}.

  Diagnostic statistics can be found with \code{\link{goodness.cca}},
  \code{\link{inertcomp}}, \code{\link{spenvcor}},
  \code{\link{intersetcor}}, \code{\link{tolerance.cca}}, and
  \code{\link{vif.cca}}.  Function \code{\link{as.mlm.cca}} refits the
  result object as a multiple \code{\link{lm}} object, and this allows
  finding influence statistics (\code{\link{lm.influence}},
  \code{\link{cooks.distance}} etc.).
  
  Permutation based significance for the overall model, single
  constraining variables or axes can be found with
  \code{\link{anova.cca}}.  Automatic model building with \R{}
  \code{\link{step}} function is possible with
  \code{\link{deviance.cca}}, \code{\link{add1.cca}} and
  \code{\link{drop1.cca}}.  Functions \code{\link{ordistep}} and
  \code{\link{ordiR2step}} (for RDA) are special functions for
  constrained ordination. Randomized data sets can be generated with
  \code{\link{simulate.cca}}.

  Separate methods based on constrained ordination model are principal
  response curves (\code{\link{prc}}) and variance partitioning between
  several components (\code{\link{varpart}}).

  Design decisions are explained in \code{\link{vignette}}
  \file{decision-vegan} which also can be accessed with
  \code{\link{vegandocs}}.

  Package \pkg{ade4} provides alternative constrained ordination
  functions \code{\link[ade4]{cca}} and \code{\link[ade4]{pcaiv}}.

}

\examples{
data(varespec)
data(varechem)
## Common but bad way: use all variables you happen to have in your
## environmental data matrix
vare.cca <- cca(varespec, varechem)
vare.cca
plot(vare.cca)
## Formula interface and a better model
vare.cca <- cca(varespec ~ Al + P*(K + Baresoil), data=varechem)
vare.cca
plot(vare.cca)
## `Partialling out' and `negative components of variance'
cca(varespec ~ Ca, varechem)
cca(varespec ~ Ca + Condition(pH), varechem)
## RDA
data(dune)
data(dune.env)
dune.Manure <- rda(dune ~ Manure, dune.env)
plot(dune.Manure) 
## For further documentation:
\dontrun{
vegandocs("decision")
}
}
\keyword{ multivariate }

