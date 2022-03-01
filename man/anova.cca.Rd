\name{anova.cca}
\alias{anova.cca}
%\alias{anova.ccanull}
%\alias{anova.ccabyaxis}
%\alias{anova.ccabyterm}
%\alias{anova.ccabymargin}
%\alias{anova.prc}
\alias{permutest}
%\alias{permutest.default}
\alias{permutest.cca}

\title{Permutation Test for Constrained Correspondence Analysis,
  Redundancy Analysis and Constrained Analysis of Principal Coordinates }

\description{
  The function performs an ANOVA like permutation test for Constrained
  Correspondence Analysis (\code{\link{cca}}), Redundancy Analysis
  (\code{\link{rda}}) or distance-based Redundancy Analysis (dbRDA,
  \code{\link{capscale}}) to assess the significance of constraints.
}

\usage{
\method{anova}{cca}(object, ..., permutations = how(nperm=999),
     by = NULL, model = c("reduced", "direct", "full"),
     parallel = getOption("mc.cores"), strata = NULL,
     cutoff = 1, scope = NULL)
\method{permutest}{cca}(x, permutations = how(nperm = 99),
     model = c("reduced", "direct", "full"), by = NULL, first = FALSE,
     strata = NULL, parallel = getOption("mc.cores"), ...)
}

\arguments{

  \item{object}{One or several result objects from \code{\link{cca}},
    \code{\link{rda}}, \code{\link{dbrda}} or \code{\link{capscale}}. If
    there are several result objects, they are compared against each
    other in the order they were supplied. For a single object, a test
    specified in \code{by} or an overall test is given.}

  \item{x}{A single ordination result object.}

  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}

  \item{by}{Setting \code{by = "axis"} will assess significance for
    each constrained axis, and setting \code{by = "terms"} will assess
    significance for each term (sequentially from first to last), and
    setting \code{by = "margin"} will assess the marginal effects of
    the terms (each marginal term analysed in a model with all other
    variables). Function \code{permutest} accepts choices
    \code{"terms"} for sequential test of terms, and \code{"onedf"}
    for sequential test of one-degree-of-freedom contrasts.}

  \item{model}{Permutation model: \code{model="direct"} permutes
    community data, \code{model="reduced"} permutes residuals
    of the community data after Conditions (partial model),
    \code{model = "full"} permutes residuals after Conditions and
    Constraints.}

  \item{parallel}{Use parallel processing with the given number of
    cores.}

  \item{strata}{An integer vector or factor specifying the strata for
    permutation. If supplied, observations are permuted only within
    the specified strata. It is an error to use this when
    \code{permutations} is a matrix, or a \code{\link[permute]{how}}
    defines \code{blocks}. This is a legacy argument that will be
    deprecated in the future: use
    \code{permutations = how(\dots, blocks)} instead. }

  \item{cutoff}{Only effective with \code{by="axis"} where stops
    permutations after an axis exceeds the \code{cutoff} \eqn{p}-value.}

  \item{scope}{Only effective with \code{by="margin"} where it can be
    used to select the marginal terms for testing. The default is to
    test all marginal terms in \code{\link{drop.scope}}.}

  \item{first}{Analyse only significance of the first axis.}

  \item{\dots}{Parameters passed to other functions.  \code{anova.cca}
    passes all arguments to \code{permutest.cca}. In \code{anova} with
    \code{by = "axis"} you can use argument \code{cutoff} (defaults
    \code{1}) which stops permutations after exceeding the given
    level. }
}

\details{

  Functions \code{anova.cca} and \code{permutest.cca} implement ANOVA
  like permutation tests for the joint effect of constraints in
  \code{\link{cca}}, \code{\link{rda}}, \code{\link{dbrda}} or
  \code{\link{capscale}}. Function \code{anova} is intended as a more
  user-friendly alternative to \code{permutest} (that is the real
  workhorse).

  Function \code{anova} can analyse a sequence of constrained
  ordination models. The analysis is based on the differences in
  residual deviance in permutations of nested models.

  The default test is for the sum of all constrained eigenvalues.
  Setting \code{first = TRUE} will perform a test for the first
  constrained eigenvalue.  Argument \code{first} can be set either in
  \code{anova.cca} or in \code{permutest.cca}.  It is also possible to
  perform significance tests for each axis or for each term
  (constraining variable) using argument \code{by} in
  \code{anova.cca}.  Setting \code{by = "axis"} will perform separate
  significance tests for each constrained axis.  All previous
  constrained axes will be used as conditions (\dQuote{partialled
  out}) and a test for the first constrained eigenvalues is
  performed (Legendre et al. 2011).
  You can stop permutation tests after exceeding a given
  significance level with argument \code{cutoff} to speed up
  calculations in large models. Setting \code{by = "terms"} will
  perform separate significance test for each term (constraining
  variable). The terms are assessed sequentially from first to last,
  and the order of the terms will influence their
  significances. Setting \code{by = "margin"} will perform separate
  significance test for each marginal term in a model with all other
  terms. The marginal test also accepts a \code{scope} argument for
  the \code{\link{drop.scope}} which can be a character vector of term
  labels that are analysed, or a fitted model of lower scope.  The
  marginal effects are also known as \dQuote{Type III} effects, but
  the current function only evaluates marginal terms. It will, for
  instance, ignore main effects that are included in interaction
  terms. In calculating pseudo-\eqn{F}, all terms are compared to the
  same residual of the full model.

  Community data are permuted with choice \code{model="direct"}, and
  residuals after partial CCA/ RDA/ dbRDA with choice
  \code{model="reduced"} (default).  If there is no partial CCA/ RDA/
  dbRDA stage, \code{model="reduced"} simply permutes the data and is
  equivalent to \code{model="direct"}.  The test statistic is
  \dQuote{pseudo-\eqn{F}}, which is the ratio of constrained and
  unconstrained total Inertia (Chi-squares, variances or something
  similar), each divided by their respective degrees of freedom.  If
  there are no conditions (\dQuote{partial} terms), the sum of all
  eigenvalues remains constant, so that pseudo-\eqn{F} and eigenvalues
  would give equal results.  In partial CCA/ RDA/ dbRDA, the effect of
  conditioning variables (\dQuote{covariables}) is removed before
  permutation, and the total Chi-square is not fixed, and test based on
  pseudo-\eqn{F} would differ from the test based on plain
  eigenvalues.

}

\value{
  The function \code{anova.cca} calls \code{permutest.cca} and fills an
  \code{\link{anova}} table. Additional attributes are
  \code{Random.seed} (the random seeds used),
  \code{control} (the permutation design, see \link[permute]{how}) and
  \code{F.perm} (the permuted test statistics).
}

\references{
  Legendre, P. and Legendre, L. (2012). \emph{Numerical Ecology}. 3rd
  English ed. Elsevier.

  Legendre, P., Oksanen, J. and ter Braak, C.J.F. (2011). Testing the
  significance of canonical axes in redundancy analysis.
  \emph{Methods in Ecology and Evolution} 2, 269--277.
}
\author{Jari Oksanen}

\seealso{\code{\link{anova.cca}}, \code{\link{cca}},
  \code{\link{rda}}, \code{\link{capscale}} to get something to
  analyse. Function \code{\link{drop1.cca}} calls \code{anova.cca}
  with \code{by = "margin"}, and \code{\link{add1.cca}} an analysis
  for single terms additions, which can be used in automatic or
  semiautomatic model building (see \code{\link{deviance.cca}}). }

\examples{
data(varespec, varechem)
mod <- cca(varespec ~ Al + P + K, varechem)
## overall test
anova(mod)
## tests for individual terms
anova(mod, by="term")
anova(mod, by="margin")
## test for adding all environmental variables
anova(mod, cca(varespec ~ ., varechem))
}
\keyword{ multivariate }
\keyword{ htest }
