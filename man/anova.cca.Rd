\name{anova.cca}
\alias{anova.cca}
\alias{anova.ccanull}
\alias{anova.ccabyaxis}
\alias{anova.ccabyterm}
\alias{anova.ccabymargin}
\alias{permutest}
\alias{permutest.default}
\alias{permutest.cca}
\alias{print.permutest.cca}

\title{Permutation Test for Constrained Correspondence Analysis,
  Redundancy Analysis and Constrained Analysis of Principal Coordinates }
\description{
  The function performs an ANOVA like permutation test for Constrained
  Correspondence Analysis (\code{\link{cca}}), Redundancy Analysis
  (\code{\link{rda}}) or Constrained Analysis of Principal Coordinates
  (\code{\link{capscale}}) to assess the significance of constraints.
}
\usage{
\method{anova}{cca}(object, alpha=0.05, beta=0.01, step=100, perm.max=9999,
      by = NULL, ...)

permutest(x, ...)

\method{permutest}{cca}(x, permutations = 100,
          model = c("reduced", "direct", "full"),
          first = FALSE, strata, ...)
}

\arguments{
  \item{object,x}{A result object from \code{\link{cca}}. }
  \item{alpha}{Targeted Type I error rate. }
  \item{beta}{Accepted Type II error rate. }
  \item{step}{Number of permutations during one step. }
  \item{perm.max}{Maximum number of permutations. }
  \item{by}{Setting \code{by = "axis"} will assess significance for each
    constrained axis, and setting \code{by = "terms"} will assess
    significance for each term (sequentially from first to last), and
    setting \code{by = "margin"} will assess the marginal effects of the
    terms (each marginal term analysed in a model with all other
    variables).}
  \item{\dots}{Parameters passed to other functions. 
    \code{anova.cca} passes all arguments to
    \code{permutest.cca}. In \code{anova} with \code{by = "axis"} you
    can use argument \code{cutoff} (defaults \code{1}) which stops
    permutations after exceeding the given level. }
  \item{permutations}{Number of permutations for assessing significance
    of constraints.}
  \item{model}{Permutation model (partial match).}
  \item{first}{Assess only the significance of the first constrained
    eigenvalue; will be passed from \code{anova.cca}. }
  \item{strata}{An integer vector or factor specifying the strata for
    permutation. If supplied, observations are permuted only within the
    specified strata.}
}
\details{
  Functions \code{anova.cca} and \code{permutest.cca} implement an ANOVA
  like permutation test for the joint effect of constraints in
  \code{\link{cca}}, \code{\link{rda}} or \code{\link{capscale}}.
  Functions \code{anova.cca} and \code{permutest.cca} differ in printout
  style and in interface.
  Function \code{permutest.cca} is the proper workhorse, but
  \code{anova.cca} passes all parameters to \code{permutest.cca}.

  In \code{anova.cca} the number of permutations is controlled by
  targeted \dQuote{critical} \eqn{P} value (\code{alpha}) and accepted Type
  II or rejection error (\code{beta}).  If the results of permutations
  differ from the targeted \code{alpha} at risk level given by
  \code{beta}, the permutations are
  terminated.  If the current estimate of \eqn{P} does not
  differ significantly from \code{alpha} of the alternative hypothesis,
  the permutations are
  continued with \code{step} new permutations (at the first step, the
  number of permutations is \code{step - 1}).  However, with \code{by =
    "terms"} a fixed number of permutations will be used, and this 
  is given by argument \code{permutations}, or if this is missing,
  by \code{step}.  
  
  The function \code{permutest.cca} implements a permutation test for
  the \dQuote{significance} of constraints in \code{\link{cca}},
  \code{\link{rda}} or \code{\link{capscale}}.  Community data are
  permuted with choice \code{model = "direct"}, residuals after
  partial CCA/RDA/CAP with choice \code{model = "reduced"} (default),
  and residuals after CCA/RDA/CAP under choice \code{model = "full"}.
  If there is no partial CCA/RDA/CAP stage, \code{model = "reduced"}
  simply permutes the data and is equivalent to \code{model = "direct"}. 
  The test statistic is ``pseudo-\eqn{F}'',
  which is the ratio of constrained and unconstrained total Inertia
  (Chi-squares, variances or something similar), each divided by their
  respective ranks.  If there are no conditions (\dQuote{partial} terms), the
  sum of all eigenvalues remains constant, so that pseudo-\eqn{F} and
  eigenvalues would give equal results.  In partial CCA/RDA/CAP, the
  effect of conditioning variables (\dQuote{covariables} is removed before
  permutation, and these residuals are added to the non-permuted fitted
  values of partial CCA (fitted values of \code{X ~ Z}).  Consequently,
  the total Chi-square is not fixed, and test based on pseudo-\eqn{F}
  would differ from the test based on plain eigenvalues. CCA is a
  weighted method, and environmental data are re-weighted at each
  permutation step using permuted weights. 

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
  performed. You can stop permutation tests after exceeding a given
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
  same residual of the full model. Permutations for all axes or terms
  will start from the same \code{\link{.Random.seed}}, and the seed
  will be advanced to the value after the longest permutation at the
  exit from the function.  
}

\value{
  Function \code{permutest.cca} returns an object of class
  \code{"permutest.cca"}, which has its own \code{print} method.  The
  function \code{anova.cca} calls \code{permutest.cca}, fills an
  \code{\link{anova}} table and uses \code{\link{print.anova}} for printing.
}
\note{
  The default permutation \code{model} changed from \code{"direct"} to
  \code{"reduced"} in \pkg{vegan} version 1.14-11 (release version
  1.15-0), and you must explicitly set \code{model = "direct"} for
  compatibility with the old version.

  Tests \code{by = "terms"} and \code{by = "margin"} are consistent
  only when \code{model = "direct"}.  
}
\references{
  Legendre, P. and Legendre, L. (1998). \emph{Numerical Ecology}. 2nd English
  ed. Elsevier.
}
\author{Jari  Oksanen}
\seealso{\code{\link{cca}}, \code{\link{rda}}, \code{\link{capscale}}
  to get something to analyse. Function \code{\link{drop1.cca}} calls
  \code{anova.cca} with \code{by = "margin"}, and
  \code{\link{add1.cca}} an analysis for single terms additions, which
  can be used in automatic or semiautomatic model building (see
  \code{\link{deviance.cca}}). }

\examples{
data(varespec)
data(varechem)
vare.cca <- cca(varespec ~ Al + P + K, varechem)
## overall test
anova(vare.cca)
## Test for axes
anova(vare.cca, by="axis", perm.max=500)
## Sequential test for terms
anova(vare.cca, by="terms", permu=200)
## Marginal or Type III effects
anova(vare.cca, by="margin")
## Marginal test knows 'scope'
anova(vare.cca, by = "m", scope="P")
}
\keyword{ multivariate }
\keyword{ htest }

