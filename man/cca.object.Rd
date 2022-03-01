\name{cca.object}
\alias{cca.object}
\alias{ordConstrained} % a common interface, not really documented
\alias{ordiYbar}
\alias{model.matrix.cca}
\alias{model.matrix.rda}
\alias{model.frame.cca}
\alias{weights.cca}
\alias{weights.rda}
\alias{weights.decorana}
\alias{print.cca}

\title{Result Object from Constrained Ordination}

\description{

  Ordination methods \code{\link{cca}}, \code{\link{rda}},
  \code{\link{dbrda}} and \code{\link{capscale}} return similar result
  objects. All these methods use the same internal function
  \code{ordConstrained}. They differ only in (1) initial
  transformation of the data and in defining inertia, (2) weighting,
  and (3) the use of rectangular rows \eqn{\times}{x} columns data or
  symmetric rows \eqn{\times}{x} rows dissimilarities:
  \code{\link{rda}} initializes data to give variance or correlations
  as inertia, \code{\link{cca}} is based on double-standardized data
  to give Chi-square inertia and uses row and column weights,
  \code{\link{capscale}} maps the real part of dissimilarities to
  rectangular data and performs RDA, and \code{\link{dbrda}} performs
  an RDA-like analysis directly on symmetric dissimilarities.

  Function \code{ordConstrained} returns the same result components
  for all these methods, and the calling function may add some more
  components to the final result. However, you should not access these
  result components directly (using \code{$}): the internal structure
  is not regarded as stable application interface (API), but it can
  change at any release. If you access the results components
  directly, you take a risk of breakage at any \pkg{vegan} release.
  The \pkg{vegan} provides a wide set of accessor functions to those
  components, and these functions are updated when the result object
  changes. This documentation gives an overview of accessor functions
  to the \code{cca} result object.

}

\usage{

ordiYbar(x, model = c("CCA", "CA", "pCCA", "partial", "initial"))
\method{model.frame}{cca}(formula, ...)
\method{model.matrix}{cca}(object, ...)
\method{weights}{cca}(object, display = "sites", ...)

}

\arguments{

  \item{object, x, formula}{A result object from \code{\link{cca}},
    \code{\link{rda}}, \code{\link{dbrda}}, or
    \code{\link{capscale}}. }

  \item{model}{Show constrained (\code{"CCA"}), unconstrained
    (\code{"CA"}) or conditioned \dQuote{partial} (\code{"pCCA"})
    results. In \code{ordiYbar} the value can also be \code{"initial"}
    for the internal working input data, and \code{"partial"} for the
    internal working input data after removing the partial effects.}

  \item{display}{Display either \code{"sites"} or \code{"species"}.}

  \item{\dots}{Other arguments passed to the the function.}

}

\details{

  The internal (\dQuote{working}) form of the dependent (community)
  data can be accessed with function \code{ordiYbar}. The form depends
  on the ordination method: for instance, in \code{\link{cca}} the
  data are weighted and Chi-square transformed, and in
  \code{\link{dbrda}} they are Gower-centred dissimilarities. The
  input data in the original (\dQuote{response}) form can be accessed
  with \code{\link{fitted.cca}} and \code{\link{residuals.cca}}.
  Function \code{\link{predict.cca}} can return either working or
  response data, and also their lower-rank approximations.

  The model matrix of independent data (\dQuote{Constraints} and
  \dQuote{Conditions}) can be extracted with \code{model.matrix}. In
  partial analysis, the function returns a list of design matrices
  called \code{Conditions} and \code{Constraints}. If either component
  was missing, a single matrix is returned. The redundant (aliased)
  terms do not appear in the model matrix. These terms can be found
  with \code{\link{alias.cca}}. Function \code{model.frame} tries to
  reconstruct the data frame from which the model matrices were
  derived. This is only possible if the original model was fitted with
  \code{formula} and \code{data} arguments, and still fails if the
  \code{data} are unavailable.

  The number of observations can be accessed with
  \code{\link{nobs.cca}}, and the residual degrees of freedom with
  \code{\link{df.residual.cca}}. The information on observations with
  missing values can be accessed with \code{\link{na.action}}.  The
  terms and formula of the fitted model can be accessed with
  \code{\link{formula}} and \code{\link{terms}}.

  The weights used in \code{\link{cca}} can be accessed with
  \code{weights}. In unweighted methods (\code{\link{rda}}) all
  weights are equal.

  The ordination results are saved in separate components for partial
  terms, constraints and residual unconstrained ordination. There is
  no guarantee that these components will have the same internal names
  as currently, and you should be cautious when developing scripts and
  functions that directly access these components.

  The constrained ordination algorithm is based on QR decomposition of
  constraints and conditions (environmental data), and the QR
  component is saved separately for partial and constrained
  components.  The QR decomposition of constraints can be accessed
  with \code{\link{qr.cca}}. This will also include the residual
  effects of partial terms (Conditions), and it should be used
  together with \code{ordiYbar(x, "partial")}. The environmental data
  are first centred in \code{rda} or weighted and centred in
  \code{cca}.  The QR decomposition is used in many functions that
  access \code{cca} results, and it can be used to find many items
  that are not directly stored in the object.  For examples, see
  \code{\link{coef.cca}}, \code{\link{coef.rda}},
  \code{\link{vif.cca}}, \code{\link{permutest.cca}},
  \code{\link{predict.cca}}, \code{\link{predict.rda}},
  \code{\link{calibrate.cca}}. See \code{\link{qr}} for other possible
  uses of this component. For instance, the rank of the constraints
  can be found from the QR decomposition.

  The eigenvalues of the solution can be accessed with
  \code{\link{eigenvals.cca}}. Eigenvalues are not evaluated for
  partial component, and they will only be available for constrained
  and residual components.

  The ordination scores are internally stored as (weighted)
  orthonormal scores matrices. These results can be accessed with
  \code{\link{scores.cca}} and \code{\link{scores.rda}} functions. The
  ordination scores are scaled when accessed with \code{\link{scores}}
  functions, but internal (weighted) orthonormal scores can be
  accessed by setting \code{scaling = FALSE}. Unconstrained residual
  component has species and site scores, and constrained component has
  also fitted site scores or linear combination scores for sites and
  biplot scores and centroids for constraint variables. The biplot
  scores correspond to the \code{model.matrix}, and centroids
  are calculated for factor variables when they were used. The scores
  can be selected by defining the axes, and there is no direct way of
  accessing all scores of a certain component. The number of dimensions
  can be assessed from \code{\link{eigenvals}}. In addition, some
  other types can be derived from the results although not saved in
  the results. For instance, regression scores and model coefficients
  can be accessed with \code{\link{scores}} and \code{\link{coef}}
  functions. Partial component will have no scores.

  Distance-based methods (\code{\link{dbrda}}, \code{\link{capscale}})
  can have negative eigenvalues and associated imaginary axis
  scores. There is no way of accessing these imaginary scores.  In
  addition, species scores are initially missing in
  \code{\link{dbrda}} and they are accessory and found after analysis
  in \code{\link{capscale}} (and may be misleading). Function
  \code{\link{sppscores}} can be used to add species scores or replace
  them with more meaningful ones.

}

\seealso{

  The core function is \code{\link{ordConstrained}} which is called by
  \code{\link{cca}}, \code{\link{rda}}, \code{\link{dbrda}} and
  \code{\link{capscale}}. The basic class is \code{"cca"} for all
  methods, and the following functions are defined for this class:
  \Sexpr[results=rd,stage=build]{require(vegan, quietly=TRUE);
  paste0("\\\\code{\\\\link{", methods(class="cca"), "}}", collapse=", ")}.
  Other functions handling \code{"cca"} objects include \code{\link{inertcomp}},
  \code{\link{intersetcor}}, \code{\link{mso}}, \code{\link{ordiresids}},
  \code{\link{ordistep}} and \code{\link{ordiR2step}}.
}

\note{

  Saving of \dQuote{working} dependent (community) data changed in
  \pkg{vegan} version 2.5-0, and you should use \code{ordiYbar}
  function instead of direct access, or your scripts and functions
  will fail (\code{ordiYbar} has been available since \pkg{vegan}
  version 2.4-3, and it works both with the old and current result
  objects).

  The \code{model.matrix} returns the unweighted model matrix also for
  \code{\link{cca}}. Prior to \pkg{vegan} version 2.5-0 it returned
  the weighted model matrix

}

\references{
  Legendre, P. and Legendre, L. (2012) \emph{Numerical Ecology}. 3rd English
  ed. Elsevier.
}
\author{ Jari Oksanen }

\keyword{ models}
\keyword{multivariate}
