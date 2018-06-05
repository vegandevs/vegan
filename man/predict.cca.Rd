\name{predict.cca}
\alias{fitted.cca}
\alias{fitted.rda}
\alias{fitted.capscale}
\alias{fitted.dbrda}
\alias{residuals.cca}
\alias{predict.cca}
\alias{predict.rda}
\alias{predict.decorana}
\alias{coef.cca}
\alias{coef.rda}
\alias{calibrate.cca}
\alias{calibrate}

\title{Prediction Tools for [Constrained] Ordination (CCA,
  RDA, DCA, CA, PCA) }
\description{
  Function \code{predict} can be used to find site and species scores or
  estimates of the response data with new data sets, Function
  \code{calibrate} estimates values of constraints with new data set.
  Functions \code{fitted} and \code{residuals} return estimates of
  response data. 
}
\usage{
\method{fitted}{cca}(object, model = c("CCA", "CA", "pCCA"),
    type =  c("response", "working"), ...)
\method{fitted}{capscale}(object, model = c("CCA", "CA", "pCCA", "Imaginary"),
    type = c("response", "working"), ...)
\method{residuals}{cca}(object, ...)
\method{predict}{cca}(object, newdata, type = c("response", "wa", "sp", "lc", "working"),
        rank = "full", model = c("CCA", "CA"), scaling = "none",
        hill = FALSE, ...)
\method{predict}{rda}(object, newdata, type = c("response", "wa", "sp", "lc", "working"),
        rank = "full", model = c("CCA", "CA"), scaling = "none",
        correlation = FALSE, const, ...)
\method{calibrate}{cca}(object, newdata, rank = "full", ...)
\method{coef}{cca}(object, norm = FALSE, ...)
\method{predict}{decorana}(object, newdata, type = c("response", "sites", "species"),
    rank = 4, ...)
}

\arguments{

  \item{object}{A result object from \code{\link{cca}},
    \code{\link{rda}}, \code{\link{dbrda}}, \code{\link{capscale}} or
    \code{\link{decorana}}. }
  
  \item{model}{Show constrained (\code{"CCA"}), unconstrained
    (\code{"CA"}) or conditioned \dQuote{partial} (\code{"pCCA"})
    results. For \code{fitted} method of \code{\link{capscale}} this can
    also be \code{"Imaginary"} for imaginary components with negative
    eigenvalues }

  \item{newdata}{New data frame to be used in prediction or in
    calibration.  Usually this a new community data frame, but with
    \code{type = "lc"} and for constrained component with \code{type =
    "response"} and \code{type = "working"} it must be a data frame of
    constraints.  The \code{newdata} must have the same number of rows
    as the original community data for a \code{\link{cca}} result with
    \code{type = "response"} or \code{type = "working"}.  If the
    original model had row or column names, then new data must contain
    rows or columns with the same names (row names for species scores,
    column names for \code{"wa"} scores and constraint names of
    \code{"lc"} scores). In other cases the rows or columns must match
    directly. }
  \item{type}{The type of prediction, fitted values or residuals:
    \code{"response"} scales results so that the same ordination gives
    the same results, and \code{"working"} gives the values used
    internally, that is after Chi-square standardization in
    \code{\link{cca}} and scaling and centring in
    \code{\link{rda}}. In \code{\link{capscale}} and
    \code{\link{dbrda}} the \code{"response"} gives the
    dissimilarities, and \code{"working"} the internal data structure
    analysed in the ordination. Alternative \code{"wa"} gives the site
    scores as weighted averages of the community data, \code{"lc"} the
    site scores as linear combinations of environmental data, and
    \code{"sp"} the species scores. In \code{predict.decorana} the
    alternatives are scores for \code{"sites"} or \code{"species"}.}
  \item{rank}{The rank or the number of axes used in the approximation.
    The default is to use all axes (full rank) of the \code{"model"} or
    all available four axes in \code{predict.decorana}.}
  \item{scaling}{logical, character, or numeric; Scaling or predicted
    scores with the same meaning as in \code{\link{cca}},
    \code{\link{rda}}, \code{\link{dbrda}}, and
    \code{\link{capscale}}. See \code{scores.cca} for further details
    on acceptable values.}
  \item{correlation, hill}{logical; correlation-like scores or Hill's
    scaling as appropriate for RDA and CCA respectively. See
    \code{\link{scores.cca}} for additional details.}
  \item{const}{Constant multiplier for  RDA scores. This will be used
    only when \code{scaling} is not \code{FALSE}, and the default value
    will give similar scaling as in \code{\link{scores.rda}}.}
  \item{norm}{Coefficients for variables that are centred and scaled
    to unit norm.}

  \item{\dots}{Other parameters to the functions.}
}
\details{

  Function \code{fitted} gives the approximation of the original data
  matrix or dissimilarities from the ordination result either in the
  scale of the response or as scaled internally by the function.
  Function \code{residuals} gives the approximation of the original data
  from the unconstrained ordination.  With argument \code{type =
  "response"} the \code{fitted.cca} and \code{residuals.cca} function
  both give the same marginal totals as the original data matrix, and
  fitted and residuals do not add up to the original data.  Functions
  \code{fitted} and \code{residuals} for \code{\link{dbrda}} and
  \code{\link{capscale}} give the dissimilarities with \code{type =
  "response"}, but these are not additive.  However, the
  \code{"working"} scores are additive for \code{\link{capscale}} (but
  not for \code{\link{dbrda}}). The \code{fitted} and \code{residuals}
  for \code{\link{capscale}} and \code{\link{dbrda}} will include the
  additive constant if that was requested in the function call.  All
  variants of \code{fitted} and \code{residuals} are defined so that for
  model \code{mod <- cca(y ~ x)}, \code{cca(fitted(mod))} is equal to
  constrained ordination, and \code{cca(residuals(mod))} is equal to
  unconstrained part of the ordination.

  Function \code{predict} can find the estimate of the original data
  matrix or dissimilarities (\code{type = "response"}) with any rank.
  With \code{rank = "full"} it is identical to \code{fitted}.  In
  addition, the function can find the species scores or site scores from
  the community data matrix for \code{\link{cca}} or \code{\link{rda}}.
  The function can be used with new data, and it can be used to add new
  species or site scores to existing ordinations.  The function returns
  (weighted) orthonormal scores by default, and you must specify
  explicit \code{scaling} to add those scores to ordination
  diagrams. With \code{type = "wa"} the function finds the site scores
  from species scores. In that case, the new data can contain new sites,
  but species must match in the original and new data.  With \code{type="sp"} 
  the function finds species scores from site constraints
  (linear combination scores). In that case the new data can contain new
  species, but sites must match in the original and new data. With
  \code{type = "lc"} the function finds the linear combination scores
  for sites from environmental data. In that case the new data frame
  must contain all constraining and conditioning environmental variables
  of the model formula. With \code{type = "response"} or 
  \code{type = "working"} the new data must contain environmental variables 
  if constrained component is desired, and community data matrix if
  residual or unconstrained component is desired.  With these types, the
  function uses \code{newdata} to find new \code{"lc"} (constrained) or
  \code{"wa"} scores (unconstrained) and then finds the response or
  working data from these new row scores and species scores.  The
  original site (row) and species (column) weights are used for
  \code{type = "response"} and \code{type = "working"} in correspondence
  analysis (\code{\link{cca}}) and therefore the number of rows must
  match in the original data and \code{newdata}.

  If a completely new data frame is created, extreme care is needed
  defining variables similarly as in the original model, in particular
  with (ordered) factors. If ordination was performed with the formula
  interface, the \code{newdata} can be a data frame or matrix, but
  extreme care is needed that the columns match in the original and
  \code{newdata}.

  Function \code{calibrate.cca} finds estimates of constraints from
  community ordination or \code{"wa"} scores from \code{\link{cca}},
  \code{\link{rda}} and \code{\link{capscale}}. This is often known as
  calibration, bioindication or environmental reconstruction.
  Basically, the method is similar to projecting site scores onto
  biplot arrows, but it uses regression coefficients.  The function
  can be called with \code{newdata} so that cross-validation is
  possible.  The \code{newdata} may contain new sites, but species
  must match in the original and new data.  The function does not work
  with \sQuote{partial} models with \code{Condition} term, and it
  cannot be used with \code{newdata} for \code{\link{capscale}} or
  \code{\link{dbrda}} results.  The results may only be interpretable
  for continuous variables.
  
  Function \code{coef} will give the regression coefficients from centred
  environmental variables (constraints and conditions) to linear
  combination scores. The coefficients are for unstandardized environmental
  variables. The coefficients will be \code{NA} for aliased effects.

  Function \code{predict.decorana} is similar to \code{predict.cca}.
  However, \code{type = "species"} is not available in detrended
  correspondence analysis  (DCA), because detrending destroys the mutual
  reciprocal averaging (except for the first axis when rescaling is not
  used). Detrended CA does not attempt to approximate the original data
  matrix, so \code{type = "response"} has no meaning in detrended
  analysis (except with \code{rank = 1}).
}

\value{
  The functions return matrices, vectors or dissimilarities as is appropriate.
}

\references{
  Greenacre, M. J. (1984). Theory and applications of correspondence
  analysis. Academic Press, London.
}

\author{Jari Oksanen.}


\seealso{\code{\link{cca}}, \code{\link{rda}}, \code{\link{dbrda}},
  \code{\link{capscale}}, \code{\link{decorana}},
  \code{\link{goodness.cca}}. }

\examples{
data(dune)
data(dune.env)
mod <- cca(dune ~ A1 + Management + Condition(Moisture), data=dune.env)
# Definition of the concepts 'fitted' and 'residuals'
mod
cca(fitted(mod))
cca(residuals(mod))
# Remove rare species (freq==1) from 'cca' and find their scores
# 'passively'.
freq <- specnumber(dune, MARGIN=2)
freq
mod <- cca(dune[, freq>1] ~ A1 + Management + Condition(Moisture), dune.env)
predict(mod, type="sp", newdata=dune[, freq==1], scaling="species")
# New sites
predict(mod, type="lc", new=data.frame(A1 = 3, Management="NM", Moisture="2"), scal=2)
# Calibration and residual plot
mod <- cca(dune ~ A1 + Moisture, dune.env)
pred <- calibrate(mod)
pred
with(dune.env, plot(A1, pred[,"A1"] - A1, ylab="Prediction Error"))
abline(h=0)
}
\keyword{ multivariate }

