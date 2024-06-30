\name{ordistep}
\alias{ordistep}
\alias{ordiR2step}

\title{
  Choose a Model by Permutation Tests in Constrained Ordination
}
\description{

  Automatic stepwise model building for constrained ordination methods
  (\code{\link{cca}}, \code{\link{rda}}, \code{\link{dbrda}},
  \code{\link{capscale}}).  The function \code{ordistep} is modelled
  after \code{\link{step}} and can do forward, backward and stepwise
  model selection using permutation tests.  Function \code{ordiR2step}
  performs forward model choice solely on adjusted \eqn{R^2}{R2} and
  \eqn{P}-value.
}

\usage{
ordistep(object, scope, direction = c("both", "backward", "forward"),
   Pin = 0.05, Pout = 0.1, permutations = how(nperm = 199), steps = 50,
   trace = TRUE, ...)
ordiR2step(object, scope, Pin = 0.05, R2scope = TRUE,
   permutations = how(nperm = 499), trace = TRUE, R2permutations = 1000, ...)
}

\arguments{
  \item{object}{
  In \code{ordistep}, an ordination object inheriting from
  \code{\link{cca}} or \code{\link{rda}}.
}

  \item{scope}{ Defines the range of models examined in the stepwise
  search.  This can be a list containing components \code{upper} and
  \code{lower}, both formulae. If it is a single item, it is interpreted
  the target scope, depending on the \code{direction}. If
  \code{direction} is \code{"forward"}, a single item is interpreted as
  the \code{upper} scope and the formula of the input \code{object} as
  the \code{lower} scope.  See \code{\link{step}} for details. In
  \code{ordiR2step}, this defines the upper scope; it can also be an
  ordination object from with the model is extracted.
  }

  \item{direction}{
  The mode of stepwise search, can be one of \code{"both"},
  \code{"backward"}, or \code{"forward"}, with a default of
  \code{"both"}.  If the \code{scope} argument is missing, the default
  for \code{direction} is \code{"backward"} in \code{ordistep} (and
  \code{ordiR2step} does not have this argument, but only works
  forward).
}
  \item{Pin, Pout}{
  Limits of permutation \eqn{P}-values for adding (\code{Pin}) a term to
  the model, or dropping (\code{Pout}) from the model. Term is added if
  \eqn{P \le}{P <=} \code{Pin}, and removed if \eqn{P >} \code{Pout}.
}

  \item{R2scope}{
  Use adjusted \eqn{R^2}{R2} as the stopping criterion: only models with
  lower adjusted \eqn{R^2}{R2} than scope are accepted.
  }

  \item{permutations}{a list of control values for the permutations as
    returned by the function \code{\link[permute]{how}}, or the number
    of permutations required, or a permutation matrix where each row
    gives the permuted indices. This is passed to
    \code{\link{anova.cca}}: see there for details.  }

  \item{steps}{
  Maximum number of iteration steps of dropping and adding terms.
}
  \item{trace}{
  If positive, information is printed during the model building. Larger
  values may give more information.
}

\item{R2permutations}{Number of permutations used in the estimation of
  adjusted \eqn{R^2}{R2} for \code{\link{cca}} using
  \code{\link{RsquareAdj}}.
}

  \item{\dots}{
  Any additional arguments to \code{\link{add1.cca}} and
  \code{\link{drop1.cca}}.
}
}
\details{
  The basic functions for model choice in constrained ordination are
  \code{\link{add1.cca}} and \code{\link{drop1.cca}}. With these functions,
  ordination models can be chosen with standard \R function
  \code{\link{step}} which bases the term choice on AIC. AIC-like
  statistics for ordination are provided by functions
  \code{\link{deviance.cca}} and \code{\link{extractAIC.cca}} (with
  similar functions for \code{\link{rda}}). Actually, constrained
  ordination methods do not have AIC, and therefore the \code{\link{step}}
  may not be trusted. This function provides an alternative using
  permutation \eqn{P}-values.

  Function \code{ordistep} defines the model, \code{scope} of models
  considered, and \code{direction} of the procedure similarly as
  \code{\link{step}}. The function alternates with \code{drop} and
  \code{add} steps and stops when the model was not changed during one
  step. The \code{-} and \code{+} signs in the summary table indicate
  which stage is performed.  It is often sensible to have \code{Pout}
  \eqn{>} \code{Pin} in stepwise models to avoid cyclic adds and drops
  of single terms.

  Function \code{ordiR2step} builds model forward so that it maximizes
  adjusted \eqn{R^2}{R2} (function \code{\link{RsquareAdj}}) at every
  step, and stopping when the adjusted \eqn{R^2}{R2} starts to decrease,
  or the adjusted \eqn{R^2}{R2} of the \code{scope} is exceeded, or the
  selected permutation \eqn{P}-value is exceeded (Blanchet et
  al. 2008). The second criterion is ignored with option \code{R2scope =
  FALSE}, and the third criterion can be ignored setting \code{Pin = 1}
  (or higher).  The function cannot be used if adjusted \eqn{R^2}{R2}
  cannot be calculated. If the number of predictors is higher than the
  number of observations, adjusted \eqn{R^2}{R2} is also unavailable.
  Such models can be analysed with \code{R2scope = FALSE}, but the
  variable selection will stop if models become overfitted and adjusted
  \eqn{R^2}{R2} cannot be calculated, and the adjusted \eqn{R^2}{R2}
  will be reported as zero. The \eqn{R^2}{R2} of \code{\link{cca}} is
  based on simulations (see \code{\link{RsquareAdj}}) and different runs
  of \code{ordiR2step} can give different results.

  Functions \code{ordistep} (based on \eqn{P} values) and \code{ordiR2step}
  (based on adjusted \eqn{R^2}{R2} and hence on eigenvalues) can select
  variables in different order.
}

\value{
  Functions return the selected model with one additional
  component, \code{anova}, which contains brief information of steps
  taken. You can suppress voluminous output during model building by
  setting \code{trace = FALSE}, and find the summary of model history
  in the \code{anova} item.
}

\references{
  Blanchet, F. G., Legendre, P. & Borcard, D. (2008) Forward selection
  of explanatory variables. \emph{Ecology} 89, 2623--2632.
}

\author{
  Jari Oksanen
}

\seealso{

  The function handles constrained ordination methods
  \code{\link{cca}}, \code{\link{rda}}, \code{\link{dbrda}} and
  \code{\link{capscale}}. The underlying functions are
  \code{\link{add1.cca}} and \code{\link{drop1.cca}}, and the function
  is modelled after standard \code{\link{step}} (which also can be
  used directly but uses AIC for model choice, see
  \code{\link{extractAIC.cca}}). Function \code{ordiR2step} builds
  upon \code{\link{RsquareAdj}}.

}
\examples{
## See add1.cca for another example

### Dune data
data(dune)
data(dune.env)
mod0 <- rda(dune ~ 1, dune.env)  # Model with intercept only
mod1 <- rda(dune ~ ., dune.env)  # Model with all explanatory variables

## With scope present, the default direction is "both"
mod <- ordistep(mod0, scope = formula(mod1))
mod
## summary table of steps
mod$anova

## Example of ordistep, forward
ordistep(mod0, scope = formula(mod1), direction="forward")

## Example of ordiR2step (always forward)
## stops because R2 of 'mod1' exceeded
ordiR2step(mod0, mod1)
}


\keyword{ multivariate }
\keyword{ models }
