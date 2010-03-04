\name{ordistep}
\Rdversion{1.1}
\alias{ordistep}

\title{
  Choose a Model by Permutation Tests in Constrained Ordination 
}
\description{
  Automatic stepwise model building for constrained ordination methods 
  (\code{\link{cca}}, \code{\link{rda}}, \code{\link{capscale}}). 
  The function is modelled after \code{\link{step}} and can do forward,
  backward and stepwise model selection.
}
\usage{
ordistep(object, scope, direction = c("both", "backward", "forward"), Pin = 0.05, Pout = 0.1, pstep = 100, perm.max = 1000, steps = 50, trace = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An ordination object inheriting from \code{\link{cca}}.
}
  \item{scope}{
  Defines the range of models examined in the stepwise search. 
  This should be either a single formula, or a list containing
  components \code{upper} and \code{lower}, both formulae. 
  See \code{\link{step}} for details. 
}
  \item{direction}{
  The mode of stepwise search, can be one of \code{"both"},
  \code{"backward"}, or \code{"forward"}, with a default of \code{"both"}.  
  If the \code{scope} argument is missing the default for \code{direction}
  is \code{"backward"}.	  
}
  \item{Pin, Pout}{
  Limits of permutation \eqn{P}-values for adding (\code{Pin}) a term to 
  the model, or dropping (\code{Pout}) from the model. Term is added if
  \eqn{P \le}{P <=} \code{Pin}, and removed if \eqn{P >} \code{Pout}. 
}
  \item{pstep}{
  Number of permutations in one step. See \code{\link{add1.cca}}.
}
  \item{perm.max}{
  Maximum number of permutation in \code{\link{anova.cca}}.
}
  \item{steps}{
  Maximum number of iteration steps of dropping and adding terms.
}
  \item{trace}{
  If positive, information is printed during the model building. Larger
  values may give more information.
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
  step. The \code{-} and \code{+} signs in the summary
  table indicate which stage is performed. The number of permutations
  is selected adaptively with respect to the defined decision limit. It
  is often sensible to have \code{Pout} \eqn{>} \code{Pin} in stepwise
  models to avoid cyclic adds and drops of single terms. 
}

\value{ 
  Function returns the selected model with one additional
  component, \code{anova}, which contains brief information of steps
  taken. You can suppress voluminous output during model building by
  setting \code{trace = FALSE}, and find the summary of model history
  in the \code{anova} item.  
}

\author{
  Jari Oksanen
}

\seealso{
  The function handles constrained ordination methods \code{\link{cca}},
  \code{\link{rda}} and \code{\link{capscale}}. The underlying functions 
  are \code{\link{add1.cca}} and \code{\link{drop1.cca}}, and the 
  function is modelled after standard \code{\link{step}} (which also can 
  be used directly but uses AIC for model choice, see 
  \code{\link{extractAIC.cca}}). 
}
\examples{
## See add1.cca for another example
data(dune)
data(dune.env)
mod1 <- rda(dune ~ ., dune.env)
ordistep(mod1)
ordistep(rda(dune ~ 1, dune.env), scope = formula(mod1))
}

\keyword{ multivariate }
\keyword{ models }
