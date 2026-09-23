\name{add1.cca}
\alias{add1.cca}
\alias{drop1.cca}

\title{Add or Drop  Single Terms to a Constrained Ordination Model }
\description{
Compute all single terms that can be added to or dropped from a
constrained ordination model.
}
\usage{
\method{add1}{cca}(object, scope, test = c("none", "permutation", "F"),
    permutations = how(nperm=199), ...)
\method{drop1}{cca}(object, scope, test = c("none", "permutation", "F"), 
    permutations = how(nperm=199), ...)
}

\arguments{
  \item{object}{ A constrained ordination object from
  \code{\link{cca}}, \code{\link{rda}}, \code{\link{dbrda}} or
    \code{\link{capscale}}. }
  \item{scope}{ A formula giving the terms to be considered for adding
    or dropping; see \code{\link{add1}} for details.}
  
  \item{test}{ Should a \code{"permutation"} or parametric \code{"F"}
    test be added using \code{\link{anova.cca}}. \code{test = "F"} can
    be only used with \code{\link{rda}} models. }
 
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}

  \item{\dots}{Other arguments passed to \code{\link{add1.default}},
  \code{\link{drop1.default}}, and \code{\link{anova.cca}}.}
}

\details{

  With argument \code{test = "none"} the functions will only call
  \code{\link{add1.default}} or \code{\link{drop1.default}}. With other
  choice the functions will add test results from
  \code{\link{anova.cca}}, either a permutation test or a parameteric
  \eqn{F}-test (only with \code{\link{rda}} models). Function
  \code{drop1.cca} will call \code{\link{anova.cca}} with argument
  \code{by = "margin"}.  Function \code{add1.cca} will implement a test
  for single term additions that is not directly available in
  \code{\link{anova.cca}}.

  Functions are used implicitly in \code{\link{step}},
  \code{\link{ordiR2step}} and \code{\link{ordistep}}. The
  \code{\link{deviance.cca}} and \code{\link{deviance.rda}} used in
  \code{\link{step}} have no firm basis, and setting argument \code{test
  = "permutation"} may help in getting useful insight into validity of
  model building. Function \code{\link{ordistep}} calls alternately
  \code{drop1.cca} and \code{add1.cca} with argument 
  \code{test = "permutation"} and selects variables by their permutation
  \eqn{P}-values.  Meticulous use of \code{add1.cca} and
  \code{drop1.cca} will allow more judicious model building.

  The default number of \code{permutations} is set to a low value, because
  permutation tests can take a long time. It should be sufficient to
  give a impression on the significances of the terms, but higher
  values of \code{permutations} should be used if \eqn{P} values really
  are important.  

}

\value{
  Returns a similar object as \code{\link{add1}} and \code{\link{drop1}}.
}

\author{Jari Oksanen}

\seealso{ \code{\link{add1}}, \code{\link{drop1}} and
  \code{\link{anova.cca}} for basic methods. Functions
  \code{\link{deviance.cca}} and \code{\link{extractAIC.cca}} are used
  to produce the other arguments than test results in the
  output. \code{add1} and \code{drop1} are the basic functions in
  automatic model building with \code{\link{step}},
  \code{\link{ordistep}} and \code{\link{ordiR2step}}.}

\examples{
data(dune, dune.env)

## Automatic model building based on AIC but with permutation tests
step(cca(dune ~  1, dune.env), formula(cca(dune ~ ., dune.env)),
  test="permutation")
## vegan ordistep does the same, but based on P-values
ordistep(cca(dune ~  1, dune.env), cca(dune ~ ., dune.env))

## Manual model building using parametric F test
## define the maximal model for scope
mbig <- rda(dune ~  ., dune.env)
## define an empty model to start with
m0 <- rda(dune ~ 1, dune.env)
## manual selection and updating
add1(m0, scope=formula(mbig), test="F")
m0 <- update(m0, . ~ . + Management)
add1(m0, scope=formula(mbig), test="F")
m0 <- update(m0, . ~ . + Moisture)
## included variables still significant?
drop1(m0, test="F")
add1(m0, scope=formula(mbig), test="F")
}
\keyword{multivariate }
\keyword{ models}
