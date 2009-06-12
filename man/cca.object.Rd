\name{cca.object}
\alias{cca.object}

\title{Result Object from Constrained Ordination with cca, rda or capscale }
\description{
  Ordination methods \code{\link{cca}}, \code{\link{rda}} and
  \code{\link{capscale}} return similar result objects.  Function
  \code{capscale} \code{\link{inherits}} from \code{rda} and \code{rda}
  inherits from \code{cca}.  This inheritance structure is due to
  historic reasons: \code{cca} was the first of these implemented in
  vegan. Hence the nomenclature in \code{cca.object} reflects
  \code{cca}.  This help page describes the internal structure of the
  \code{cca} object for programmers.
}

\value{
  A \code{cca} object has the following elements:
  \item{call }{the function call.}
  \item{colsum, rowsum }{Column and row sums in \code{cca}.  In
    \code{rda}, item \code{colsum} contains standard deviations of
    species and \code{rowsum} is \code{NA}.}
  \item{grand.total}{Grand total of community data in \code{cca} and
    \code{NA} in \code{rda}.}
  \item{inertia}{Text used as the name of inertia.}
  \item{method}{Text used as the name of the ordination method.}
  \item{terms}{The \code{\link{terms}} component of the
    \code{\link{formula}}. This is missing if the ordination was not called
    with \code{formula}.}
  \item{terminfo}{Further information on terms with three subitems:
    \code{terms} which is like the \code{terms} component above, but
    lists conditions and constrainst similarly;  \code{xlev}
    which lists the factor levels, and \code{ordered} which is
    \code{TRUE} to ordered factors.
    This is produced by \pkg{vegan} internal function
    \code{\link{ordiTerminfo}}, and it is needed in
    \code{\link{predict.cca}} with \code{newdata}.  This is missing if
    the ordination was not called with \code{formula}.}
  \item{tot.chi}{Total inertia or the sum of all eigenvalues.}
  \item{pCCA, CCA, CA}{Actual ordination results for conditioned
    (partial), constrained and unconstrained components of the
    model. Any of these can be \code{NULL} if there is no corresponding
    component.
    Items \code{pCCA}, \code{CCA} and \code{CA} have similar
    structure, and contain following items:
    \describe{
    \item{\code{alias}}{The names of the aliased constraints or conditions.
      Function \code{\link{alias.cca}} does not access this item
      directly, but it finds the aliased variables and their defining
      equations from the \code{QR} item.}
    \item{\code{biplot}}{Biplot scores of constraints.  Only in \code{CCA}.}
    \item{\code{centroids}}{(Weighted) centroids of factor levels of
      constraints. Only in \code{CCA}. Missing if the ordination was not
    called with \code{formula}.}
    \item{\code{eig}}{Eigenvalues of axes. In \code{CCA} and \code{CA}.}
    \item{\code{envcentre}}{(Weighted) means of the original constraining or
      conditioning variables. In \code{pCCA} and in \code{CCA}.}
    \item{\code{Fit}}{The fitted values of standardized data matrix after
      fitting conditions. Only in \code{pCCA}.}
    \item{\code{QR}}{The QR decomposition of explanatory variables as produced
      by \code{\link{qr}}. 
      The constrained ordination 
      algorithm is based on QR decomposition of constraints and
      conditions (environmental data).  The environmental data
      are first centred in \code{rda} or weighted and centred in
      \code{cca}.  The QR decomposition is used in many functions that
      access \code{cca} results, and it can be used to find many items
      that are not directly stored in the object.  For examples, see
      \code{\link{coef.cca}}, \code{\link{coef.rda}},
      \code{\link{vif.cca}}, \code{\link{permutest.cca}},
      \code{\link{predict.cca}}, \code{\link{predict.rda}},
      \code{\link{calibrate.cca}}.  For possible uses of this component,
      see \code{\link{qr}}. In \code{pCCA} and \code{CCA}.} 
    \item{\code{rank}}{The rank of the ordination component.}
    \item{\code{qrank}}{The rank of the constraints which is the
     difference of the ranks of QR decompositions in \code{pCCA} and
     \code{CCA} components. Only in \code{CCA}.}
    \item{\code{tot.chi}}{Total inertia or the sum of all eigenvalues of the
      component.}
    \item{\code{imaginary.chi}, \code{imaginary.rank}}{The sum and
     rank (number) of negative eigenvalues in
     \code{\link{capscale}}. Only in \code{CA} and only if negative
     eigenvalues were found in \code{\link{capscale}}.}
    \item{\code{u}}{(Weighted) orthonormal site scores.  Please note that
      scaled scores are not stored in the \code{cca} object, but they
      are made when the object is accessed with functions like
      \code{\link{scores.cca}}, \code{\link{summary.cca}} or
      \code{\link{plot.cca}}, or their \code{rda} variants.   Only in
      \code{CCA} and \code{CA}.  In the \code{CCA} component these are
      the so-called linear combination scores. }
    \item{\code{u.eig}}{\code{u} scaled by eigenvalues.  There is no
      guarantee that any \code{.eig} variants of scores will be kept in
      the future releases.}
    \item{\code{v}}{(Weighted) orthonormal species scores.  If missing species
      were omitted from the analysis, this will contain
      attribute \code{\link{na.action}} that lists the
      omitted species. Only in \code{CCA} and \code{CA}.}
    \item{\code{v.eig}}{\code{v} weighted by eigenvalues.}
    \item{\code{wa}}{Site scores found as weighted averages (\code{cca}) or
      weighted sums (\code{rda}) of 
      \code{v} with weights \code{Xbar}, but the multiplying effect of
      eigenvalues  removed. These often are known as WA scores in
      \code{cca}. Only in  \code{CCA}.}
    \item{\code{wa.eig}}{The direct result of weighted avaraging or weighted
      summation  (matrix multiplication)
      with the resulting eigenvalue inflation.}
    \item{\code{Xbar}}{The standardized data matrix after previous stages of
      analysis. In \code{CCA} this is after possible \code{pCCA} or
      after partialling out the effects of conditions, and in \code{CA}
      after both \code{pCCA} and \code{CCA}. In \code{\link{cca}} the
      standardization is Chi-square, and in \code{\link{rda}} centring
      and optional scaling by species standard deviations using function
      \code{\link{scale}}.} }
  }
}
  
\seealso{The description here provides a hacker's interface.  For more
  user friendly acces to the \code{cca} object see
  \code{\link{alias.cca}}, \code{\link{coef.cca}},
  \code{\link{deviance.cca}}, \code{\link{predict.cca}},
  \code{\link{scores.cca}}, 
  \code{\link{summary.cca}},  \code{\link{vif.cca}},
  \code{\link{weights.cca}}, \code{\link{spenvcor}} or \code{rda}
  variants of these functions.
  You can use \code{\link{as.mlm}} to cast a \code{cca.object} into
  result of multiple response
  linear model (\code{\link{lm}}) in order to more easily find some
  statistics (which in principle could be directly found from the
  \code{cca.object} as well). 
}
\references{
  Legendre, P. and Legendre, L. (1998) \emph{Numerical Ecology}. 2nd English
  ed. Elsevier.
}
\author{ Jari Oksanen }
\examples{
# Some species will be missing in the analysis, because only a subset
# of sites is used below.
data(dune)
data(dune.env)
mod <- cca(dune[1:15,] ~ ., dune.env[1:15,])
# Look at the names of missing species
attr(mod$CCA$v, "na.action")
# Look at the names of the aliased variables:
mod$CCA$alias
# Access directly constrained weighted orthonormal species and site
# scores, constrained eigenvalues and margin sums.
spec <- mod$CCA$v
sites <- mod$CCA$u
eig <- mod$CCA$eig
rsum <- mod$rowsum
csum <- mod$colsum
}
\keyword{ models}
\keyword{multivariate}
