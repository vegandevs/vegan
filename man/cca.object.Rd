\name{cca.object}
\alias{cca.object}
\alias{ordConstrained} % a common interface, not really documented

\title{Result Object from Constrained Ordination with cca, rda or capscale }

\description{

  Ordination methods \code{\link{cca}}, \code{\link{rda}},
  \code{\link{dbrda}} and \code{\link{capscale}} return similar result
  objects. All these methods use the same internal function
  \code{ordConstrained}. They differ only in (1) initial transformation
  of the data and inertia, (2) weighting, and (3) the use of rectangular
  rows \eqn{\times}{x} columns data or symmetric rows \eqn{\times}{x}
  rows dissimilarities: \code{\link{rda}} initializes data to give
  variance or correlations as inertia, \code{\link{cca}} is based on
  double-standardized data to give Chi-square related inertia, and uses
  row and column weights, and \code{\link{dbrda}} performs an RDA-like
  analysis based on symmetric dissimilarities. Function
  \code{ordConstrained} returns the same result components for all these
  methods, and the calling function may set some more components for the
  final result. This help page describes the internal structure of the
  \code{cca} result object for programmers. This knowledge is not needed
  in normal use, but it is better to use support functions to access the
  result, because these should work even if the internal structure of
  the result object changes.}

\value{
  A \code{cca} object has the following elements:
  \item{call }{the function call.}
  \item{colsum, rowsum }{Column and row sums in \code{cca}.  In
    \code{rda}, item \code{colsum} contains standard deviations of
    species and \code{rowsum} is \code{NA}. }
  \item{grand.total}{Grand total of community data in \code{cca} and
    \code{NA} in \code{rda}.}
  \item{inertia}{Text used as the name of inertia.}
  \item{method}{Text used as the name of the ordination method.}
  \item{terms}{The \code{\link{terms}} component of the
    \code{\link{formula}}. This is missing if the ordination was not called
    with \code{formula}.}
  \item{terminfo}{Further information on terms with three subitems:
    \code{terms} which is like the \code{terms} component above, but
    lists conditions and constraints similarly;  \code{xlev}
    which lists the factor levels, and \code{ordered} which is
    \code{TRUE} to ordered factors.
    This is produced by \pkg{vegan} internal function
    \code{\link{ordiTerminfo}}, and it is needed in
    \code{\link{predict.cca}} with \code{newdata}.  This is missing if
    the ordination was not called with \code{formula}.}
  \item{tot.chi}{Total inertia or the sum of all eigenvalues.}
  \item{na.action}{The result of \code{\link{na.action}} if missing
    values in constraints were handled by \code{\link{na.omit}} or
    \code{\link{na.exclude}} (or \code{NULL} if there were no missing
    values). This is a vector of indices of missing value rows in the
    original data and a class of the action, usually either
    \code{"omit"} or \code{"exclude"}.}

  \item{pCCA, CCA, CA}{Actual ordination results for conditioned
    (partial), constrained and unconstrained components of the
    model. The components are called \code{pCCA}, \code{CCA} and
    \code{CA} for historic reasons, because \code{\link{cca}} was the
    first, and for long the only, constrained ordination method in
    \pkg{vegan}. This naming may change in the future and should not be
    trusted: it is safer to always use accessor functions instead of
    directly using these names in your own code. If constraints or
    conditions are not given, the corresponding components \code{CCA}
    and \code{pCCA} are \code{NULL}. If they are specified but have zero
    rank and zero eigenvalue (e.g., due to aliasing), they have a
    standard structure like described below, but the result scores have
    zero columns, but the correct number of rows. The residual component
    is never \code{NULL}, and if there is no residual variation (like in
    overdefined model), its scores have zero columns. The standard
    \code{print} command does not show \code{NULL} components, but it
    prints zeros for zeroed components. Items \code{pCCA}, \code{CCA}
    and \code{CA} contain following items:

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
    \item{\code{poseig}}{Number of positive eigenvalues in the
      ordination component. In most methods this is equal to
      \code{rank}, but \code{\link{dbrda}} can have negative eigenvalues.}
    \item{\code{tot.chi}}{Total inertia or the sum of all eigenvalues of the
      component.}
    \item{\code{imaginary.chi}, \code{imaginary.rank},
     \code{imaginary.u.eig}}{The sum, rank (number) of negative
     eigenvalues and scaled site scores for imaginary axes in
     \code{\link{capscale}}. Only in \code{CA} item and only if
     negative eigenvalues were found in \code{\link{capscale}}.}
    \item{\code{u}}{(Weighted) orthonormal site scores.  Please note that
      scaled scores are not stored in the \code{cca} object, but they
      are made when the object is accessed with functions like
      \code{\link{scores.cca}}, \code{\link{summary.cca}} or
      \code{\link{plot.cca}}, or their \code{rda} variants.   Only in
      \code{CCA} and \code{CA}.  In the \code{CCA} component these are
      the so-called linear combination scores. }
    \item{\code{v}}{(Weighted) orthonormal species scores.  If missing
      species were omitted from the analysis, this will contain
      attribute \code{\link{na.action}} that lists the omitted
      species. Only in \code{CCA} and \code{CA}. This is a zero-row
      matrix in \code{\link{dbrda}} which has no information on
      species.}
    \item{\code{wa}}{Site scores found as weighted averages (\code{cca}) or
      weighted sums (\code{rda}) of 
      \code{v} with weights \code{Xbar}, but the multiplying effect of
      eigenvalues  removed. These often are known as WA scores in
      \code{cca}. Only in  \code{CCA}.}
    \item{\code{Xbar}}{The initialized working response matrix after
      previous stages of analysis. In \code{CCA} this is after possible
      \code{pCCA} or after partialling out the effects of conditions,
      and in \code{CA} after both \code{pCCA} and \code{CCA}.}  } } }

\section{NA Action and Subset}{
  If the constraints had missing values or subsets, and \code{\link{na.action}}
  was set to  \code{\link{na.exclude}} or \code{\link{na.omit}}, the
  result will have some extra items:
  \describe{
    \item{\code{subset}}{subset evaluated as a logical vector
      (\code{TRUE} for included cases).}
    \item{\code{na.action}}{The object returned by
      \code{\link{na.action}} which is a named vector of indices of
      removed items. The class of the vector is either \code{"omit"} or
      \code{"exclude"} as set by \code{\link{na.action}}. The \code{na.action}
      is applied after \code{subset} so that the indices refer to the subset
      data.}
  \item{\code{residuals.zombie}}{A zombie vector of the length of
    number of rows in the residual ordination. \R versions before 2.13.0
    may use this vector to find the number of valid observations,
    and it is provided for their use although this is useless in \R 2.13.0
    and in \pkg{vegan}. Currently \R uses \code{\link{nobs.cca}} to find
    the number of observations.}
  \item{\code{rowsum.excluded}}{Row sums of removed observations. Only
    in \code{\link{cca}}.}
  \item{\code{CCA$wa.excluded}}{The WA scores for sites (found from
    community data) in constrained ordination if
    \code{\link{na.action}}  was \code{\link{na.exclude}} and the
    scores could be calculated. The scores cannot be found for 
    \code{\link{capscale}} and in partial ordination.}
  \item{\code{CA$u.excluded}}{Row scores for sites in unconstrained
    ordination with identical conditions as above.}
  }
}

\section{capscale and dbrda}{

  \pkg{vegan} has two functions for distance-based Redundancy
  analysis: \code{\link{capscale}} and \code{\link{dbrda}}.  Function
  \code{\link{capscale}} uses \code{\link{rda}} and returns its result
  object, but it may add some items depending on its arguments:
  
  \describe{
    \item{\code{metaMDSdist}}{The data set name if 
       \code{metaMDSdist = TRUE}.}
    \item{\code{sqrt.dist}}{Logical value, \code{TRUE} if squareroots of
       dissimilarities were used.}
    \item{\code{ac}}{Additive constant used if \code{add = TRUE}.}
    \item{\code{add}}{The adjustment method to find \code{ac}, either
       \code{"lingoes"} or \code{"cailliez"} (Legendre & Legendre
       2012).}
    \item{\code{adjust}}{Adjustment to base analysis on variance instead
      of sum of squares. For \eqn{n} observations this is
      \eqn{\sqrt{n-1}}{sqrt(n-1)} or \eqn{1} if adjustment was not
      done.}
  }

  Function \code{\link{dbrda}} does not use \code{\link{rda}} but
  provides a parallel implementation for dissimilarities. Its result
  output is very similar to \code{\link{capscale}} described above
  with the following differences:

  \describe{
  
   \item{\code{u}}{or row scores only give real axes with positive
      eigenvalues. The imaginary scores (if any) are in
      \code{imaginary.u}. The number of columns of real scores (positive
      eigenvalues) is given in item \code{poseig}. There is no
      \code{imaginary.u.eig}. The \code{ordConstrained} function returns
      all eigen vectors in one matrix, and the split is only made in the
      \code{\link{dbrda}} function that calls \code{ordConstrained}.} }
      
  }

  
\seealso{The description here provides a hacker's interface.  User
  level functions for further analysis and handling of \code{cca}
  objects are described in this section in \code{\link{cca}}. Also for
  a hacker interface, it may be better to use following low level
  functions to access the results:
  \code{\link{scores.cca}} (which also scales results), 
  \code{\link{predict.cca}} (which can also use \code{newdata}), 
  \code{\link{fitted.cca}}, \code{\link{residuals.cca}},
  \code{\link{alias.cca}}, \code{\link{coef.cca}},
  \code{\link{model.frame.cca}}, \code{\link{model.matrix.cca}},
  \code{\link{deviance.cca}}, \code{\link{eigenvals.cca}},
  \code{\link{RsquareAdj.cca}}, 
  \code{\link{weights.cca}}, \code{\link{nobs.cca}}, or \code{rda}
  variants of these functions.
  You can use \code{\link{as.mlm}} to cast a \code{cca.object} into
  result of multiple response
  linear model (\code{\link{lm}}) in order to more easily find some
  statistics (which in principle could be directly found from the
  \code{cca} object as well). 

  This section in \code{\link{cca}} gives a more complete list of
  methods to handle the constrained ordination result object.

}
\references{
  Legendre, P. and Legendre, L. (2012) \emph{Numerical Ecology}. 3rd English
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
