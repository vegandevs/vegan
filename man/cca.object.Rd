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
  \item{colsum, rowsum, rowsum.excluded }{Column and row sums in
    \code{cca}.  In \code{rda}, item \code{colsum} contains standard
    deviations of species and \code{rowsum} is \code{NA}. If some data
    were removed in \code{\link{na.action}}, the row sums of excluded
    observations are in item \code{rowsum.excluded} in \code{cca} (but
    not in \code{rda}). The \code{rowsum.excluded} add to the total
    (one) of \code{rowsum}. }
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
    model. If constraints or conditions are not given, the
    corresponding components \code{CCA} and \code{pCCA} are
    \code{NULL}. If they are specified but have zero rank and zero
    eigenvalue (e.g., due to aliasing), they have a standard structure
    like described below, but the result scores have zero columns, but
    the correct number of rows. The residual component is never
    \code{NULL}, and if there is no residual variation (like in
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
    \item{\code{tot.chi}}{Total inertia or the sum of all eigenvalues of the
      component.}
    \item{\code{real.tot.chi}}{If there are negative eigenvalues in
      \code{\link{capscale}}, these will be included in \code{tot.chi},
      and the sum of positive eigenvalues will be given in these items.}
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
    
    \item{\code{v}}{(Weighted) orthonormal species scores.  If missing species
      were omitted from the analysis, this will contain
      attribute \code{\link{na.action}} that lists the
      omitted species. Only in \code{CCA} and \code{CA}.}
    
    \item{\code{wa}}{Site scores found as weighted averages (\code{cca}) or
      weighted sums (\code{rda}) of 
      \code{v} with weights \code{Xbar}, but the multiplying effect of
      eigenvalues  removed. These often are known as WA scores in
      \code{cca}. Only in  \code{CCA}.}

    \item{\code{wa.excluded, u.excluded}}{WA scores for rows removed by
      \code{na.action = na.exclude} in \code{CCA} and \code{CA}
      components if these could be calculated.}
    \item{\code{Xbar}}{The standardized data matrix after previous stages of
      analysis. In \code{CCA} this is after possible \code{pCCA} or
      after partialling out the effects of conditions, and in \code{CA}
      after both \code{pCCA} and \code{CCA}. In \code{\link{cca}} the
      standardization is Chi-square, and in \code{\link{rda}} centring
      and optional scaling by species standard deviations using function
      \code{\link{scale}}.}
  }
}
}

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

  \pkg{Vegan} has two functions for distance-based Redundancy
  analysis: \code{\link{capscale}} and \code{\link{dbrda}}.  Function
  \code{\link{capscale}} uses \code{\link{rda}} and returns its result
  object, but it may add some items depending on its arguments:
  
  \describe{
    \item{\code{real.tot.chi}}{Sum of positive eigenvalues if there are
       negative eigenvalues. The item \code{tot.chi} gives the total
       inertia with negative eigenvalues. This item is given for the
       whole model and for each component \code{pCCA}, \code{CCA} and
       \code{CA} if there are negative eigenvalues.}
    \item{\code{metaMDSdist}}{The data set name if 
       \code{metaMDSdist = TRUE}.}
    \item{\code{sqrt.dist}}{Logical value, \code{TRUE} if squareroots of
       dissimilarities were used.}
    \item{\code{ac}}{Additive constant used if \code{add = TRUE}.}
    \item{\code{add}}{The adjustment method to find \code{ac}, either
       \code{"lingoes"} or \code{"cailliez"} (Legendre & Legendre
       2012).}
    \item{\code{adjust}}{Adjustment of dissimilarities: see
      \code{\link{capscale}}, section \dQuote{Notes}.}
    \item{\code{G}}{The working structure of Gower transformed
      dissimilarities defined as \eqn{-(D^2 - M)/2}, where \eqn{D} are
      the dissimilarities and \eqn{M} is the centring matrix. This
      structure is used to asses the total inertia, and it will be used
      also in permutation tests. This is given for items \code{pCCA} and
      \code{CCA}, and for \code{CCA} it is the residual \eqn{G} after
      \code{pCCA}.}
  }

  Function \code{\link{dbrda}} does not use \code{\link{rda}} but
  provides a parallel implementation for dissimilarities. Its result
  output is very similar to \code{\link{capscale}} described above
  with the following differences:

  \describe{
  
    \item{\code{Xbar}, \code{v}}{are \code{NA} because they cannot be
      calculated from dissimilarities.}

    \item{\code{Fit}}{ of \code{pCCA} is from Gower double centred
      dissimilarities \code{G} instead of \code{Xbar} (that does not
      exist).}
      
    \item{\code{G}}{ is returned with \code{pCCA}, \code{CCA} and
      \code{CA} components. It always gives the transformed
      dissimilarities as they enter the stage of analysis, i.e.,
      before applying conditions or constraints.}

    \item{\code{eig}}{lists also negative eigenvalues for \code{CCA}
      and \code{pCCA}.}

   \item{\code{u}}{or row scores only give real axes with positive
      eigenvalues. The imaginary scores (if any) are in
      \code{imaginary.u}. The number of columns of real scores
      (positive eigenvalues) is given in item \code{poseig}. There is
      no \code{imaginary.u.eig}.}
      
  }
}


\note{
  In old versions of \pkg{vegan} the object also included scores
  scaled by eigenvalues (\code{u.eig}, \code{v.eig} and \code{wa.eig}),
  but these were removed in \pkg{vegan} 2.2-0. The scores are scaled
  when they are accessed with \code{\link{scores}} function. It is
  advisable to always use \code{\link{scores}} in accessing the
  results instead of directly accessing the elements of the the
  \code{cca} object.
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
