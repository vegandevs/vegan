\name{goodness.cca}
\alias{goodness}
\alias{goodness.rda}
\alias{goodness.cca}
\alias{inertcomp}
\alias{spenvcor}
\alias{intersetcor}
\alias{vif.cca}
\alias{alias.cca}

\title{Diagnostic Tools for [Constrained] Ordination (CCA,
  RDA, DCA, CA, PCA) }
\description{
  Functions \code{goodness} and \code{inertcomp} can
  be used to assess the goodness of fit for individual sites or
  species. Function \code{vif.cca} and \code{alias.cca} can be used to
  analyse linear dependencies among constraints and conditions. In
  addition, there are some other diagnostic tools (see 'Details').
  
}
\usage{
\method{goodness}{cca}(object, display = c("species", "sites"), choices,
    model = c("CCA", "CA"), statistic = c("explained", "distance"),
    summarize = FALSE, ...)
inertcomp(object, display = c("species", "sites"),
    statistic = c("explained", "distance"), proportional = FALSE)
spenvcor(object)
intersetcor(object)
vif.cca(object)
\method{alias}{cca}(object, names.only = FALSE, ...)
}

\arguments{
  \item{object}{A result object from \code{\link{cca}},
    \code{\link{rda}}, \code{\link{capscale}} or \code{\link{decorana}}. }
  \item{display}{Display \code{"species"} or \code{"sites"}. }
  \item{choices}{Axes shown. Default is to show all axes of the \code{"model"}. }
  \item{model}{Show constrained (\code{"CCA"}) or unconstrained
    (\code{"CA"}) results. }
  \item{statistic}{Statistic used: \code{"explained"} gives the cumulative
  percentage accounted for, \code{"distance"} shows the residual
  distances. Distances are not available for sites in constrained or
  partial analyses. }
  \item{summarize}{Show only the accumulated total. }
  \item{proportional}{Give the inertia components as proportional for
    the corresponding total.}
  \item{names.only}{Return only names of aliased variable(s) instead of
    defining equations.}
  \item{\dots}{Other parameters to the functions.}
}
\details{
  Function \code{goodness} gives the diagnostic statistics for species
  or sites. The alternative statistics are the cumulative proportion of
  inertia accounted for by the axes, and the residual distance left
  unaccounted for.  The conditional (``partialled out'') constraints are
  always regarded as explained and included in the statistics.

  Function \code{inertcomp} decomposes the inertia into partial,
  constrained and unconstrained components for each site or
  species. Instead of inertia, the function can give the total
  dispersion or distances from the centroid for each component.

  Function \code{spenvcor} finds the so-called \dQuote{species --
    environment correlation} or (weighted) correlation of 
  weighted average scores and linear combination scores.  This is a bad
  measure of goodness of ordination, because it is sensitive to extreme
  scores (like correlations are), and very sensitive to overfitting or
  using too many constraints. Better models often have poorer
  correlations. Function \code{\link{ordispider}} can show the same
  graphically.

  Function \code{intersetcor} finds the so-called \dQuote{interset
    correlation} or (weighted) correlation of weighted averages scores
  and constraints.  The defined contrasts are used for factor
  variables.  This is a bad measure since it is a correlation.  Further,
  it focuses on correlations between single contrasts and single axes
  instead of looking at the multivariate relationship.  Fitted vectors
  (\code{\link{envfit}}) provide a better alternative.  Biplot scores
  (see \code{\link{scores.cca}}) are a multivariate alternative for
  (weighted) correlation between linear combination scores and
  constraints. 
  
  Function \code{vif.cca} gives the variance inflation factors for each
  constraint or contrast in factor constraints. In partial ordination,
  conditioning variables are analysed together with constraints. Variance
  inflation is a diagnostic tool to identify useless constraints. A
  common rule is that values over 10 indicate redundant
  constraints. If later constraints are complete linear combinations of
  conditions or previous constraints, they will be completely removed
  from the estimation, and no biplot scores or centroids are calculated
  for these aliased constraints. A note will be printed with default
  output if there are aliased constraints. Function \code{alias} will
  give the linear coefficients defining the aliased constraints, or 
  only their names with argument \code{names.only = TRUE}.
}

\value{
  The functions return matrices or vectors as is appropriate.
}
\references{
  Greenacre, M. J. (1984). Theory and applications of correspondence
  analysis. Academic Press, London.

  Gross, J. (2003). Variance inflation factors. \emph{R News} 3(1),
  13--15. 
}
\author{Jari Oksanen. The \code{vif.cca} relies heavily on the code by
  W. N. Venables. \code{alias.cca} is a simplified version of
  \code{\link{alias.lm}}.}

\note{

  It is a common practise to use \code{goodness} statistics to remove
  species from ordination plots, but this may not be a good idea, as the
  total inertia is not a meaningful concept in \code{cca}, in particular
  for rare species.

  Function \code{vif} is defined as generic in package \pkg{car}
  (\code{\link[car]{vif}}), but if you have not loaded that package
  you must specify the call as \code{vif.cca}.  Variance inflation
  factor is useful diagnostic tool for  detecting nearly collinear
  constraints, but these are not a problem with algorithm used in this
  package to fit a constrained ordination.
}

\seealso{\code{\link{cca}}, \code{\link{rda}}, \code{\link{capscale}},
  \code{\link{decorana}}, \code{\link[car]{vif}}. }
\examples{
data(dune)
data(dune.env)
mod <- cca(dune ~ A1 + Management + Condition(Moisture), data=dune.env)
goodness(mod)
goodness(mod, summ = TRUE)
# Inertia components
inertcomp(mod, prop = TRUE)
inertcomp(mod, stat="d")
# vif.cca 
vif.cca(mod)
# Aliased constraints
mod <- cca(dune ~ ., dune.env)
mod
vif.cca(mod)
alias(mod)
with(dune.env, table(Management, Manure))
# The standard correlations (not recommended)
spenvcor(mod)
intersetcor(mod)
}
\keyword{ multivariate }

