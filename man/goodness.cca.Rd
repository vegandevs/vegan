\encoding{UTF-8}
\name{goodness.cca}
\alias{goodness}
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
\method{goodness}{cca}(object, choices, display = c("species", "sites"),
    model = c("CCA", "CA"), summarize = FALSE, addprevious = FALSE, ...)
inertcomp(object, display = c("species", "sites"),
    unity = FALSE, proportional = FALSE)
spenvcor(object)
intersetcor(object)
vif.cca(object)
\method{alias}{cca}(object, names.only = FALSE, ...)
}

\arguments{
  \item{object}{A result object from \code{\link{cca}},
    \code{\link{rda}}, \code{\link{dbrda}} or \code{\link{capscale}}. }

  \item{display}{Display \code{"species"} or \code{"sites"}. Species
    are not available in \code{\link{dbrda}} and \code{\link{capscale}}. }

  \item{choices}{Axes shown. Default is to show all axes of the
    \code{"model"}. }

  \item{model}{Show constrained (\code{"CCA"}) or unconstrained
    (\code{"CA"}) results. }

  \item{summarize}{Show only the accumulated total.}

  \item{addprevious}{Add the variation explained by previous components
     when \code{statistic="explained"}. For \code{model = "CCA"} add
     conditioned (partialled out) variation, and for \code{model = "CA"}
     add both conditioned and constrained variation. This will give
     cumulative explanation with previous components.
   }

  \item{unity}{Scale inertia components to unit sum (sum of all items is
    1).}
  \item{proportional}{Give the inertia components as proportional for
    the corresponding total of the item (sum of each row is 1). This
    option takes precedence over \code{unity}.}
  \item{names.only}{Return only names of aliased variable(s) instead of
    defining equations.}
  \item{\dots}{Other parameters to the functions.}
}
\details{

  Function \code{goodness} gives cumulative proportion of inertia
  accounted by species up to chosen axes. The proportions can be
  assessed either by species or by sites depending on the argument
  \code{display}, but species are not available in distance-based
  \code{\link{dbrda}}. The function is not implemented for
  \code{\link{capscale}}.

  Function \code{inertcomp} decomposes the inertia into partial,
  constrained and unconstrained components for each site or species.
  Legendre & De \enc{Cáceres}{Caceres} (2012) called these inertia
  components as local contributions to beta-diversity (LCBD) and
  species contributions to beta-diversity (SCBD), and they give these
  as relative contributions summing up to unity (argument
  \code{unity = TRUE}). For this interpretation, appropriate dissimilarity
  measures should be used in \code{\link{dbrda}} or appropriate
  standardization in \code{\link{rda}} (Legendre & De
  \enc{Cáceres}{Caceres} 2012). The function is not implemented for
  \code{\link{capscale}}.

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

  Legendre, P. & De \enc{Cáceres}{Caceres}, M. (2012). Beta diversity as
  the variance of community data: dissimilarity coefficients and
  partitioning. \emph{Ecology Letters} 16, 951--963.
  \doi{10.1111/ele.12141}
}
\author{Jari Oksanen. The \code{vif.cca} relies heavily on the code by
  W. N. Venables. \code{alias.cca} is a simplified version of
  \code{\link{alias.lm}}.}

\seealso{\code{\link{cca}}, \code{\link{rda}}, \code{\link{dbrda}},
  \code{\link{capscale}}. }

\examples{
data(dune)
data(dune.env)
mod <- cca(dune ~ A1 + Management + Condition(Moisture), data=dune.env)
goodness(mod, addprevious = TRUE)
goodness(mod, addprevious = TRUE, summ = TRUE)
# Inertia components
inertcomp(mod, prop = TRUE)
inertcomp(mod)
# vif.cca
vif.cca(mod)
# Aliased constraints
mod <- cca(dune ~ ., dune.env)
mod
vif.cca(mod)
alias(mod)
with(dune.env, table(Management, Manure))
# The standard correlations (not recommended)
## IGNORE_RDIFF_BEGIN
spenvcor(mod)
intersetcor(mod)
## IGNORE_RDIFF_END
}
\keyword{ multivariate }

