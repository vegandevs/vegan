\name{ordicontrasts}
\alias{ordicontrasts}
\alias{ordicontrasts.default}
\alias{ordicontrasts.cca}
\alias{print.vegan_contrasts}
\title{Construct Reference Profiles for Ordination Contrasts}
\description{
Construct adjusted factor-level profiles from which ordination contrasts
are formed, optionally within groups. The displayed rows are labelled \code{b1},
\code{b2}, and so on, and can be used directly in \code{ordicomparisons}
expressions.
}
\usage{
ordicontrasts(object, ...)

\method{ordicontrasts}{cca}(object, variable, newdata = NULL,
    by = NULL, wts = NULL, ...)

\method{print}{vegan_contrasts}(x, ...)
}
\arguments{
\item{object}{A constrained ordination object fitted by the formula
interface, such as a result from \code{rda}, \code{cca}, or \code{dbrda}.}
\item{variable}{The name of one focal factor in the fitted model.}
\item{newdata}{The reference population used to construct profiles.
The default uses the fitted analysis data. Character value \code{"balanced"}
constructs a balanced reference grid, \code{"mean"} constructs a mean-or-mode
reference profile, and a data frame supplies an explicit population or
setting.}
\item{by}{Optional model variables defining groups within which profiles
are averaged.}
\item{wts}{Optional non-negative averaging weights, or the name of a weight
column in \code{newdata}.}
\item{x}{A \code{"vegan_contrasts"} object.}
\item{\dots}{Arguments passed to methods or to \code{print.data.frame}.}
}
\details{
Use \code{\link{ordicomparisons}} as the main interface for specifying and
testing comparisons. Use \code{ordicontrasts} to inspect the adjusted
profiles or refer to their rows in custom expressions. Individual rows
are quantities to be contrasted; they are not themselves necessarily
zero-sum contrasts.

A comparison profile represents one level of the focal factor with the other
variables averaged or held fixed. The printed table identifies what is being
compared but does not display a single response estimate for an individual
row. \code{ordicomparisons} tests differences among
the retained profile rows using constrained inertia and a permutation
pseudo-\eqn{F} statistic.

With \code{newdata = NULL}, every fitted analysis row is copied once for each
focal-factor level. Model-matrix rows are then averaged over the observed
values of the remaining variables. The printed \sQuote{Averaged over} field
states which variables vary in that reference population. If \code{by} is
used, averaging is performed separately within each displayed group. A
variable shown in a \code{by} column therefore does not also appear under
\sQuote{Averaged over}.

With \code{newdata = "balanced"}, non-focal factors form an equally weighted
grid and numeric variables are fixed at their means. With
\code{newdata = "mean"}, numeric variables are fixed at their arithmetic
means and other variables at their modes. This gives one representative
profile, or one profile within each \code{by} group. Tied factor modes are
resolved using fitted factor-level order; other ties use order of first
appearance. These comparisons are conditional on the resulting typical
profile rather than averaged over the observed sample or over factor levels.

With an explicit data frame, variables that have one supplied value are printed under
\sQuote{Fixed at}; variables that vary are printed under
\sQuote{Averaged over}. \code{N.reference} is the number of rows used to
construct a profile, not a sample size for a separate fitted model.

The stable \code{Profile} identifiers make otherwise cumbersome questions
easy to state. For example, \code{"b2 - b1 = 0"} compares two displayed rows,
and \code{"(b6 - b4) - (b3 - b1) = 0"} compares two simple effects across
groups. Because the identifiers depend on the displayed order, inspect the
profile table before writing a cross-group expression.

The method reconstructs the complete fitted model matrix, including
interactions and transformations, with the original factor levels and
contrast coding. Fitted transformation settings, such as the coefficients
used by \code{poly()} and the centre and scale used by \code{scale()},
are retained when the model is fitted and reused for the reference data.
They are not recalculated from the reference population. Settings are
retained before subset selection and missing-data removal, matching the
original fit. Older fitted objects without this information must be refitted
before making comparisons involving transformed terms.
User-defined transformations that learn settings from the data need to
supply fitted prediction information, usually through a
\code{makepredictcall()} method, when used as model terms.

Attributes retain these averaged design rows, the fitted
model, reference data, weights, grouping, and extrapolation diagnostics.
Because the displayed rows and retained design rows form one object, do not
edit or subset a profile table before passing it to \code{ordicomparisons};
reconstruct the desired profiles instead.
}
\value{
A data frame of class \code{"vegan_contrasts"}. Columns contain the profile
identifier in \code{Profile}, focal-factor value, optional \code{by} values,
\code{N.reference}, and an \code{Extrapolated} flag. Attributes retain the
model-matrix representation in \code{profile.matrix}, together with the
construction metadata used by \code{ordicomparisons}. The profile matrix
has one row per displayed profile, labelled by its \code{Profile} identifier.
}
\examples{
data(dune)
data(dune.env)

mod <- rda(dune ~ Management + A1, data = dune.env)
sc <- ordicontrasts(mod, "Management")
sc
ordicomparisons(sc, comparison = "b2 - b1 = 0", permutations = 99)

## Compare management types at the mean observed value of A1
mean.sc <- ordicontrasts(mod, "Management", newdata = "mean")
mean.sc
ordicomparisons(mean.sc, comparison = "b2 - b1 = 0", permutations = 99)

## Keep the fitted polynomial basis when comparing at the mean of A1
pmod <- rda(dune ~ Management + poly(A1, 2), data = dune.env)
ordicomparisons(pmod, "Management", comparison = "HF - BF = 0",
               newdata = "mean", permutations = 99)

set.seed(42)
idat <- expand.grid(treatment = factor(c("control", "low", "high"),
                                      levels = c("control", "low", "high")),
                    habitat = factor(c("grass", "heath")), replicate = 1:3)
iY <- matrix(rpois(nrow(idat) * 4, 5), ncol = 4)
imod <- rda(iY ~ treatment * habitat, data = idat)
isc <- ordicontrasts(imod, "treatment", by = "habitat")
isc
ordicomparisons(isc, comparison = "pairwise", permutations = 99)

fixed <- ordicontrasts(imod, "treatment", newdata = data.frame(
    habitat = factor("heath", levels = levels(idat$habitat))))
fixed
}
\seealso{\code{\link{ordicomparisons}}, \code{\link{rda}},
\code{\link{cca}}, \code{\link{dbrda}}}
\keyword{models}
\keyword{multivariate}
